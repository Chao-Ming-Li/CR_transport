#pragma once
#include "CR_diffusion.hpp"  
#include <algorithm>
#include <cmath>
#include <span>
#include <stdexcept>


struct CNWorkspace
{
    // ----------------------------
    // Precomputed coefficients
    // ----------------------------
    Field2D temp;
    Array1D alpha_R;
    Array1D alpha_Z;

    // ----------------------------
    // Thomas matrix
    // ----------------------------
    Array1D Ra, Rb, Rb1, Rc, Rd, Rx, Rcp, Rdp;
    Array1D Za, Zb, Zb1, Zc, Zd, Zx, Zcp, Zdp;
 

    CNWorkspace()
        : alpha_R(NR), alpha_Z(NZ),
          Ra(NR), Rb(NR), Rb1(NR), Rc(NR), Rd(NR), Rx(NR),
          Za(NZ), Zb(NZ), Zb1(NZ), Zc(NZ), Zd(NZ), Zx(NZ),
          Rcp(NR), Rdp(NR), Zcp(NZ), Zdp(NZ)
    {
    }

    //----------------------------------------------------
    // Thomas solver
    //----------------------------------------------------
    void solve_tridiagonal_R(int n)
    {
        auto RA  = std::span(Ra.arr ).subspan(1,n);
        auto RB  = std::span(Rb1.arr ).subspan(1,n);
        auto RC  = std::span(Rc.arr ).subspan(1,n);
        auto RD  = std::span(Rd.arr ).subspan(1,n);
        auto RX  = std::span(Rx.arr ).subspan(1,n);
        auto RCP = std::span(Rcp.arr).subspan(1,n);
        auto RDP = std::span(Rdp.arr).subspan(1,n);

        if(n==0) return;

        //----------------------------
        // forward sweep
        //----------------------------
        RCP[0]=RC[0]/RB[0];
        RDP[0]=RD[0]/RB[0];

        for(int i=1;i<n;i++)
        {
            double denom=RB[i]-RA[i]*RCP[i-1];

            if(std::abs(denom)<1e-14)
            {
                throw std::runtime_error(
                    "Zero pivot in Thomas solver.");
            }

            double inv=1.0/denom;

            RCP[i]=(i<n-1)?RC[i]*inv:0.0;
            RDP[i]=(RD[i]-RA[i]*RDP[i-1])*inv;
        }

        //----------------------------
        // backward substitution
        //----------------------------
        RX[n-1]=RDP[n-1];

        for(int i=n-2;i>=0;i--)
        {
            RX[i]=RDP[i]-RCP[i]*RX[i+1];
        }
    }

    void solve_tridiagonal_Z(int n)
    {
        auto ZA  = std::span(Za.arr ).subspan(1,n);
        auto ZB  = std::span(Zb1.arr ).subspan(1,n);
        auto ZC  = std::span(Zc.arr ).subspan(1,n);
        auto ZD  = std::span(Zd.arr ).subspan(1,n);
        auto ZX  = std::span(Zx.arr ).subspan(1,n);
        auto ZCP = std::span(Zcp.arr).subspan(1,n);
        auto ZDP = std::span(Zdp.arr).subspan(1,n);

        if(n==0) return;

        //----------------------------
        // forward sweep
        //----------------------------
        ZCP[0]=ZC[0]/ZB[0];
        ZDP[0]=ZD[0]/ZB[0];

        for(int i=1;i<n;i++)
        {
            double denom=ZB[i]-ZA[i]*ZCP[i-1];

            if(std::abs(denom)<1e-14)
            {
                throw std::runtime_error(
                    "Zero pivot in Thomas solver.");
            }

            double inv=1.0/denom;

            ZCP[i]=(i<n-1)?ZC[i]*inv:0.0;
            ZDP[i]=(ZD[i]-ZA[i]*ZDP[i-1])*inv;
        }

        //----------------------------
        // backward substitution
        //----------------------------
        ZX[n-1]=ZDP[n-1];

        for(int i=n-2;i>=0;i--)
        {
            ZX[i]=ZDP[i]-ZCP[i]*ZX[i+1];
        }
    }


    //----------------------------------------------------
    // Crank-Nicolson + ADI solver
    //----------------------------------------------------  
    void CN_scheme_2nd_order(Field2D& ndis, 
                const Field2D& src, 
                const Field2D& ndis_H, 
                Field2D& dT_tau) 
    {

        // ==========================================
        // Step 1: Implicit R, Explicit Z
        // ==========================================
        for (int k = 0; k < NZ - 1; k++) {     
            for (int j = 0; j < NR - 1; j++) {
                Rd(j) = - Za(k) * ndis(j, k - 1) + (2.0 - Zb(k) - dT_tau(j, k) / 4.0) * ndis(j, k) - Zc(k) * ndis(j, k + 1);
            }
            
            for (int j = 0; j < NR - 1; j++) {Rb1(j) = Rb(j) + dT_tau(j, k) / 4.0;}
            Rb1(-1) = 1.0;
            solve_tridiagonal_R(NR);
            
            Rx(NR - 1) = 0.0; // enforce absorbing boundary at Rmax
            
            for (int j = 0; j < NR; j++) { temp(j, k) = Rx(j) + src(j, k) / 2.0; }
        }
        apply_boundary_conditions(temp);

        // ==========================================
        // Step 2: Implicit Z, Explicit R
        // ==========================================
        for (int j = 0; j < NR - 1; j++){
            for (int k = 0; k < NZ - 1; k++){
                Zd(k) = - Ra(j) * temp(j - 1, k) + (2.0 - Rb(j) - dT_tau(j, k) / 4.0) * temp(j, k) - Rc(j) * temp(j + 1, k);
            }    

            for (int k = 0; k < NZ - 1; k++) { Zb1(k) = Zb(k) + dT_tau(j, k) / 4.0;}
            Zb1(-1) = 1.0;

            solve_tridiagonal_Z(NZ);

            Zx(NZ - 1) = 0.0; // enforce absorbing boundary at Zmax
                
            for (int k = 0; k < NZ; k++) { ndis(j, k) = Zx(k) + src(j, k) / 2.0; }
        
        }
        apply_boundary_conditions(ndis);
    }
};


void initialize_CN_workspace( CNWorkspace& CN, Field2D& dT_tau, double D, double dt);