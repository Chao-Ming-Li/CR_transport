#include <iostream>
#include <vector>

int main() {
    std::vector<double_t> dZ;
    std::vector<double_t> Z;
    double_t Z_temp = 0.0;
    double_t dZ_temp = 0.02;
    while (Z_temp <= 200.0) {
        if (Z_temp < 0.2){ dZ_temp = 0.02; }
        else if (Z_temp < 20.0 && Z_temp > 0.2){
            double x = (Z_temp - 0.2) / (20.0 - 0.2); // Normalize Z_temp to [0, 1]
            dZ_temp = 0.02 + (1.0 - 0.02) * (6. * pow(x, 5.) - 15. * pow(x, 4.) + 10. * pow(x, 3.)); 
        }
        else { dZ_temp = 1.0; }
        
        dZ.push_back(dZ_temp);
        Z_temp += dZ_temp;
        Z.push_back(Z_temp);
    }

    // Output the results
    std::cout << "Z values:\n";
    for (const auto& z : Z) {
        std::cout << z << " ";
    }
    std::cout << "\n\nCorresponding dZ values:\n";
    for (const auto& dz : dZ) {
        std::cout << dz << " ";
    }
    std::cout << std::endl;

    return 0;
}