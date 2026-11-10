from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

folder = Path(__file__).resolve().parent
data = np.genfromtxt(folder / 'pp_powerlaw_spectrum.csv', delimiter=',', names=True)
edges = np.append(data['E_low_GeV'], data['E_high_GeV'][-1])
# Bin-averaged E*dN/dE equals counts/(N_events * delta ln E).
# This normalization integrates over d ln E to the proton multiplicity.
dlogE = np.log(data['E_high_GeV'] / data['E_low_GeV'])
N_events = np.sum(data['N_injected'])
injected = data['N_injected'] / (N_events * dlogE)
protons = data['N_protons'] / (N_events * dlogE)
particles = np.genfromtxt(folder / 'pp_powerlaw_protons.csv', delimiter=',', names=True)
events = np.genfromtxt(folder / 'pp_powerlaw_events.csv', delimiter=',', names=True)
# Keep one leading proton per event AFTER all enabled hadron decays.
# This includes neutron/hyperon decay daughters. It approximates the outgoing
# projectile branch; identical protons cannot be uniquely assigned to beam/target.
# Events without a proton contribute no entry, retaining per-injected-event normalization.
leading_energy = np.full(int(N_events), np.nan)
for particle in particles:
    event_id = int(particle['event_id'])
    E = particle['E_proton_GeV']
    if np.isnan(leading_energy[event_id]) or E > leading_energy[event_id]:
        leading_energy[event_id] = E
has_proton = np.isfinite(leading_energy)
assert np.array_equal(has_proton, events['N_protons'] > 0)
leading_counts, _ = np.histogram(leading_energy[has_proton], bins=edges)
assert np.sum(leading_counts) == np.sum(has_proton)
leading_spectrum = leading_counts / (N_events * dlogE)
fig, ax = plt.subplots(figsize=(8, 5))
ax.stairs(injected, edges, label='Injected protons')
ax.stairs(protons, edges, label='All final protons, including decay daughters')
ax.stairs(leading_spectrum, edges, label='Leading final proton per event (after decays)')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Laboratory total proton energy [GeV]')
ax.set_ylabel(r'$(1/N_{\rm events})\,dN/d\ln E$')
ax.legend(fontsize=9)
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(folder / 'pp_powerlaw_spectrum.png', dpi=180)
fig.savefig(folder / 'pp_powerlaw_spectrum.pdf')
print('Events:', int(N_events))
print('Mean final proton multiplicity:', np.sum(data['N_protons']) / N_events)
print('Events retaining a proton:', int(np.sum(has_proton)))
print('Leading spectrum integral over d ln E:', np.sum(leading_spectrum * dlogE))
