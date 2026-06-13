import numpy as np
import matplotlib.pyplot as plt
from pr_solver import MultiComponentPR

fluid_library = {
    'H2':  {'Tc': 59.94,  'Pc': 188.1,  'omega': -0.216},
    'CH4': {'Tc': 343.00, 'Pc': 667.0,  'omega': 0.011}
}

eos = MultiComponentPR(fluid_library)
pressures = np.linspace(500, 5000, 50)
reservoir_temp_r = 130.0 + 459.67  

z_pure_h2 = []
z_pure_ch4 = []
z_mixture = []

bip_matrix = np.zeros((2, 2))
bip_matrix[0, 1] = 0.088
bip_matrix[1, 0] = 0.088

for P in pressures:
    z_pure_h2.append(eos.calculate_mixture_z({'H2': 1.0, 'CH4': 0.0}, P, reservoir_temp_r))
    z_pure_ch4.append(eos.calculate_mixture_z({'H2': 0.0, 'CH4': 1.0}, P, reservoir_temp_r))
    z_mixture.append(eos.calculate_mixture_z({'H2': 0.5, 'CH4': 0.5}, P, reservoir_temp_r, bip_matrix))

plt.figure(figsize=(10, 6))
plt.plot(pressures, z_pure_h2, 'r-', linewidth=2, label='Pure Hydrogen (H2) - Resists Compression')
plt.plot(pressures, z_mixture, 'g--', linewidth=2, label='50% H2 / 50% CH4 Mixture')
plt.plot(pressures, z_pure_ch4, 'b-', linewidth=2, label='Pure Methane (CH4) - Highly Compressible')
plt.axhline(y=1.0, color='gray', linestyle=':', label='Ideal Gas Reference (Z = 1.0)')

plt.title('Gas Deviation Factor (Z) vs. Reservoir Pressure at 130°F', fontsize=14, fontweight='bold')
plt.xlabel('Reservoir Pressure (psi)', fontsize=12)
plt.ylabel('Compressibility Factor (Z)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=11, loc='upper left')

plt.savefig('z_factor_chart.png', dpi=300, bbox_inches='tight')
print("Chart successfully generated and saved as 'z_factor_chart.png'!")
