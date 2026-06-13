import numpy as np

class MultiComponentPR:
    """
    Advanced Peng-Robinson Equation of State Solver for Multi-Component Gas Mixtures.
    Specifically tailored for Underground Hydrogen Storage (UHS) screening.
    """
    R = 10.73146  # Gas constant: ft3 * psi / (lb-mol * R)

    def __init__(self, components_data):
        self.comp_data = components_data
        self.names = list(components_data.keys())
        self.num_comps = len(self.names)
        
        self.omega_a = 0.45724
        self.omega_b = 0.07780
        
        self.m = {}
        for name in self.names:
            w = components_data[name]['omega']
            if w <= 0.491:
                self.m[name] = 0.37464 + 1.54226 * w - 0.26992 * (w ** 2)
            else:
                self.m[name] = 0.379642 + 1.48503 * w - 0.164423 * (w ** 2)

    def calculate_mixture_z(self, mole_fractions, p_psi, t_r, bips=None):
        if isinstance(mole_fractions, dict):
            x = [mole_fractions[name] for name in self.names]
        else:
            x = mole_fractions
            
        x = np.array(x) / np.sum(x)
        if bips is None:
            bips = np.zeros((self.num_comps, self.num_comps))

        a_pure = np.zeros(self.num_comps)
        b_pure = np.zeros(self.num_comps)
        
        for i, name in enumerate(self.names):
            Tc = self.comp_data[name]['Tc']
            Pc = self.comp_data[name]['Pc']
            tr = t_r / Tc
            alpha = (1 + self.m[name] * (1 - np.sqrt(tr))) ** 2
            
            ac = self.omega_a * (self.R ** 2) * (Tc ** 2) / Pc
            a_pure[i] = ac * alpha
            b_pure[i] = self.omega_b * self.R * Tc / Pc

        b_mix = np.sum(x * b_pure)
        a_mix = 0.0
        for i in range(self.num_comps):
            for j in range(self.num_comps):
                a_ij = np.sqrt(a_pure[i] * a_pure[j]) * (1.0 - bips[i, j])
                a_mix += x[i] * x[j] * a_ij

        A = (a_mix * p_psi) / ((self.R ** 2) * (t_r ** 2))
        B = (b_mix * p_psi) / (self.R * t_r)
        
        c2 = -(1.0 - B)
        c1 = A - 3.0 * (B ** 2) - 2.0 * B
        c0 = -(A * B - (B ** 2) - (B ** 3))
        
        roots = np.roots([1.0, c2, c1, c0])
        real_roots = [r.real for r in roots if np.isreal(r) and r.real > 0]
        
        if not real_roots:
            raise ValueError("Thermodynamic state yield no valid real roots for Z.")
            
        return max(real_roots)

    def evaluate_reservoir_capacity(self, pv_cf, mole_fractions, p_psi, t_r, bips=None):
        z_mix = self.calculate_mixture_z(mole_fractions, p_psi, t_r, bips)
        bg = z_mix * (t_r / 519.67) * (14.7 / p_psi)
        total_gas_scf = pv_cf / bg
        return total_gas_scf, z_mix, bg

if __name__ == "__main__":
    fluid_library = {
        'H2':  {'Tc': 59.94,  'Pc': 188.1,  'omega': -0.216},
        'CH4': {'Tc': 343.00, 'Pc': 667.0,  'omega': 0.011},
        'CO2': {'Tc': 547.40, 'Pc': 1071.0, 'omega': 0.225}
    }
    eos = MultiComponentPR(fluid_library)
    net_pore_volume = 43560 * 2000 * 50 * 0.15 * (1 - 0.25)
    composition = {'H2': 0.80, 'CH4': 0.20, 'CO2': 0.00}
    P_reservoir = 3500.0 
    T_reservoir = 130.0 + 459.67
    
    bip_matrix = np.zeros((3, 3))
    bip_matrix[0, 1] = 0.088
    bip_matrix[1, 0] = 0.088
    
    gas_scf, z, bg = eos.evaluate_reservoir_capacity(net_pore_volume, composition, P_reservoir, T_reservoir, bip_matrix)
    print("--- GitHub Portfolio Component Validation ---")
    print(f"Mixture Composition: 80% H2 / 20% CH4")
    print(f"Gas Deviation Factor (Z-mix): {z:.4f}")
    print(f"Gas Formation Volume Factor (Bg): {bg:.6f} ft3/SCF")
    print(f"Total Stored Mixture Volume: {gas_scf / 1e9:.3f} BSCF")
