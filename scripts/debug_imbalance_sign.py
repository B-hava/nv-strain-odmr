import numpy as np
from src.hamiltonian import hamiltonian
from src.observables import imbalance

D, GAMMA = 2.870, 0.028
Mz, Mx, My = 0.05, 0.3, 0.1
phi_str = np.arctan2(My, Mx)

H = hamiltonian(D, GAMMA, Bz=0, Mz=Mz, Mx=Mx, My=My, Nx=0, Ny=0)
eigvals, eigvecs = np.linalg.eigh(H)

print(f"{'phi_MW':>10} {'yours':>12} {'-cos(...)':>12} {'+cos(...)':>12}")
for phi_MW in [0, np.pi/4, np.pi/2, 3*np.pi/4]:
    I_yours = imbalance(eigvals, eigvecs, phi_MW=phi_MW)
    closed_neg = -np.cos(2*phi_MW + phi_str)
    closed_pos =  np.cos(2*phi_MW + phi_str)
    print(f"{phi_MW:10.4f} {I_yours:12.6f} {closed_neg:12.6f} {closed_pos:12.6f}")