# scripts/check_imbalance_N_scaling.py
import numpy as np
from src.hamiltonian import hamiltonian
from src.observables import imbalance
from src.constants import D, GAMMA

# Fixed, generic M background (nonzero, non-axis-aligned)
Mz, Mx, My = 0.0, 0.02, 0.015

def imbalance_at(Nx, Ny):
    H = hamiltonian(D, GAMMA, Bz=0, Mz=Mz, Mx=Mx, My=My, Nx=Nx, Ny=Ny)
    eigvals, eigvecs = np.linalg.eigh(H)
    return imbalance(eigvals, eigvecs)

I0 = imbalance_at(0.0, 0.0)  # baseline, N=0

Nx0, Ny0 = 0.01, 0.007  # base perturbation scale
print(f"I(N=0) = {I0:.10f}\n")
print(f"{'scale':>8} {'|I(N)-I(0)|':>16} {'ratio to prev':>14}")

prev_dev = None
for scale in [1, 0.5, 0.25, 0.125, 0.0625]:
    I_N = imbalance_at(scale * Nx0, scale * Ny0)
    dev = abs(I_N - I0)
    ratio = (prev_dev / dev) if prev_dev else float('nan')
    print(f"{scale:8.4f} {dev:16.3e} {ratio:14.3f}")
    prev_dev = dev