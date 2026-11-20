# scripts/case3_shear_xy.py
import numpy as np
import matplotlib.pyplot as plt
from src.hamiltonian import hamiltonian, make_strain_tensor, strain_params
from src.observables import shift_splitting, imbalance
from src.constants import D, GAMMA, H_IJ

EXX_BG = 1e-4
EYY_BG = -1e-4

exy_range = np.linspace(-50e-5, 50e-5, 51) # match paper's Fig. 3(III) x-axis range
shifts, splittings, imbalances = [], [], []

for exy in exy_range:
    eps = make_strain_tensor(exx=EXX_BG, eyy=EYY_BG, exy=exy)
    Mz, Mx, My, Nx, Ny = strain_params(eps, **H_IJ)
    H = hamiltonian(D, GAMMA, Bz=0, Mz=Mz, Mx=Mx, My=My, Nx=Nx, Ny=Ny)
    eigvals, eigvecs = np.linalg.eigh(H)

    s, sp = shift_splitting(eigvals, D)
    I = imbalance(eigvals, eigvecs)
    shifts.append(s * 1000)
    splittings.append(sp * 1000)
    imbalances.append(I)

shifts, splittings, imbalances = map(np.array, (shifts, splittings, imbalances))

# Sanity checks
h15 = H_IJ["h15"]
Mx_bg = -0.25 * h15 * (EXX_BG - EYY_BG)  # constant background contribution
floor_splitting_mhz = 2 * abs(Mx_bg) * 1000  # expected splitting at exy=0, ignoring N-terms

mid = len(exy_range) // 2
print(f"max |shift| (should be ~0): {np.max(np.abs(shifts)):.6f} MHz")
print(f"splitting at exy=0 (expected floor ~{floor_splitting_mhz:.4f} MHz): {splittings[mid]:.4f} MHz")
print(f"splitting range: [{splittings.min():.4f}, {splittings.max():.4f}] MHz")
print(f"max |imbalance|: {np.max(np.abs(imbalances)):.6e}")

fig, axes = plt.subplots(1, 3, figsize=(12, 3))
axes[0].plot(exy_range * 1e5, shifts);      axes[0].set_title("shift Δν (MHz)")
axes[1].plot(exy_range * 1e5, splittings);  axes[1].set_title("splitting δν (MHz)")
axes[2].plot(exy_range * 1e5, imbalances);  axes[2].set_title("imbalance I")
plt.tight_layout()
plt.savefig("case3_shear_xy.png", dpi=150)