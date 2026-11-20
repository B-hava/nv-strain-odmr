# scripts/case2_shear_yz.py
import numpy as np
import matplotlib.pyplot as plt
from src.hamiltonian import hamiltonian, make_strain_tensor, strain_params
from src.observables import shift_splitting, imbalance
from src.constants import D, GAMMA, H_IJ

eyz_range = np.linspace(-50e-5, 50e-5, 50)  # match paper's Fig. 3(II)-(b/c/d) x-axis range
shifts, splittings, imbalances = [], [], []

for eyz in eyz_range:
    eps = make_strain_tensor(eyz=eyz)  # exx=eyy=ezz=exy=exz=0
    Mz, Mx, My, Nx, Ny = strain_params(eps, **H_IJ)
    H = hamiltonian(D, GAMMA, Bz=0, Mz=Mz, Mx=Mx, My=My, Nx=Nx, Ny=Ny)
    eigvals, eigvecs = np.linalg.eigh(H)

    s, sp = shift_splitting(eigvals, D)
    I = imbalance(eigvals, eigvecs)
    shifts.append(s * 1000)
    splittings.append(sp * 1000)
    imbalances.append(I)

shifts, splittings, imbalances = map(np.array, (shifts, splittings, imbalances))

# Sanity check against the paper's N=0 approximation: delta_nu ~ |h16 * eyz|
h16 = H_IJ["h16"]
approx_splitting_mhz = np.abs(h16 * eyz_range) * 1000  # GHz -> MHz
rel_dev = np.abs(splittings - approx_splitting_mhz) / np.maximum(approx_splitting_mhz, 1e-12)
print(f"max |shift| (should be ~0): {np.max(np.abs(shifts)):.6f} MHz")
print(f"max |imbalance| (paper predicts 0 under N=0 approx): {np.max(np.abs(imbalances)):.6e}")
print(f"max relative deviation of exact splitting from N=0 approx: {np.max(rel_dev):.4%}")

fig, axes = plt.subplots(1, 3, figsize=(12, 3))
axes[0].plot(eyz_range * 1e5, shifts);      axes[0].set_title("shift Δν (MHz)")
axes[1].plot(eyz_range * 1e5, splittings);  axes[1].set_title("splitting δν (MHz)")
axes[2].plot(eyz_range * 1e5, imbalances);  axes[2].set_title("imbalance I")
plt.tight_layout()
plt.savefig("case2_shear_yz.png", dpi=150)