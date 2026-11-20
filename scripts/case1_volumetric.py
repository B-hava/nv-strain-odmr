# scripts/case1_volumetric.py
import numpy as np
import matplotlib.pyplot as plt
from src.hamiltonian import hamiltonian, make_strain_tensor, strain_params
from src.observables import shift_splitting, imbalance
from src.constants import D, GAMMA, H_IJ

ev_range = np.linspace(-50e-5, 50e-5, 50)  # match paper's Fig. 3(I-b) x-axis range
shifts, splittings, imbalances = [], [], []

for ev in ev_range:
    eps = make_strain_tensor(exx=ev, eyy=ev, ezz=ev)
    Mz, Mx, My, Nx, Ny = strain_params(eps, **H_IJ)
    H = hamiltonian(D, GAMMA, Bz=0, Mz=Mz, Mx=Mx, My=My, Nx=Nx, Ny=Ny)
    eigvals, eigvecs = np.linalg.eigh(H)  # already ascending

    s, sp = shift_splitting(eigvals, D)
    I = imbalance(eigvals, eigvecs)
    shifts.append(s*1000)       # GHz -> MHz, matches their axis units
    splittings.append(sp*1000)
    imbalances.append(I)

fig, axes = plt.subplots(1, 3, figsize=(12,3))
axes[0].plot(ev_range*1e5, shifts);     axes[0].set_title("shift Δν (MHz)")
axes[1].plot(ev_range*1e5, splittings); axes[1].set_title("splitting δν (MHz)")
axes[2].plot(ev_range*1e5, imbalances); axes[2].set_title("imbalance I")
plt.tight_layout()
plt.savefig("case1_volumetric.png", dpi=150)
