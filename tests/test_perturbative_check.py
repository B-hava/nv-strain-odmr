# tests/test_perturbative_check.py
import numpy as np
from src.hamiltonian import hamiltonian
from src.observables import shift_splitting
from src.constants import D, GAMMA


def perturbative_shift_splitting(D, Mz, Mx, My, Nx, Ny):
    """Leading-order (2nd-order) non-degenerate perturbation theory,
    treating Nx,Ny as a small perturbation on the exactly-solvable
    Nx=Ny=0 problem (Eq. 8). Hand-derived independently of
    src.hamiltonian's own algebra — see NOTES.md for the derivation.
    Valid only for |Mx,My,Nx,Ny| << D.
    """
    Mperp2 = Mx**2 + My**2
    Nperp2 = Nx**2 + Ny**2
    phi_M = np.arctan2(My, Mx)
    phi_N = np.arctan2(Ny, Nx)
    shift = Mz + 1.5 * Nperp2 / D
    splitting = 2*np.sqrt(Mperp2) + Nperp2 * np.cos(phi_M + 2*phi_N) / D
    return shift, splitting


def test_perturbative_matches_exact_small_N():
    # Deliberately generic / non-axis-aligned values, so the Mperp-Nperp
    # cross term in the splitting formula is actually exercised rather
    # than accidentally zeroed by a special phase relationship (as
    # happens for pure eyz shear, see Case 2).
    Mz, Mx, My, Nx, Ny = 0.01, 0.02, 0.015, 0.003, 0.002

    H = hamiltonian(D, GAMMA, Bz=0, Mz=Mz, Mx=Mx, My=My, Nx=Nx, Ny=Ny)
    eigvals = np.sort(np.linalg.eigvalsh(H))
    shift_exact, splitting_exact = shift_splitting(eigvals, D)
    shift_pert, splitting_pert = perturbative_shift_splitting(D, Mz, Mx, My, Nx, Ny)

    rel_shift_err = abs(shift_exact - shift_pert) / abs(shift_pert)
    rel_split_err = abs(splitting_exact - splitting_pert) / abs(splitting_pert)

    print(f"shift: exact={shift_exact:.8f} GHz, pert={shift_pert:.8f} GHz, rel_err={rel_shift_err:.4%}")
    print(f"splitting: exact={splitting_exact:.8f} GHz, pert={splitting_pert:.8f} GHz, rel_err={rel_split_err:.4%}")

    # Tolerance is provisional — PT accuracy depends on how small N/D
    # actually is here. Run once, read the printed errors, then
    # tighten/loosen this before trusting it as a permanent check.
    assert rel_shift_err < 0.001
    assert rel_split_err < 0.001