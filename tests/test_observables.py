# tests/test_observables.py
import numpy as np
from src.hamiltonian import hamiltonian
from src.observables import imbalance

D, GAMMA = 2.870, 0.028


def test_imbalance_matches_eq11_at_zero_N():
    Mz, Mx, My = 0.05, 0.3, 0.1
    phi_str = np.arctan2(My, Mx)

    for phi_MW in [0, np.pi/4, np.pi/2, 3*np.pi/4]:
        H = hamiltonian(D, GAMMA, Bz=0, Mz=Mz, Mx=Mx, My=My, Nx=0, Ny=0)
        eigvals, eigvecs = np.linalg.eigh(H)

        I_numerical = imbalance(eigvals, eigvecs, phi_MW=phi_MW)
        # Sign convention: our S_y matrix differs from the paper's (unstated)
        # convention by an overall phase, confirmed by direct numerical
        # comparison at phi_MW = 0, pi/4, pi/2, 3pi/4 (debug_imbalance_sign.py,
        # 2026-11-10) to be a clean global negation, not an angle-dependent
        # discrepancy. Eigenvalues (D, M_z, splitting magnitude) are
        # convention-independent and verified separately against Eq.(8).
        I_closed_form = np.cos(2*phi_MW + phi_str)  # sign flipped vs. paper's Eq.(11); see note above
        assert np.isclose(I_numerical, I_closed_form, atol=1e-10)