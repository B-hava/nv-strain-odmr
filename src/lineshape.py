# src/lineshape.py
import numpy as np
from src.observables import transition_amplitudes
from src.transformations import ALL_T, rotate_strain, rotate_field
from src.hamiltonian import hamiltonian, strain_params
from src.constants import D, GAMMA, H_IJ

GAMMA_LINEWIDTH = 2.5e-3  # GHz (2.5 MHz), fixed value from Section V


def lorentzian(nu, gamma):
    """L(nu) = gamma / (nu^2 + gamma^2) — Eq. 6's building block."""
    return gamma / (nu**2 + gamma**2)


def single_nv_lineshape(nu_array, eigvals, eigvecs, B_vec, A=1.0, gamma=GAMMA_LINEWIDTH):
    """Eq. 6: f(nu) = 1 - A[alpha_- L(nu-nu_-) + alpha_+ L(nu-nu_+)],
    with alpha_-, alpha_+ normalized to sum to 1 (their relative share
    of this orientation's total transition strength — note the overall
    |B_perp| magnitude cancels in this ratio, so only B_vec's direction
    matters here, not its magnitude)."""
    E_gs, E_minus, E_plus = eigvals
    nu_minus, nu_plus = E_minus - E_gs, E_plus - E_gs

    a_minus, a_plus = transition_amplitudes(eigvecs, B_vec)
    total = a_minus + a_plus
    a_minus, a_plus = a_minus / total, a_plus / total

    return 1 - A * (a_minus * lorentzian(nu_array - nu_minus, gamma)
                     + a_plus * lorentzian(nu_array - nu_plus, gamma))


def ensemble_lineshape(nu_array, eps_lab, B_lab, A_list, Bz=0.0, gamma=GAMMA_LINEWIDTH, D_val=D):
    """Corrected Eq. 6+7 (baseline restored once, outside the sum, per
    the Supplementary Information):
        F(nu) = 1 - sum_{n=1}^4 A_n [alpha_n- L(nu-nu_n-) + alpha_n+ L(nu-nu_n+)]

    eps_lab : 3x3 strain tensor in the crystal/lab (T,L,V) frame.
    B_lab   : 3-vector, MW field direction in the crystal/lab (T,L,V) frame.
    A_list  : length-4 sequence of per-orientation amplitudes [A_1..A_4],
              matching Table II.
    D_val   : zero-field splitting (GHz), defaults to constants.D.

    Returns (F, components) where components is a list of the four
    individual "1 - A_n[...]" curves (for reproducing Fig. 4's
    orange/red/black/blue lines), and F is the actual summed total.
    """
    if len(A_list) != 4:
        raise ValueError(f"A_list must have length 4, got {len(A_list)}")
    B_lab = np.asarray(B_lab, dtype=float)
    if np.linalg.norm(B_lab) == 0.0:
        raise ValueError("B_lab must be a non-zero vector")

    dip_sum = np.zeros_like(nu_array, dtype=float)
    components = []

    for T, A_n in zip(ALL_T, A_list):
        eps_n = rotate_strain(eps_lab, T)
        B_n = rotate_field(B_lab, T)
        Mz, Mx, My, Nx, Ny = strain_params(eps_n, **H_IJ)
        H = hamiltonian(D_val, GAMMA, Bz, Mz, Mx, My, Nx, Ny)
        eigvals, eigvecs = np.linalg.eigh(H)

        a_minus, a_plus = transition_amplitudes(eigvecs, B_n)
        total = a_minus + a_plus
        a_minus, a_plus = a_minus / total, a_plus / total

        E_gs, E_minus, E_plus = eigvals
        nu_minus, nu_plus = E_minus - E_gs, E_plus - E_gs

        dip_n = A_n * (a_minus * lorentzian(nu_array - nu_minus, gamma)
                        + a_plus * lorentzian(nu_array - nu_plus, gamma))
        dip_sum += dip_n
        components.append(1 - dip_n)  # matches single_nv_lineshape's convention

    F = 1 - dip_sum
    return F, components
