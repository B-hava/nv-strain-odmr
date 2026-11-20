# src/observables.py
import numpy as np
from src.hamiltonian import Sx, Sy, Sz

def shift_splitting(eigvals, D):
    """eigvals: sorted array of 3 eigenvalues [E_gs, E_-, E_+] (ascending).
    Returns (shift, splitting) per paper's definitions, Eqs. (12)-(13)."""
    E_gs, E_minus, E_plus = eigvals
    nu_minus = E_minus - E_gs
    nu_plus  = E_plus  - E_gs
    shift = (nu_plus + nu_minus)/2 - D
    splitting = nu_plus - nu_minus
    return shift, splitting


def transition_amplitudes(eigvecs, B_vec):
    """Raw (unnormalized) squared transition matrix elements:
    alpha_minus, alpha_plus = |<gs|S.B|minus/plus>|^2, for a general
    B_vec = (Bx, By, Bz) expressed in the same frame as Sx,Sy,Sz."""
    S_dot_B = B_vec[0]*Sx + B_vec[1]*Sy + B_vec[2]*Sz
    gs, minus, plus = eigvecs[:, 0], eigvecs[:, 1], eigvecs[:, 2]
    alpha_minus = np.abs(np.vdot(gs, S_dot_B @ minus))**2
    alpha_plus = np.abs(np.vdot(gs, S_dot_B @ plus))**2
    return alpha_minus, alpha_plus


def imbalance(eigvals, eigvecs, phi_MW=np.pi/2):
    B_vec = (np.cos(phi_MW), np.sin(phi_MW), 0.0)
    alpha_minus, alpha_plus = transition_amplitudes(eigvecs, B_vec)
    return (alpha_minus - alpha_plus) / (alpha_minus + alpha_plus)