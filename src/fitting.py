# src/fitting.py
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, nnls

from src.constants import GAMMA, H_IJ
from src.hamiltonian import hamiltonian, strain_params
from src.lineshape import GAMMA_LINEWIDTH, ensemble_lineshape, lorentzian
from src.observables import transition_amplitudes
from src.transformations import ALL_T, rotate_field, rotate_strain

DEFAULT_PREPROCESSED_PATH = (
    Path(__file__).resolve().parent.parent / "zenodo-data-set" / "Preprocessed_data.txt"
)


@dataclass
class FitResult:
    eps_TT: float
    eps_VV: float
    eps_TV: float
    D: float
    A_list: np.ndarray
    rmse: float
    y_fit: np.ndarray
    components: list
    success: bool
    message: str


def load_preprocessed_data(filepath=DEFAULT_PREPROCESSED_PATH):
    """Load Preprocessed_data.txt.

    Returns
    -------
    freqs : np.ndarray, shape (n_freq,)
        Microwave frequency grid in GHz.
    spectra : np.ndarray, shape (n_freq, n_positions)
        Normalized ODMR spectra across all spatial positions.
    """
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, 1:]


def make_tlv_strain_tensor(eps_TT, eps_VV, eps_TV):
    """Build 3x3 symmetric 2D plane-strain tensor in the (T, L, V) frame
    (with eps_LL = eps_TL = eps_LV = 0)."""
    return np.array(
        [
            [eps_TT, 0.0, eps_TV],
            [0.0, 0.0, 0.0],
            [eps_TV, 0.0, eps_VV],
        ],
        dtype=float,
    )


def _orientation_dip_basis(freqs, eps_lab, B_lab, D_val, Bz=0.0, gamma=GAMMA_LINEWIDTH):
    """Return (n_freq, 4) matrix whose columns are the unit-amplitude dip
    profiles alpha_- L(nu-nu_-) + alpha_+ L(nu-nu_+) for orientations 1..4."""
    basis = []
    for T in ALL_T:
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
        basis.append(
            a_minus * lorentzian(freqs - nu_minus, gamma)
            + a_plus * lorentzian(freqs - nu_plus, gamma)
        )
    return np.column_stack(basis)


def _solve_amplitudes(X, y_dip, equal_A):
    if equal_A:
        x_sum = np.sum(X, axis=1)
        denom = float(np.dot(x_sum, x_sum))
        a = max(0.0, float(np.dot(x_sum, y_dip)) / denom) if denom > 0 else 0.0
        return np.full(4, a, dtype=float)
    A_fit, _ = nnls(X, y_dip)
    return A_fit


def fit_spectrum(
    freqs,
    y_obs,
    B_lab=(0.0, 0.0, 1.0),
    D_init=2.86670,
    fit_D=True,
    equal_A=True,
    Bz=0.0,
    gamma=GAMMA_LINEWIDTH,
    init_strain=(-100e-5, -30e-5, 10e-5),
):
    """Fit the 4-orientation ensemble ODMR lineshape to an observed spectrum
    using variable-projection nonlinear least squares.

    Parameters
    ----------
    freqs : array-like
        Frequency array in GHz.
    y_obs : array-like
        Normalized ODMR spectrum (baseline ~1.0).
    B_lab : 3-vector
        MW magnetic field direction in the (T, L, V) frame.
    D_init : float
        Initial (or fixed, if fit_D=False) zero-field splitting D in GHz.
    fit_D : bool
        Whether to optimize D alongside (eps_TT, eps_VV, eps_TV).
    equal_A : bool
        If True, fit a single shared amplitude A across all 4 orientations;
        if False, fit per-orientation amplitudes A_1..A_4 >= 0 via NNLS.
    Bz : float
        Static axial magnetic field in mT (default 0.0).
    gamma : float
        Lorentzian HWHM linewidth in GHz.
    init_strain : 3-tuple
        Initial guess (eps_TT, eps_VV, eps_TV) in dimensionless strain.
    """
    freqs = np.asarray(freqs, dtype=float)
    y_obs = np.asarray(y_obs, dtype=float)
    B_lab = np.asarray(B_lab, dtype=float)
    if np.linalg.norm(B_lab) == 0.0:
        raise ValueError("B_lab must be a non-zero vector")
    B_lab = B_lab / np.linalg.norm(B_lab)

    y_dip = 1.0 - y_obs
    s0 = np.asarray(init_strain, dtype=float) * 1e5  # scale to 10^-5 units

    if fit_D:
        p0 = np.array([s0[0], s0[1], s0[2], D_init], dtype=float)
        bounds = ([-500.0, -500.0, -500.0, 2.84], [500.0, 500.0, 500.0, 2.90])
    else:
        p0 = s0.copy()
        bounds = ([-500.0, -500.0, -500.0], [500.0, 500.0, 500.0])

    def residual(p):
        eps_TT, eps_VV, eps_TV = p[0] * 1e-5, p[1] * 1e-5, p[2] * 1e-5
        D_val = p[3] if fit_D else D_init
        eps_lab = make_tlv_strain_tensor(eps_TT, eps_VV, eps_TV)
        X = _orientation_dip_basis(freqs, eps_lab, B_lab, D_val, Bz=Bz, gamma=gamma)
        A_fit = _solve_amplitudes(X, y_dip, equal_A=equal_A)
        return X @ A_fit - y_dip

    res = least_squares(residual, p0, bounds=bounds)

    eps_TT, eps_VV, eps_TV = res.x[0] * 1e-5, res.x[1] * 1e-5, res.x[2] * 1e-5
    D_fit = float(res.x[3]) if fit_D else float(D_init)
    eps_lab = make_tlv_strain_tensor(eps_TT, eps_VV, eps_TV)
    X = _orientation_dip_basis(freqs, eps_lab, B_lab, D_fit, Bz=Bz, gamma=gamma)
    A_list = _solve_amplitudes(X, y_dip, equal_A=equal_A)
    y_fit, components = ensemble_lineshape(
        freqs, eps_lab, B_lab, A_list, Bz=Bz, gamma=gamma, D_val=D_fit
    )
    rmse = float(np.sqrt(np.mean((y_fit - y_obs) ** 2)))

    return FitResult(
        eps_TT=float(eps_TT),
        eps_VV=float(eps_VV),
        eps_TV=float(eps_TV),
        D=D_fit,
        A_list=A_list,
        rmse=rmse,
        y_fit=y_fit,
        components=components,
        success=bool(res.success),
        message=str(res.message),
    )
