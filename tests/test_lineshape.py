# tests/test_lineshape.py
import numpy as np
from src.constants import D, GAMMA, H_IJ
from src.hamiltonian import hamiltonian, make_strain_tensor, strain_params
from src.lineshape import (
    GAMMA_LINEWIDTH,
    ensemble_lineshape,
    lorentzian,
    single_nv_lineshape,
)
from src.observables import imbalance


def test_lorentzian_properties():
    gamma = 2.5e-3
    assert np.isclose(lorentzian(0.0, gamma), 1.0 / gamma)
    assert np.isclose(lorentzian(gamma, gamma), 0.5 / gamma)
    assert np.isclose(lorentzian(-gamma, gamma), 0.5 / gamma)


def test_single_nv_lineshape_baseline_and_asymmetry():
    Mz, Mx, My = 0.0, 0.01, 0.005
    H = hamiltonian(D, GAMMA, Bz=0.0, Mz=Mz, Mx=Mx, My=My)
    eigvals, eigvecs = np.linalg.eigh(H)
    E_gs, E_m, E_p = eigvals
    nu_m, nu_p = E_m - E_gs, E_p - E_gs

    far = np.array([2.50, 3.20])
    f_far = single_nv_lineshape(far, eigvals, eigvecs, B_vec=(1.0, 0.0, 0.0), A=1e-4)
    np.testing.assert_allclose(f_far, 1.0, atol=1e-5)

    phi_MW = np.pi / 3
    B_vec = (np.cos(phi_MW), np.sin(phi_MW), 0.0)
    I_val = imbalance(eigvals, eigvecs, phi_MW=phi_MW)
    f_res = single_nv_lineshape(np.array([nu_m, nu_p]), eigvals, eigvecs, B_vec=B_vec, A=1e-4)
    dip_m, dip_p = 1.0 - f_res[0], 1.0 - f_res[1]
    assert np.sign(dip_m - dip_p) == np.sign(I_val)


def test_ensemble_lineshape_baseline_is_one():
    nu = np.array([2.50, 3.20])
    eps = make_strain_tensor(exx=-50e-5, ezz=-20e-5, exz=10e-5)
    F, _ = ensemble_lineshape(nu, eps, B_lab=(0.0, 0.0, 1.0), A_list=[5e-5] * 4)
    np.testing.assert_allclose(F, 1.0, atol=1e-5)


def test_ensemble_lineshape_zero_strain():
    nu = np.linspace(2.85, 2.89, 81)
    eps = make_strain_tensor()
    A_list = [2e-5, 3e-5, 4e-5, 5e-5]
    F, components = ensemble_lineshape(nu, eps, B_lab=(0.0, 0.0, 1.0), A_list=A_list)
    expected = 1.0 - sum(A_list) * lorentzian(nu - D, GAMMA_LINEWIDTH)
    np.testing.assert_allclose(F, expected, atol=1e-12)
    assert len(components) == 4


def test_ensemble_lineshape_volumetric_strain():
    ev = -20e-5
    eps = make_strain_tensor(exx=ev, eyy=ev, ezz=ev)
    Mz, _, _, _, _ = strain_params(eps, **H_IJ)
    nu = np.linspace(2.85, 2.89, 81)
    A_list = [4e-5, 4e-5, 4e-5, 4e-5]
    F, _ = ensemble_lineshape(nu, eps, B_lab=(0.0, 0.0, 1.0), A_list=A_list)
    expected = 1.0 - sum(A_list) * lorentzian(nu - (D + Mz), GAMMA_LINEWIDTH)
    np.testing.assert_allclose(F, expected, atol=1e-12)


def test_ensemble_lineshape_components_sum():
    nu = np.linspace(2.84, 2.90, 61)
    eps = make_strain_tensor(exx=-100e-5, ezz=-25e-5, exz=8e-5)
    A_list = [3e-5, 5e-5, 4e-5, 6e-5]
    F, components = ensemble_lineshape(nu, eps, B_lab=(0.2, 0.1, 1.0), A_list=A_list)
    reconstructed = 1.0 - sum(1.0 - c for c in components)
    np.testing.assert_allclose(F, reconstructed, atol=1e-12)
