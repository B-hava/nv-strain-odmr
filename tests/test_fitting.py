# tests/test_fitting.py
import numpy as np
from src.fitting import (
    fit_spectrum,
    load_preprocessed_data,
    make_tlv_strain_tensor,
)
from src.lineshape import ensemble_lineshape


def test_make_tlv_strain_tensor_symmetric():
    eps = make_tlv_strain_tensor(-100e-5, -25e-5, 8e-5)
    np.testing.assert_allclose(eps, eps.T)
    assert eps[0, 0] == -100e-5
    assert eps[2, 2] == -25e-5
    assert eps[0, 2] == 8e-5
    assert eps[1, 1] == 0.0


def test_fit_spectrum_synthetic_roundtrip():
    freqs = np.linspace(2.84, 2.9045, 130)
    true_eTT, true_eVV, true_eTV = -103.9e-5, -27.1e-5, 8.4e-5
    true_D = 2.86670
    eps_true = make_tlv_strain_tensor(true_eTT, true_eVV, true_eTV)
    y_synth, _ = ensemble_lineshape(
        freqs, eps_true, B_lab=(0.0, 0.0, 1.0), A_list=[4.6e-5] * 4, D_val=true_D
    )

    res_fixed = fit_spectrum(
        freqs, y_synth, B_lab=(0.0, 0.0, 1.0), D_init=true_D, fit_D=False, equal_A=True
    )
    assert res_fixed.success
    assert res_fixed.rmse < 1e-10
    assert np.isclose(res_fixed.eps_TT, true_eTT, atol=1e-7)
    assert np.isclose(res_fixed.eps_VV, true_eVV, atol=1e-7)
    assert np.isclose(abs(res_fixed.eps_TV), abs(true_eTV), atol=1e-7)

    res_free = fit_spectrum(
        freqs, y_synth, B_lab=(0.0, 0.0, 1.0), D_init=2.868, fit_D=True, equal_A=True
    )
    assert res_free.success
    assert res_free.rmse < 1e-9
    assert np.isclose(res_free.D, true_D, atol=1e-6)
    assert np.isclose(res_free.eps_TT, true_eTT, atol=1e-6)
    assert np.isclose(res_free.eps_VV, true_eVV, atol=1e-6)
    assert np.isclose(abs(res_free.eps_TV), abs(true_eTV), atol=1e-6)


def test_fit_zenodo_milestone_position25():
    freqs, spectra = load_preprocessed_data()
    assert spectra.shape[1] == 47
    y25 = spectra[:, 24]  # Position 25 (2.49 um center)

    res = fit_spectrum(
        freqs, y25, B_lab=(0.0, 0.0, 1.0), D_init=2.86670, fit_D=False, equal_A=True
    )
    assert res.success
    assert res.rmse < 0.003
    assert abs(res.eps_TT - (-103.9e-5)) < 3e-5
    assert abs(res.eps_VV - (-27.1e-5)) < 4e-5
    assert abs(abs(res.eps_TV) - 8.4e-5) < 2e-5
