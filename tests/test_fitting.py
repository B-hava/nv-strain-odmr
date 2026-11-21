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


def test_fit_zenodo_mapped_positions():
    freqs, spectra = load_preprocessed_data()
    assert spectra.shape[1] == 47

    # Position 24 (col index 23) = Table II position no. 23 (0.0 um, Fig. 4 center)
    res24 = fit_spectrum(
        freqs, spectra[:, 23], B_lab=(0.0, 0.0, 1.0), D_init=2.8670, fit_D=True, equal_A=True
    )
    assert res24.success
    assert res24.rmse < 0.003
    assert abs(res24.eps_TT - (-103.94e-5)) < 4e-5
    assert abs(res24.eps_VV - (-27.13e-5)) < 3e-5
    assert abs(abs(res24.eps_TV) - 8.4e-5) < 3e-5
    assert abs(res24.D - 2.867) < 0.001

    # Position 25 (col index 24) = Table II position no. 24 (0.83 um)
    res25 = fit_spectrum(
        freqs, spectra[:, 24], B_lab=(0.0, 0.0, 1.0), D_init=2.8670, fit_D=True, equal_A=True
    )
    assert res25.success
    assert res25.rmse < 0.003
    assert abs(res25.eps_TT - (-98.83e-5)) < 2e-5
    assert abs(res25.eps_VV - (-17.28e-5)) < 2e-5
    assert abs(res25.D - 2.8676) < 0.0005

    # Position 27 (col index 26) = Table II position no. 26 (2.49 um)
    res27 = fit_spectrum(
        freqs, spectra[:, 26], B_lab=(0.0, 0.0, 1.0), D_init=2.8670, fit_D=True, equal_A=True
    )
    assert res27.success
    assert res27.rmse < 0.003
    assert abs(res27.eps_TT - (-100.98e-5)) < 4e-5
    assert abs(res27.eps_VV - (-11.38e-5)) < 3e-5
    assert abs(abs(res27.eps_TV) - 15.09e-5) < 3e-5
    assert abs(res27.D - 2.8686) < 0.0005
