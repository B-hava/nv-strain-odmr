import numpy as np
import pytest
from src.hamiltonian import hamiltonian

D = 2.870 #GHz
gamma = 0.028 #GHz


#testing if hermitian 
def test_hermititian():
    H = hamiltonian(D, gamma, Bz=1.0, Mx=0.3, My=0.1, Nx=0.2, Ny=0.05, Mz=0.1)
    assert np.allclose(H, H.conj().T)

def test_zero_strain_reduces_to_H0():
    H = hamiltonian(D, gamma, Bz=2.0)
    eigs = np.sort(np.linalg.eigvalsh(H))
    expected = np.sort([D - gamma*2.0, 0.0, D + gamma*2.0])
    assert np.allclose(eigs, expected)

def test_matches_closed_form_eq8():
    # N_x = N_y = 0, B_z = 0: compare against the papers obtained eigenvalues 
    Mz, Mx, My = 0.05, 0.3, 0.1
    H = hamiltonian(D, gamma, Bz=0, Mz=Mz, Mx=Mx, My=My)
    eigs = np.sort(np.linalg.eigvalsh(H))
    E_gs = 0.0
    E_minus = D + Mz - np.sqrt(Mx**2 + My**2)
    E_plus  = D + Mz + np.sqrt(Mx**2 + My**2)
    expected = np.sort([E_gs, E_minus, E_plus])
    assert np.allclose(eigs, expected)