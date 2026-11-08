import numpy as np

Sz = np.diag([1, 0, -1]).astype(complex)
Sx = (1/np.sqrt(2)) * np.array([[0,1,0],[1,0,1],[0,1,0]], dtype=complex)
Sy = (1/np.sqrt(2)) * np.array([[0,-1j,0],[1j,0,-1j],[0,1j,0]], dtype=complex)

def hamiltonian(D, gamma, Bz, Mz=0.0, Mx=0.0, My=0.0, Nx=0.0, Ny=0.0):
    """general ground state H inlcuding the strain terms with basis |1>,  |0>, |-1>. ALl paramters are in GHz, inlcuding gamma*Bz.
    """

    H0 = D * (Sz @ Sz) + gamma * Bz * Sz
    Heps = (Mz * (Sz @ Sz)
            + Mx * (Sy @ Sy - Sx @ Sx)
            + My * (Sx @ Sy + Sy @ Sx)
            + Nx * (Sx @ Sz + Sz @ Sx)
            + Ny * (Sy @ Sz + Sz @ Sy))

    return H0 + Heps
