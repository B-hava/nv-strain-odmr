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


# defines the strain tensor e_ij
def make_strain_tensor(exx=0.0,  eyy=0.0, ezz=0.0, exy=0.0, exz=0.0, eyz=0.0):
    # symmetric 3x3 strain tensor in the NV-centered (x,y,z) frame
    return np.array([
        [exx, exy, exz],
        [exy, eyy, eyz],
        [exz, eyz, ezz]
    ])


# takes in the 6 independent strain tensor components and the hij table and returns (M_z,M_x,M_y,N_x,N_y) via Eqs. (3a–3e)
def strain_params(eps_tensor, h41, h43, h15, h16, h25, h26):
    exx, eyy, ezz = eps_tensor[0,0], eps_tensor[1,1], eps_tensor[2,2]
    exy, exz, eyz = eps_tensor[0,1], eps_tensor[0,2], eps_tensor[1,2]
    
    Mz = h41*(exx + eyy) + h43*ezz
    Mx = 0.5*(h16*exz - 0.5*h15*(exx - eyy))
    My = 0.5*(h16*eyz + h15*exy)
    Nx = 0.5*(h26*exz - 0.5*h25*(exx - eyy))
    Ny = 0.5*(h26*eyz + h25*exy)
    return Mz, Mx, My, Nx, Ny