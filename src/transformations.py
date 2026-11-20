# src/transformations.py
import numpy as np

T1 = np.array([
    [-1/np.sqrt(6),  1/np.sqrt(2), 1/np.sqrt(3)],
    [-1/np.sqrt(6), -1/np.sqrt(2), 1/np.sqrt(3)],
    [ np.sqrt(2/3),  0,            1/np.sqrt(3)],
])

T2 = np.array([
    [-1/np.sqrt(6),  1/np.sqrt(2),  1/np.sqrt(3)],
    [ 1/np.sqrt(6),  1/np.sqrt(2), -1/np.sqrt(3)],
    [-np.sqrt(2/3),  0,            -1/np.sqrt(3)],
])

T3 = np.array([
    [ 1/np.sqrt(6), -1/np.sqrt(2), -1/np.sqrt(3)],
    [-1/np.sqrt(6), -1/np.sqrt(2),  1/np.sqrt(3)],
    [-np.sqrt(2/3),  0,            -1/np.sqrt(3)],
])

T4 = np.array([
    [1/np.sqrt(6), -1/np.sqrt(2), -1/np.sqrt(3)],
    [1/np.sqrt(6),  1/np.sqrt(2), -1/np.sqrt(3)],
    [np.sqrt(2/3),  0,             1/np.sqrt(3)],
])

ALL_T = [T1, T2, T3, T4]


def rotate_strain(eps_cryst, T):
    """eps^(n) = T^(n).T @ eps^(cryst) @ T^(n), Appendix B."""
    return T.T @ eps_cryst @ T


def rotate_field(B_cryst, T):
    """B_MW^(n) = T^(n).T @ B_MW^(cryst)."""
    return T.T @ B_cryst