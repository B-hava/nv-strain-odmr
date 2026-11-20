# tests/test_transformations.py
import numpy as np
from src.transformations import ALL_T

def test_all_T_orthogonal():
    for i, T in enumerate(ALL_T, start=1):
        np.testing.assert_allclose(T.T @ T, np.eye(3), atol=1e-12,
            err_msg=f"T^({i}) is not orthogonal")

def test_all_T_determinant_plus_one():
    # proper rotations, not reflections
    for i, T in enumerate(ALL_T, start=1):
        det = np.linalg.det(T)
        assert np.isclose(det, 1.0, atol=1e-12), f"T^({i}) det={det}, expected +1"