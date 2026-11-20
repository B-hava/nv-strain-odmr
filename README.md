# nv-strain-odmr

Simulation code for strain- and magnetic-field-dependent ODMR (optically detected magnetic resonance) spectra of nitrogen-vacancy (NV) centers in diamond, built to reproduce the model and results of Alam et al.

## Final goal

Replicate the paper's full pipeline end to end:

1. **Single-NV spin Hamiltonian** — build the 3×3 ground-state spin-1 Hamiltonian including zero-field splitting (D), Zeeman coupling (γ), and the strain-induced M_z, M_x, M_y, N_x, N_y terms, and diagonalize it exactly.
2. **Four-orientation ensemble model** — rotate the crystal-frame strain and microwave (MW) field into each of the four NV symmetry axes using the paper's Appendix B transformation matrices, compute each orientation's resonance, and sum the four single-NV Lorentzian lineshapes into one ensemble ODMR spectrum.
3. **Data fitting** — load the paper's published Zenodo dataset and do a nonlinear least-squares fit of the ensemble model to real ODMR data, reproducing their reported strain-tensor and D values at the 2.49 µm sample position as the final replication milestone.

## Quickstart

```bash
source venv/bin/activate
pytest -v                        # 16 passed
python -m scripts.fit_milestone  # fits Position 25 (2.49 um) -> milestone_fit_pos25.png
```

## Current state (as of 20/11/26)

### Repo layout
```text
.
├── NOTES.md
├── README.md
├── milestone_fit_pos25.png
├── requirements.txt
├── scripts/
│   ├── case1_volumetric.py
│   ├── case2_shear_yz.py
│   ├── case3_shear_xy.py
│   ├── check_imbalance_N_scaling.py
│   ├── debug_imbalance_sign.py
│   └── fit_milestone.py
├── src/
│   ├── constants.py
│   ├── fitting.py
│   ├── hamiltonian.py
│   ├── lineshape.py
│   ├── observables.py
│   └── transformations.py
├── tests/
│   ├── test_fitting.py
│   ├── test_hamiltonian.py
│   ├── test_lineshape.py
│   ├── test_observables.py
│   ├── test_perturbative_check.py
│   └── test_transformations.py
└── zenodo-data-set/
    ├── Preprocessed_data.txt
    └── Raw_data.txt
```

### Validated (Phase 0 — Single-NV Hamiltonian & observables)
- Hamiltonian matches the hand-derived form (D, γ, M_z, M_x, M_y, N_x, N_y terms) and is verified Hermitian; zero-strain limit reduces to H₀ = D S_z² + γB_z S_z.
- With N_x = N_y = 0 and B_z = 0, eigenvalues match Alam et al. Eq. (8) to machine precision.
- **Full 5-parameter case (N_x, N_y ≠ 0):**
  - Shift and splitting validated against hand-derived 2nd-order non-degenerate perturbation theory in `tests/test_perturbative_check.py` (0.0001% relative error vs. exact `eigh()` diagonalization; pass threshold tightened from 5% to 0.1%).
  - General `imbalance()` validated via ΔM = ±1 selection rules (zero 1st-order N-correction) and verified numerically to scale cleanly as O(N²) in `scripts/check_imbalance_N_scaling.py`.

### Sign convention (resolved, documented)
- Our S_y matrix uses the standard spin-1 convention, yielding a global sign flip on imbalance I relative to the paper's Eq. (11) (confirmed at 4 values of φ_MW), consistent with a y → −y / S_y phase convention difference.
- Does not affect eigenvalue-only quantities (D, M_z, splitting magnitude) — only the sign of φ_str.
- Checked against Position 25 ODMR data: with `B_lab` in the `(T, V)` plane (along `V`), the predicted lower-frequency dip is deeper than the upper-frequency dip, matching the measured spectrum directly.

### Validated (Cases 1–3 — single-NV strain sweeps)
- **Case 1 — volumetric strain (`scripts/case1_volumetric.py`):** Swept ε_xx = ε_yy = ε_zz = e_v over ±50×10⁻⁵ (B_z = 0). Linear shift (~±5.2 MHz, compressive strain → positive shift) matches hand-calculated slope (2·h41 + h43 = −10.54 GHz/strain); splitting and imbalance are identically zero. *(Minor cosmetic caveat: x-axis `ev` is per-component strain rather than trace 3·e_v = Tr(ε).)*
- **Case 2 — pure shear, ε_yz (`scripts/case2_shear_yz.py`):** Swept ε_yz over ±50×10⁻⁵. Splitting matches the N = 0 formula to 0.0000% relative deviation; tiny shift (max ~262 Hz) and imbalance (max ~1.8×10⁻⁵) residuals match expected 2nd-order N_y coupling.
- **Case 3 — in-plane shear, ε_xy (`scripts/case3_shear_xy.py`):** Fixed background ε_xx = −ε_yy = 10⁻⁴, swept ε_xy over ±50×10⁻⁵ (51-point odd grid so ε_xy = 0 is sampled exactly). Splitting follows 2·√(M_x² + M_y²) with exact 0.5700 MHz floor at ε_xy = 0; imbalance peaks at |I| ≈ 1.0 at ε_xy = 0 and drops to ≈ 0.2 at sweep edges.
- Direct rasterized Fig. 3 visual overlay explicitly deferred as redundant given analytical and numerical precision.

### Four-orientation ensemble model (`src/transformations.py`, `src/lineshape.py`, `tests/test_lineshape.py`)
- **Appendix B transformations (`src/transformations.py`):** Coded T⁽¹⁾–T⁽⁴⁾ along with `rotate_strain()` and `rotate_field()`; verified orthogonal and proper rotations (det = +1) in `tests/test_transformations.py`.
- **Lorentzian lineshapes (`src/lineshape.py`):** Implemented `single_nv_lineshape()` (Eq. 6) and `ensemble_lineshape()` (Eq. 7, with optional `D_val`) using per-orientation amplitudes (A₁…A₄) and a single outer baseline restoration F(ν) = 1 − Σ_n A_n [α_{n−} L(ν − ν_{n−}) + α_{n+} L(ν − ν_{n+})]. Validated in `tests/test_lineshape.py` (baseline = 1.0, zero-strain and volumetric-strain collapse, component summation, and single-NV asymmetry).

### Data fitting & 2.49 µm milestone (`src/fitting.py`, `tests/test_fitting.py`, `scripts/fit_milestone.py`)
- Loaded `zenodo-data-set/Preprocessed_data.txt` (47 positions, 2.84–2.9045 GHz) and implemented variable-projection nonlinear least squares (`fit_spectrum`) for `(ε_TT, ε_VV, ε_TV)` (and optional `D`).
- Verified `(T, L, V)` aligned with crystal `(X, Y, Z)` and `B_lab` in the `(T, V)` plane (`B_lab = (0, 0, 1)` along `V`); under this geometry `ε_TV → −ε_TV` swaps `(NV1, NV3) ↔ (NV2, NV4)` identically, explaining the paper's `±8.4×10⁻⁵` sign ambiguity.
- **Position 25 (2.49 µm) comparison:**
  - Paper reported: `ε_TT = −103.9×10⁻⁵`, `ε_VV = −27.1×10⁻⁵`, `ε_TV = ±8.4×10⁻⁵`, `D = 2.86670 GHz`
  - Fit (`D = 2.86670` fixed): `ε_TT = −103.29×10⁻⁵`, `ε_VV = −29.61×10⁻⁵`, `ε_TV = ±8.20×10⁻⁵`, `rmse = 0.00240`
  - Fit (`D` free): `ε_TT = −98.02×10⁻⁵`, `ε_VV = −17.34×10⁻⁵`, `ε_TV = ±9.51×10⁻⁵`, `D = 2.86760 GHz`, `rmse = 0.00213`

## Status checklist

### Ensemble validation & MW field orientation
- [x] Validate ensemble lineshapes (`tests/test_lineshape.py`) and resolve `B_MW` orientation against measured data (`B_lab` in `(T, V)` plane / along `V`).

### Data-fitting milestone
- [x] Load the Zenodo dataset (`zenodo-data-set/Preprocessed_data.txt`, DOI 10.5281/zenodo.10629041).
- [x] Write nonlinear least-squares fitting code (`src/fitting.py`, `tests/test_fitting.py`) and verify asymmetry against Position 25 data.
- [x] Reproduce the paper's reported fit at the 2.49 µm position: `ε_TT = −103.9×10⁻⁵`, `ε_VV = −27.1×10⁻⁵`, `ε_TV = ±8.4×10⁻⁵`, `D = 2.86670 GHz`.

### Optional cleanup & open questions
- [ ] Relabel Case 1's x-axis (`ev` vs. `3*ev` trace definition) if needed for direct slope comparisons.
- [ ] Resolve whether the "small splitting" claimed in the text of Fig. 5.1.0.1 (lab manual) is visible in the plotted dip or only present in unplotted raw data.