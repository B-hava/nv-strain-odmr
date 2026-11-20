## Directory Structure 
(venv) bubba@bubba-ThinkPad-T490:~/Desktop/nv-strain-odmr$ tree -L 2
.
├── case1_volumetric.png
├── case2_shear_yz.png
├── NOTES.md
├── requirements.txt
├── scripts
│   ├── case1_volumetric.png
│   ├── case1_volumetric.py
│   ├── case2_shear_yz.py
│   ├── debug_imbalance_sign.py
│   ├── __init__.py
│   └── __pycache__
├── src
│   ├── constants.py
│   ├── hamiltonian.py
│   ├── __init__.py
│   ├── observables.py
│   └── __pycache__
├── tests
│   ├── __init__.py
│   ├── __pycache__
│   ├── test_hamiltonian.py
│   └── test_observables.py
└── venv
    ├── bin
    ├── include
    ├── lib
    ├── lib64 -> lib
    ├── pyvenv.cfg
    └── share


## Phase 0 validation — 8/11/26
- Hamiltonian matches hand-derived form (D, γ, M_z, M_x, M_y, N_x, N_y terms), verified Hermitian.
- Zero-strain limit reduces exactly to H0 = D Sz² + γBz Sz.
- N_x=N_y=0, B_z=0 case matches Alam et al. Eq. (8) closed-form eigenvalues to machine precision.
- Test parameter values (Mx=0.3 etc.) are arbitrary, not physical — chosen only to exercise the algebra.
- Not yet validated: full 5-parameter case against anything (no closed form exists for N_x,N_y≠0 to check against).


## Sign convention note — imbalance / phi_str
- Our S_x, S_y spin-1 matrices use the standard convention (S_y has the
  i/-i off-diagonal structure), but the paper never gives explicit
  matrices, only the operator algebra.
- Our computed imbalance I matches the paper's Eq.(11) functional form
  (|cos(2*phi_MW + phi_str)|) exactly, but with the opposite overall
  sign, confirmed by direct comparison at 4 values of phi_MW
  (0, pi/4, pi/2, 3pi/4) — a clean global negation, not a bug.
- Likely cause: relative phase/sign convention for S_y (equivalent to
  a y -> -y redefinition) differing between our matrices and whatever
  the paper implicitly used. This does not affect any eigenvalue-only
  result (D, M_z, splitting magnitude Delta) — only sign of phi_str.
- Resolution: adopted our own internally consistent convention;
  documented rather than "fixed," since there's no way to verify
  which convention the paper intended without their source code.
- Implication for later: if/when fitting real ODMR imbalance data
  (Phase 1), be aware the sign of predicted imbalance vs. measured
  asymmetry may need this same flip applied — check against a known
  measured sign if possible rather than assuming.

  ## Case 1 validation (volumetric strain) — 18/11/26
- Swept exx=eyy=ezz=ev over ±50e-5, Bz=0, via scripts/case1_volumetric.py.
- Shift: linear in strain, ~±5.2 MHz over the range. Sign confirmed against
  Alam et al.: compressive (negative) strain → positive shift (higher
  energy), matching their Fig. 3(I-a)/(I-b) description.
- Magnitude check: dMz/dev = 2*h41 + h43 = -10.54 GHz/strain (hand calc from
  H_IJ table) → predicts ~5.27 MHz at ev=5e-4, matches plot.
- Splitting and imbalance: exactly zero across the full range, as expected
  (no shear content in volumetric strain → no M_x, M_y → no coupling
  between |±1>).
- Caveat (not yet fixed, cosmetic only): script's x-axis variable `ev` is
  the individual diagonal strain component (exx=eyy=ezz=ev), so the actual
  trace being swept is 3*ev, not ev. Paper's ε_v = Tr(ε) directly. Doesn't
  affect sign/shape/magnitude checks above, but would need correcting
  (relabel axis as 3*ev, or redefine ev=trace/3) before overlaying directly
  on Fig. 3(I-b) for a slope comparison.

## Case 2 validation (pure shear, ε_yz only) — 18/11/26
- Swept eyz over ±50e-5, all other strain components zero, via
  scripts/case2_shear_yz.py.
- Note: unlike Case 1, this case has Nx=0 but Ny = 0.5*h26*eyz ≠ 0 — the
  paper's closed-form Eqs. (8)/(11)/(13) assume Nx=Ny=0 and treat N-terms
  as a dropped second-order correction. Our hamiltonian() includes them
  exactly, so small deviations from the paper's approximate formulas are
  expected, not bugs.
- Shift: max |shift| = 0.00026 MHz (262 Hz), not exactly zero. Attributed
  to the second-order N_y-mediated coupling of |0> into the |±1> manifold,
  which the paper's Eq. (8) neglects. ~1e-7 relative to D — consistent with
  a genuine second-order perturbative effect, not an error. 
- Splitting: 0.0000% relative deviation from the N=0 approximate formula
  (δν ≈ |h16 * eyz|). Splitting is essentially unaffected by the N-terms
  at this order — makes sense, as N perturbs E_gs/imbalance well before it
  perturbs the E_+ - E_- splitting itself.
- Imbalance: max |I| = 1.8e-5, not exactly zero as the paper's N=0 symmetry
  argument predicts for φ_MW=π/2. Same root cause as the shift residual —
  a real, tiny symmetry-breaking from N_y that the paper's simplified
  treatment doesn't capture. Negligible for any practical comparison to
  data, but worth remembering if N-term residuals ever compound
  noticeably across the four-orientation ensemble sum later.
- Open: plot (case2_shear_yz.png) not yet visually cross-checked against
  Fig. 3(II-a/c) for linearity shape — numbers alone were conclusive enough
  to proceed, visual check deferred.

  ## Case 3 validation (shear strain in xy-plane, fixed background) — 19/11/26
- Fixed εxx=-εyy=1e-4 (background), swept εxy over ±50e-5, via
  scripts/case3_shear_xy.py.
- Shift: max |shift| = 0.00023 MHz, same order as Case 2's N-term residual
  — consistent, not a new effect.
- Splitting floor at εxy=0: predicted 0.5700 MHz (from Mx background alone,
  N=0 approx), measured 0.5730 MHz — within 0.5%, gap attributable to the
  Nx contribution from the same fixed background (same mechanism as Case
  2's deviation, scaled up by the larger background here). Splitting curve
  shape (rounded-V / hyperbola-like, not linear) matches expected
  δν ≈ 2*sqrt(Mx² + My²) with Mx pinned and My sweeping.
- Imbalance: NOT negligible here, unlike Cases 1-2. Peaks at |I|≈1.0 at
  εxy=0, falls to ≈0.2 at sweep edges. Confirmed consistent with the
  Eq.(11)-style functional form (I ∝ cos(2*phi_MW + phi_str)) validated
  earlier: at εxy=0, My=0 so phi_str is pinned at 0/pi, giving |I|≈1 at
  phi_MW=pi/2; as εxy grows, My becomes comparable to Mx, phi_str rotates,
  cosine drops. Not a bug — first case where imbalance is physically
  significant rather than near-zero.
- Important flag: this is the FIRST case where the documented sign-
  convention ambiguity (global negation of phi_str vs. the paper, see
  "Sign convention note" above) could actually matter for the *predicted
  sign* of a physically significant imbalance, not just its magnitude.
  Cases 1-2 couldn't expose this since their imbalance was ~0 regardless
  of sign convention. Relevant for the ensemble model and any later fit
  to real imbalance data — check sign against a known reference then.

## Decision: skipping visual Fig. 3 overlay for now — 19/11/26
- Considered pulling the actual Fig. 3 panels from the paper PDF to
  visually/numerically overlay against Cases 1-3's plots.
- Judged not worth it at this stage: we've already validated against the
  closed-form equations to machine precision (Case 1) and cross-checked
  magnitudes/signs by hand for all three cases. Reading values off a
  rasterized PDF figure would have *lower* precision than what we already
  have, so the main value would just be catching a gross error (wrong
  axis/sign/shape) — already ruled out another way.
- Would become worth it if we specifically needed to check whether our
  small N-term residuals (shift ~260 Hz, imbalance ~1.8e-5 in Case 2) are
  visible in the paper's own plotted curves, i.e. whether their figure was
  generated by exact diagonalization or by the simplified N=0 formulas.
  Deferred as a stretch goal, not a blocker.
- Proceeding instead to the two genuinely unvalidated/unimplemented items:
  test_perturbative_check.py and the four-orientation ensemble model.

  ## Independent perturbative cross-check (N_x, N_y ≠ 0) — 19/11/26
- Hand-derived leading-order (2nd-order) non-degenerate perturbation
  theory formulas for shift/splitting with Nx,Ny≠0, treating N-coupling
  as a perturbation on the exactly-solvable N=0 problem (Eq. 8):
      shift     ≈ Mz + (3/2)|N_perp|^2 / D
      splitting ≈ 2|M_perp| + |N_perp|^2 * cos(phi_M + 2*phi_N) / D
  where M_perp = Mx + iMy, N_perp = Nx + iNy, phi_M/phi_N their phases.
  Derivation independent of src/hamiltonian.py's own algebra — see
  tests/test_perturbative_check.py docstring for the formula, full
  derivation in chat history.
- This formula retroactively explains two earlier empirical results
  without having been fit to them:
  - Case 2 (pure eyz shear) showed exactly 0.0000% splitting deviation
    from the N=0 approx. The formula predicts this exactly: My,Ny are
    both purely imaginary there (phi_M=phi_N=pi/2), so
    cos(phi_M+2*phi_N)=cos(3pi/2)=0 — the cross term vanishes identically,
    not by coincidence.
  - Case 3 (shear in xy, fixed background) showed a ~0.5% splitting
    deviation at exy=0. There Mx,Nx are purely real (My=Ny=0 at that
    point), so phi_M=phi_N=0, cos(0)=1 — maximal, consistent with the
    nonzero deviation observed.
- Test (tests/test_perturbative_check.py::test_perturbative_matches_exact_small_N)
  compares this formula against exact eigh() diagonalization at generic,
  non-axis-aligned Mz,Mx,My,Nx,Ny values (deliberately chosen so the
  phase-dependent cross term is actually exercised).
- Result: rel_err = 0.0001% for both shift and splitting — far better
  than the 5% provisional tolerance. Consistent with 2nd-order PT error
  scaling as (N/D)^2 ~ 1e-6 given N~0.002-0.003 GHz vs D~2.87 GHz.
- Tightened the test's pass threshold from 5% -> 0.1% (1000x margin
  above observed error) now that the real accuracy is known — safe to
  tighten since the test uses fixed deterministic parameters, no
  flakiness risk.
- Status: this closes the "test_perturbative_check.py ... explicitly
  not implemented" item. The general imbalance() formula for Nx,Ny≠0
  still has no independent check (separate, still-open item).

  ## Correction to Case 3 notes (19/11/26) — supersedes the "splitting floor" entry above
- The ~0.5% gap between the predicted N=0 floor (0.5700 MHz) and the
  measured value (0.5730 MHz) was INCORRECTLY attributed to the Nx
  background contribution in the original Case 3 entry above. That
  attribution was wrong.
- Real cause: exy_range = linspace(-50e-5, 50e-5, 50) has 50 (even)
  points, so no sample lands exactly at εxy=0 — the nearest grid point
  is at εxy ≈ +1.02e-5, not 0. At that point My ≈ 0.0291 MHz is already
  nonzero, and 2*sqrt(Mx²+My²) ≈ 0.573 MHz — which is the entire "gap,"
  with no N-physics involved.
- Confirmed directly: evaluating exactly at εxy=0 gives splitting =
  0.5700 MHz (matches the N=0 floor prediction exactly) and shift =
  0.000009 MHz, which matches the hand-derived perturbative formula's
  prediction (1.5*Nx²/D ≈ 0.0000088 MHz using h25=-2.6 GHz) almost
  exactly. True N-term effect here is ~Hz-scale, negligible, as expected
  from N/D ~ 1e-4-1e-5 in this regime.
- Fix applied: exy_range changed to linspace(-50e-5, 50e-5, 51) (odd
  point count) so the sweep grid includes an exact εxy=0 sample.
- Lesson: don't attribute a numerical discrepancy to new physics before
  checking whether it's a grid/sampling artifact — should have verified
  order of magnitude against the perturbative formula before writing the
  original explanation into these notes.


  ## General imbalance() validation for Nx,Ny≠0 — 19/11/26
- Derived (by hand) that the first-order eigenvector correction to the
  ground state from N-mixing produces NO first-order correction to the
  imbalance matrix elements <gs|S.B|a>, <gs|S.B|b>. Structural reason:
  Sx,Sy only connect |0> <-> |+1>/|-1> (ΔM=±1 selection rule built into
  spin-1 operators), so every first-order correction term requires
  <a|S.B|a> or <b|S.B|a>, both of which vanish identically for all
  phi_MW, since |a>,|b> have zero |0>-component and Sx,Sy map them only
  into the |0> direction.
- This predicts the leading N-correction to imbalance is second order
  (~N²), not first order.
- Naive estimate of the coefficient (treating it as purely (N/D)²)
  under-predicted Case 2's observed deviation (1.8e-5) by ~300x,
  suggesting the full second-order coefficient has additional structure
  (likely extra 1/|M_perp| factors near the Mx,My-split a/b states) —
  not derived in full; judged not worth the additional algebra given the
  effect is already shown negligible by direct exact-diagonalization
  measurement.
- Instead validated the claim that actually matters — that the
  correction is higher-order, not first-order — via a direct numerical
  scaling test (scripts/check_imbalance_N_scaling.py): fixed generic
  Mx,My background, varied N by halving repeatedly, measured
  |I(N)-I(0)|.
- Result: ratio of successive deviations = 4.000 at every halving,
  exactly matching quadratic scaling, confirming the structural
  derivation and ruling out a first-order leak.
- Status: closes the "general imbalance() formula ... has no independent
  check" item. Confirmed via scaling argument + structural derivation,
  not a closed-form second-order formula (judged not worth deriving
  given diminishing returns — exact magnitudes already in hand from
  Cases 2-3's direct numerics).

  ## Everything between here and the fitting milestone, in dependency order:

1. Ensemble (four-orientation) model — none of it exists yet:

Code the four Appendix B transformation matrices T⁽¹⁾, T⁽²⁾, T⁽³⁾, T⁽⁴⁾ (rotate crystal-frame strain/MW field into each NV symmetry axis)
Implement the per-orientation rotation: ε⁽ⁿ⁾ = T⁽ⁿ⁾† ε^(cryst) T⁽ⁿ⁾ (and the equivalent for the MW field direction)
Implement the single-NV Lorentzian lineshape f(ν) (paper’s Eq. 6) — you have eigenvalues and transition rates/imbalance from the single-NV model, but nothing that turns those into an actual spectral lineshape (linewidth, amplitude, position) yet
Implement the ensemble sum F(ν) = Σᵢ fᵢ(ν) (Eq. 7) — summing the four orientations’ lineshapes into one simulated spectrum

2. Validation of the ensemble model — not listed in your original notes explicitly, but implied by the pattern we’ve followed so far: before trusting the four-orientation sum enough to fit real data with it, we’d want to sanity-check it against something — at minimum confirm that with all four orientations given identical strain (should reduce to a known single-NV-like limit), or against any ensemble-level closed-form/figure the paper gives, the same way we validated the single-NV model against Eq. 8 and Fig. 3. Worth deciding how rigorously to do this before moving on, rather than shipping straight to fitting with an unchecked model.

3. Data fitting infrastructure:

Fetch/load the actual Zenodo dataset (DOI 10.5281/zenodo.10629041) — not yet touched
Write the nonlinear least-squares fitting code (likely scipy.optimize.curve_fit or least_squares, fitting F(ν) against the real spectrum, with εxx,εyy,εzz,εxy,εxz,εyz,D — or however they parameterize it — as free parameters)

4. The actual milestone: run that fit at the 2.49 µm position and check whether it reproduces the paper’s reported values (ε_TT=-103.9×10⁻⁵, ε_VV=-27.1×10⁻⁵, ε_TV=±8.4×10⁻⁵, D=2.86670 GHz).

## Appendix B transformation matrices — coded & validated — 19/11/26
- Transcribed T^(1)-T^(4) from the paper's Appendix B into
  src/transformations.py.
- Caught and corrected a transcription slip during transcription: my
  first reading of T^(2) row 1, middle entry, was -1/√2, which made row
  1 exactly -1x row 0 (linearly dependent — impossible for a rotation
  matrix). Cross-checked two ways: (a) rederiving T^(2) via the stated
  rule T^(n)=R^(n)T^(1) (π-rotation about the lab X-axis, flipping rows
  1,2 of T^(1)), and (b) same rule applied to T^(3)/T^(4) (π about Y/Z)
  reproduced the image exactly for those two. Both pointed to +1/√2.
  User confirmed against the actual PDF: entry is indeed +1/√2.
- tests/test_transformations.py checks all four are orthogonal
  (T^T T = I) and proper rotations (det = +1), not reflections. Both
  pass for all four matrices.
- Also added rotate_strain(eps_cryst, T) = T.T @ eps_cryst @ T and
  rotate_field(B_cryst, T) = T.T @ B_cryst, implementing the Appendix B
  transformation rules directly.

  ## Eq. 6/7 baseline resolution — 19/11/26
- Main text's f(nu)=1-A[...] (Eq.6) summed directly as F(nu)=sum f_i(nu)
  (Eq.7) would give a baseline of 4 far from resonance, contradicting
  Fig.4's ~1.0 baseline. Confirmed with user (checked SI): this is a
  notational slip in the main text. Actual fitting form shifts baseline
  to 0 before fitting, sums the four dip terms, then restores baseline
  to 1 ONCE outside the sum:
      F(nu) = 1 - sum_n A_n[alpha_n- L(nu-nu_n-) + alpha_n+ L(nu-nu_n+)]
  Also: A is per-orientation (A_1..A_4, Table II), not one shared
  amplitude across all four NVs.
- Implemented in src/lineshape.py: single_nv_lineshape() (individual
  orientation's own "1-A[...]" curve, matches Fig.4's colored lines) and
  ensemble_lineshape() (correct baseline-once-restored total, matches
  Fig.4's green curve).
- Per-orientation A_n being a free/fit parameter means we don't need the
  true physical |B_perp| magnitude per orientation — only B_MW's
  direction matters, to set how each orientation's amplitude splits
  between its +/- lines (the alpha_-,alpha_+ ratio). Open question
  (MW field direction: TV-plane vs along L) still unresolved; next step
  is forward-simulating with Fig.4's known strain values under a few
  candidate directions and comparing against the published curves.