# scripts/fit_milestone.py
import matplotlib.pyplot as plt

from src.fitting import fit_spectrum, load_preprocessed_data

freqs, spectra = load_preprocessed_data()

# Preprocessed_data.txt uses 1-indexed headers (Position N = Table II position no. N-1)
y24 = spectra[:, 23]  # Position 24 -> position no. 23 (0.0 um, Fig. 4)
y25 = spectra[:, 24]  # Position 25 -> position no. 24 (0.83 um)
y27 = spectra[:, 26]  # Position 27 -> position no. 26 (2.49 um)

res24_free = fit_spectrum(freqs, y24, B_lab=(0.0, 0.0, 1.0), D_init=2.8670, fit_D=True, equal_A=True)
res24_fix = fit_spectrum(freqs, y24, B_lab=(0.0, 0.0, 1.0), D_init=2.86670, fit_D=False, equal_A=True)
res25_free = fit_spectrum(freqs, y25, B_lab=(0.0, 0.0, 1.0), D_init=2.8670, fit_D=True, equal_A=True)
res25_fix = fit_spectrum(freqs, y25, B_lab=(0.0, 0.0, 1.0), D_init=2.86670, fit_D=False, equal_A=True)
res27_free = fit_spectrum(freqs, y27, B_lab=(0.0, 0.0, 1.0), D_init=2.8670, fit_D=True, equal_A=True)

for label, r in [
    ("Pos 24 (0.0 um, D free)", res24_free),
    ("Pos 24 (0.0 um, D=2.86670 fixed)", res24_fix),
    ("Pos 25 (0.83 um, D free)", res25_free),
    ("Pos 25 (0.83 um, D=2.86670 fixed)", res25_fix),
    ("Pos 27 (2.49 um, D free)", res27_free),
]:
    print(
        f"{label:34s}: eps_TT={r.eps_TT*1e5:7.2f}e-5, "
        f"eps_VV={r.eps_VV*1e5:6.2f}e-5, eps_TV=+/-{abs(r.eps_TV)*1e5:5.2f}e-5, "
        f"D={r.D:.5f} GHz, rmse={r.rmse:.5f}"
    )

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(freqs, y24, "o", ms=3, alpha=0.6, label="Preprocessed data (Pos 24 / 0.0 um)")
ax.plot(freqs, res24_free.y_fit, "-", lw=2, color="green", label="Ensemble fit F(nu)")
colors = ["orange", "red", "black", "blue"]
for i, (comp, col) in enumerate(zip(res24_free.components, colors), start=1):
    ax.plot(freqs, comp, "--", lw=1, color=col, alpha=0.8, label=f"NV{i}")
ax.set_xlabel("MW frequency (GHz)")
ax.set_ylabel("Normalized PL")
ax.set_title("Ensemble Fit at Position 24 (0.0 um, Fig. 4)")
ax.legend(ncol=2, fontsize=8)
plt.tight_layout()
plt.savefig("milestone_fit_pos24.png", dpi=150)
