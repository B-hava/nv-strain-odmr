# scripts/fit_milestone.py
import matplotlib.pyplot as plt

from src.fitting import fit_spectrum, load_preprocessed_data

freqs, spectra = load_preprocessed_data()

# Position 25 (0-indexed 24) corresponds to the 2.49 um center position
y25 = spectra[:, 24]

res_fixed = fit_spectrum(
    freqs, y25, B_lab=(0.0, 0.0, 1.0), D_init=2.86670, fit_D=False, equal_A=True
)
res_free = fit_spectrum(
    freqs, y25, B_lab=(0.0, 0.0, 1.0), D_init=2.86670, fit_D=True, equal_A=True
)

print("Paper reported (2.49 um): eps_TT=-103.9e-5, eps_VV=-27.1e-5, eps_TV=+/-8.4e-5, D=2.86670 GHz")
print(
    f"Fit (D=2.86670 fixed):    eps_TT={res_fixed.eps_TT*1e5:.2f}e-5, "
    f"eps_VV={res_fixed.eps_VV*1e5:.2f}e-5, eps_TV=+/-{abs(res_fixed.eps_TV)*1e5:.2f}e-5, "
    f"rmse={res_fixed.rmse:.5f}"
)
print(
    f"Fit (D free):             eps_TT={res_free.eps_TT*1e5:.2f}e-5, "
    f"eps_VV={res_free.eps_VV*1e5:.2f}e-5, eps_TV=+/-{abs(res_free.eps_TV)*1e5:.2f}e-5, "
    f"D={res_free.D:.5f} GHz, rmse={res_free.rmse:.5f}"
)

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(freqs, y25, "o", ms=3, alpha=0.6, label="Preprocessed data (Pos 25)")
ax.plot(freqs, res_fixed.y_fit, "-", lw=2, color="green", label="Ensemble fit F(nu)")
colors = ["orange", "red", "black", "blue"]
for i, (comp, col) in enumerate(zip(res_fixed.components, colors), start=1):
    ax.plot(freqs, comp, "--", lw=1, color=col, alpha=0.8, label=f"NV{i}")
ax.set_xlabel("MW frequency (GHz)")
ax.set_ylabel("Normalized PL")
ax.set_title("Milestone Fit at Position 25 (2.49 um)")
ax.legend(ncol=2, fontsize=8)
plt.tight_layout()
plt.savefig("milestone_fit_pos25.png", dpi=150)
