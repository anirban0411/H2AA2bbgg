import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

# Define the sum of three Gaussians function
def triple_gaussian(x, a1, b1, c1, a2, b2, c2, a3, b3, c3):
    return (a1 * np.exp(-((x - b1) ** 2) / (2 * c1 ** 2)) +
            a2 * np.exp(-((x - b2) ** 2) / (2 * c2 ** 2)) +
            a3 * np.exp(-((x - b3) ** 2) / (2 * c3 ** 2)))

# List of ROOT files and corresponding mass points
original_mass_points = list(range(20, 65, 5))
files = [f'output_WHToAA2B2G_M{mass}_pythia8.root' for mass in original_mass_points]

# Store fitted parameters for interpolation
params_list = []

# Plot setup
plt.figure(figsize=(12, 8))

# Process each file and mass point
for file, mass in zip(files, original_mass_points):
    with uproot.open(file) as f:
        tree_name = f"ggh_{mass}_13TeV_H2AA2bbgg_Tag0"
        tree = f[tree_name]
        x_data = tree["CMS_hgg_mass"].array(library="np")
        weights = tree["weight"].array(library="np")

    hist, bin_edges = np.histogram(x_data, bins=550, range=(15, 70), weights=weights)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    p0 = [10, mass - 5, 2, 15, mass, 2, 7, mass + 5, 2]
    popt, _ = curve_fit(triple_gaussian, bin_centers, hist, p0=p0, maxfev=10000)
    params_list.append(popt)

    # Plot original mass points data and fit
    plt.scatter(bin_centers, hist, s=30, marker='o', label=f'Data Mass {mass}')
    plt.plot(bin_centers, triple_gaussian(bin_centers, *popt), label=f'Fit Mass {mass}', linewidth=2, linestyle='-')

# Interpolation for intermediate mass points
interp_mass_points = [38]
params_array = np.array(params_list)
interpolated_params = interp1d(original_mass_points, params_array, axis=0, kind='linear')

# Plot interpolated fits
for mass in interp_mass_points:
    popt = interpolated_params(mass)
    x_smooth = np.linspace(15, 70, 550)
    plt.plot(x_smooth, triple_gaussian(x_smooth, *popt), label=f'Interpolated Mass {mass}', linewidth=2, linestyle='--')

# Finalize plot
plt.ylim(0, 0.02)
plt.xlim(15, 70)
plt.xlabel(r'M$_{\gamma \gamma}$ [GeV]')
plt.ylabel('Events/(0.1 GeV)')
plt.title('Simultaneous Fit')
plt.legend()
plt.grid(True)

# Save the figure as a PNG file
plt.savefig('fit_interpolated_plot.png', dpi=300)
plt.show()
