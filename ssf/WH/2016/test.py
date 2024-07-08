import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the sum of three Gaussians function
def triple_gaussian(x, a1, b1, c1, a2, b2, c2, a3, b3, c3):
    return (a1 * np.exp(-((x - b1) ** 2) / (2 * c1 ** 2)) +
            a2 * np.exp(-((x - b2) ** 2) / (2 * c2 ** 2)) +
            a3 * np.exp(-((x - b3) ** 2) / (2 * c3 ** 2)))

# List of ROOT files and corresponding mass points
files = [f'output_VHToAA2B2G_M{mass}_pythia8_ggh.root' for mass in range(20, 65, 5)]
mass_points = list(range(20, 65, 5))

# Plot setup
plt.figure(figsize=(12, 8))

# Process each file and mass point
for file, mass in zip(files, mass_points):
    # Open ROOT file and extract data and weights
    with uproot.open(file) as f:
        tree_name = f"ggh_{mass}_13TeV_AC_Tag0"  # Adjust with the format of your tree name
        tree = f[tree_name]
        x_data = tree["CMS_hgg_mass"].array(library="np")  # Adjust with your branch name for data
        weights = tree["weight"].array(library="np")  # Adjust with your branch name for weights

    # Generate weighted histogram from data
    hist, bin_edges = np.histogram(x_data, bins=550, range=(15, 70), weights=weights*36)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Fit the data
    p0 = [10, mass - 5, 2, 15, mass, 2, 7, mass + 5, 2]  # Initial parameter guesses
    popt, _ = curve_fit(triple_gaussian, bin_centers, hist, p0=p0, maxfev=10000)

    # Plot data and fit
    plt.scatter(bin_centers, hist, s=30, marker='o', label=f'Data Mass {mass}')
    plt.plot(bin_centers, triple_gaussian(bin_centers, *popt), label=f'Fit Mass {mass}', linewidth=2)

# Finalize plot
plt.ylim(0, 1)
plt.xlim(15,70)
plt.xlabel(r'M$_{\gamma \gamma}$ [GeV]')
plt.ylabel('Events/(0.1 GeV)')
plt.title('Simultaneous Fit')
plt.legend()
#plt.grid(True)

# Save the figure as a PNG file
plt.savefig('fit_plot.png', dpi=300)
plt.show()
