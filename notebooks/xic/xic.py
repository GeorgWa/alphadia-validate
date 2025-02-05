import numpy as np
from matplotlib import pyplot as plt

from .utils import normalize_profiles, median_axis
# intensity_slice with dimensions
# 0: mz
# 1: rt


def plot_median_xic(spectrum_slice):
    intensity_slice = spectrum_slice[0].sum(axis=1).sum(axis=1)
    normalized_intensity_slice = normalize_profiles(intensity_slice)
    median_profile = median_axis(normalized_intensity_slice, axis=0)
    # corr_list = correlation_coefficient(median_profile, intensity_slice)

    for i in range(normalized_intensity_slice.shape[0]):
        plt.plot(normalized_intensity_slice[i])

    plt.show()

    for fragment in range(normalized_intensity_slice.shape[0]):
        corr = np.corrcoef(normalized_intensity_slice[fragment], median_profile)[0, 1]
        print(corr)
        plt.plot(normalized_intensity_slice[fragment], color="gray", alpha=0.5)

    plt.plot(median_profile, color="red", linewidth=2)
    plt.show()
