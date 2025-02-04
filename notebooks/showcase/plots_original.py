from matplotlib import pyplot as plt


def plot_simple_mirror(dense, intensity_library, mz_library):
    """First we will select the intensity dimension and sum over all other dimensions but the fragment mz dimension."""
    intensity_observed = dense[0].sum(axis=(1, 2, 3))
    intensity_observed_normalized = intensity_observed / intensity_observed.max()
    intensity_library_normalized = intensity_library / intensity_library.max()

    plt.stem(mz_library, intensity_observed_normalized)
    plt.stem(mz_library, -intensity_library_normalized)
    plt.show()


def plot_simple_xic(dense):
    """
    Finally, we will visualize the Precusor ion chromatogram.
    We will again select the intensity dimension and sum over ion mobility and observation but leave the retention time dimension.

    :param dense:
    :return:
    """
    xic_observed = dense[0].sum(axis=(1, 2))
    for i in range(xic_observed.shape[0]):
        plt.plot(xic_observed[i])
    plt.show()
