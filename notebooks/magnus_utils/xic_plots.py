from matplotlib import pyplot as plt


def plot_all_xics(dense, n_cols=4):
    """Plot the XICs of the observed precursors in separate panels."""
    default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color'] * 10

    xic_observed = dense[0].sum(axis=(1, 2))

    n_xics = xic_observed.shape[0]

    n_rows = (n_xics + 1) // n_cols + 1
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 10))
    axes = axes.flatten()

    # all in one
    ax = axes[0]
    for i in range(n_xics):
        ax.plot(xic_observed[i], color=default_colors[i])

    # each in a separate subplot
    for i in range(n_xics):
        ax = axes[i + 1]
        ax.plot(xic_observed[i], color=default_colors[i])