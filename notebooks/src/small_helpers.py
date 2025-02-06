import numpy as np
import pandas as pd


def get_random_precursor_hash(
    precursor_df: pd.DataFrame, quantile: int, *, n_quantiles: int = 100
) -> int:
    """Get a random precursor hash, defined by quantile of q-value.
    quantile = 0, 1, 2, .. -> "pretty ones"
    quantile = .., 97, 98, 99 -> "pretty ones"
    """

    precursor_df_sorted = precursor_df.sort_values("qval")

    # Calculate the quantiles
    quantiles = np.array_split(precursor_df_sorted, n_quantiles)

    # Select a random precursor from the 3rd quantile
    random_precursor = quantiles[quantile].sample(n=1)
    print(
        "picked",
        random_precursor.index.values[0],
        random_precursor["sequence"].values[0],
        "qval=",
        random_precursor["qval"].values[0],
    )
    return random_precursor["mod_seq_charge_hash"].values[0]
