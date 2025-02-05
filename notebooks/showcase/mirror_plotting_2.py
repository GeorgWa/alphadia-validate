import altair as alt
import numpy as np
import pandas as pd


def get_library_entry_by_hash(speclib, _hash, min_intensity=0.01):
    speclib_entry = speclib.precursor_df[
        speclib.precursor_df["mod_seq_charge_hash"] == _hash
    ].iloc[0]

    fragment_mz = (
        speclib.fragment_mz_df.iloc[
            speclib_entry.frag_start_idx : speclib_entry.frag_stop_idx
        ]
        .to_numpy()
        .flatten()
    )
    fragment_intensity = (
        speclib.fragment_intensity_df.iloc[
            speclib_entry.frag_start_idx : speclib_entry.frag_stop_idx
        ]
        .to_numpy()
        .flatten()
    )
    fragment_mask = fragment_intensity > min_intensity

    fragment_mz = fragment_mz[fragment_mask]
    fragment_intensity = fragment_intensity[fragment_mask]

    # sort both by mz
    fragment_order = np.argsort(fragment_mz)
    fragment_mz = fragment_mz[fragment_order]
    fragment_intensity = fragment_intensity[fragment_order]

    return speclib_entry, fragment_mz, fragment_intensity


def get_observed_intensity_dense(
    dia_data, mz_library, speclib_entry, precursor_entry, tolerance=30
):
    jit_data = dia_data.jitclass()
    precursor_query = np.array(
        [[speclib_entry.precursor_mz, speclib_entry.precursor_mz]], dtype=np.float32
    )
    scan_limits = np.array(
        [[precursor_entry.scan_start, precursor_entry.scan_stop, 1]], dtype=np.int64
    )
    frame_limits = np.array(
        [[precursor_entry.frame_start, precursor_entry.frame_stop, 1]], dtype=np.int64
    )

    dense, _ = jit_data.get_dense(
        frame_limits,
        scan_limits,
        mz_library,
        tolerance,
        precursor_query,
    )
    return dense[0]


def get_mirrorplot_df(
    mz_library, intensity_library, dia_data, speclib_entry, precursor_df, _hash
):
    df_lib = pd.DataFrame(
        {
            "m/z": mz_library,
            "Intensity": -intensity_library / intensity_library.max(),
            "Scan": -1,
        }
    )

    precursor_entry = precursor_df[precursor_df["mod_seq_charge_hash"] == _hash].iloc[0]

    dense = get_observed_intensity_dense(
        dia_data, mz_library, speclib_entry, precursor_entry
    )
    intensity_by_scan = dense.sum(axis=(1, 2))
    max_intensity = dense.sum(axis=(1, 2, 3)).max()

    dfs_obs = []
    for i in range(dense.shape[3]):
        intensities = intensity_by_scan[:, i]
        df_ = pd.DataFrame(
            {
                "m/z": mz_library,
                "Intensity": intensities / max_intensity,
                "Scan": int(i),
            }
        )
        dfs_obs.append(df_)

    df = pd.concat([df_lib, *dfs_obs])
    return df


def _plot_mirror(df_mirrorplot):
    chart = (
        alt.Chart(df_mirrorplot)
        .mark_bar(size=2)
        .encode(
            x=alt.X("m/z:Q", bin=False),
            y=alt.Y("Intensity:Q", stack=True),
            color=alt.Color("Scan:N"),
        )
    )
    chart
    return chart


def plot_mirror(spectral_library, precursor_df, dia_data, _hash):
    speclib_entry, mz_library, intensity_library = get_library_entry_by_hash(
        spectral_library, _hash
    )

    df_ = get_mirrorplot_df(
        mz_library, intensity_library, dia_data, speclib_entry, precursor_df, _hash
    )
    return _plot_mirror(df_).properties(width=700, height=300)
