import altair as alt
import pandas as pd


def _get_mirrorplot_df(
    mz_library,
    intensity_library,
    spectrum_slice,
):
    df_lib = pd.DataFrame(
        {
            "m/z": mz_library,
            "Intensity": -intensity_library / intensity_library.max(),
            "Scan": -1,
        }
    )
    dense = spectrum_slice[0]
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
    return chart


def plot_mirror_2(spectrum_slice, mz_library, intensity_library, width=800, height=300):
    df_ = _get_mirrorplot_df(mz_library, intensity_library, spectrum_slice)
    return _plot_mirror(df_).properties(width=width, height=height)
