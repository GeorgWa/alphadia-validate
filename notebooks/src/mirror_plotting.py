import altair
import numpy as np
import pandas as pd
import seaborn as sns

from .xic_utils import normalize_profiles, median_axis, correlation_coefficient


fragment_colos = {
    "a": "#388E3C",
    "b": "#1976D2",
    "c": "#00796B",
    "x": "#7B1FA2",
    "y": "#D32F2F",
    "z": "#F57C00",
    "unknown": "#212121",
    None: "#808080",
}


def _get_colors(n_colors, palette_name=None):
    if palette_name is not None:
        pal = sns.color_palette(palette_name, n_colors)
        colors = list(pal.as_hex())
    else:
        pal = sns.color_palette("Spectral", n_colors + 3)  # .as_hex
        colors = list(pal.as_hex())
        remove_middle_elements(colors)
    return colors


def remove_middle_elements(lst):
    middle_index = len(lst) // 2
    lst.pop(middle_index)
    lst.pop(middle_index - 1)
    lst.pop(middle_index - 2)
    # lst.pop(0) ## delete first 3 elements
    # lst.pop(0)
    # lst.pop(0)
    return lst


def split_dense_byRT(dense, mz_library, label_library):
    intensity_observed_full = dense[0].sum(
        axis=(1, 2)
    )  ## do not sum over rt have shape (int, rt)

    intensity_flat = intensity_observed_full.ravel()
    mz_flat = np.tile(mz_library, dense[0].shape[3])
    # intensity_flat.shape, mz_flat.shape

    # Recover the original (mz, rt) index from the flattened array
    # Flatten mz library
    reconstructed_indices = np.unravel_index(
        np.arange(intensity_flat.size), intensity_observed_full.shape
    )
    mz_flat = mz_library[reconstructed_indices[0]]
    label_flat = label_library[reconstructed_indices[0]]

    indices = np.indices(intensity_observed_full.shape)
    rt_indices = indices[1].ravel()  # Retention time indices
    # cylce_indices = ...

    df_obs = pd.DataFrame(
        {
            "mz": mz_flat,
            "frag_label": label_flat,
            "intensity": intensity_flat,
            "rt_idx": rt_indices,
            # 'cycle_idx': cycle_indices
        }
    )
    df_obs["rt_cat"] = pd.Categorical(df_obs.rt_idx)
    max_int = df_obs.groupby(["mz"])["intensity"].sum().reset_index()["intensity"].max()
    df_obs["norm_intensity"] = df_obs["intensity"] / max_int

    df_obs["frag_label"] = df_obs.apply(
        lambda x: "" if x.intensity <= 0 else x.frag_label, axis=1
    )

    return df_obs


def convert_library_to_df(mz_library, intensity_library):
    df_library = pd.DataFrame({"mz": mz_library, "intensity": intensity_library})
    df_library["norm_intensity"] = df_library["intensity"] / max(
        df_library["intensity"]
    )
    return df_library


def interpolate_colors_red_to_green(num_colors):
    colors = []
    for i in range(num_colors):
        red = int(255 * (i / (num_colors - 1)) ** 0.7)  # Interpolating red channel
        green = int(
            255 * (1 - (i / (num_colors - 1)) ** 0.3)
        )  # Interpolating green channel
        colors.append(f"#{red:02X}{green:02X}00")  # Keep blue at 0
    return colors[::-1]


def interpolate_colors_orange_yellow_green(num_colors):
    colors = []
    midpoint = num_colors // 2  # Halfway point (orange → yellow, then yellow → green)

    for i in range(num_colors):
        if i < midpoint:
            # Interpolate from Orange (#FFA500) to Yellow (#FFFF00)
            red = 255
            green = int(
                165 + (90 * (i / (midpoint - 1)) ** 0.5)
            )  # From 165 (orange) to 255 (yellow)
        else:
            # Interpolate from Yellow (#FFFF00) to Green (#00FF00)
            red = int(
                255 * (1 - ((i - midpoint) / (num_colors - midpoint - 1)) ** 0.5)
            )  # Red decreases
            green = 255  # Green stays at 255

        colors.append(f"#{red:02X}{green:02X}00")  # Keep blue at 0

    return colors


# Function to interpolate colors from red to yellow (existing function)
def interpolate_colors_red_to_yellow(num_colors):
    colors = []
    for i in range(num_colors):
        green = int(
            255 * (i / (num_colors)) ** 0.7
        )  # Linear interpolation for the green channel
        colors.append(
            f"#FF{green:02X}00"
        )  # Format as HEX, keeping red at 255 and blue at 0
    return colors


def interpolate_colors_green_to_yellow(num_colors):
    colors = []
    for i in range(num_colors):
        red = int(
            255 * (1 - (i / (num_colors - 1)) ** 0.7)
        )  # Interpolating red channel
        colors.append(f"#{red:02X}FF00")  # Keep green at 255, blue at 0
    return colors[::-1]


def interpolate_complex_colormap(num_colors):
    original_num_colors = num_colors
    num_colors = max(12, num_colors)
    colors = []
    segments = [
        ("#FFFF00", "#CCCC00"),  # Yellow → Dark Yellow
        ("#CCCC00", "#FF8C00"),  # Dark Yellow → Orange
        ("#FF8C00", "#00FF00"),  # Orange → Green (Middle)
        ("#00FF00", "#8B4513"),  # Green → Brown
        ("#8B4513", "#FFC0CB"),  # Brown → Pink
        ("#FFC0CB", "#ff0000"),  # Pink → Red
    ]

    num_segments = len(segments)
    colors_per_segment = num_colors // num_segments

    def interpolate_color(color1, color2, t):
        """Linearly interpolate between two hex colors"""
        c1 = [int(color1[i : i + 2], 16) for i in (1, 3, 5)]
        c2 = [int(color2[i : i + 2], 16) for i in (1, 3, 5)]
        return f"#{int(c1[0] + (c2[0] - c1[0]) * t):02X}{int(c1[1] + (c2[1] - c1[1]) * t):02X}{int(c1[2] + (c2[2] - c1[2]) * t):02X}"

    for start_color, end_color in segments:
        for i in range(colors_per_segment):
            t = i / (colors_per_segment - 1 if colors_per_segment > 1 else 1)
            colors.append(interpolate_color(start_color, end_color, t))

    # print(num_colors)
    # print(len(colors), colors)
    colors = list(set(colors))[:original_num_colors]
    return colors


def get_palette_for_RT(num_rt_scans):
    if num_rt_scans == 1:
        return ["#FF0000"]
    # left_cols = interpolate_colors_red_to_green(math.ceil(num_rt_scans/2.))[:-1]
    # right_cols = interpolate_colors_orange_yellow_green(math.floor(num_rt_scans/2.)+1)[::-1]
    # return (left_cols + right_cols)[::-1]
    return interpolate_complex_colormap(num_rt_scans)


def format_prec_entry(prec):
    if type(prec.mods) == "str":
        mod_label_mapping = {"Carbamidomethyl": "", "Oxidation": "(Ox)"}
        prec = precursor_df.iloc[9]
        seq = list(prec.sequence)
        mods = [m.split("@")[0] for m in prec.mods.split(";")]
        sites = np.array(prec.mod_sites.split(";"), dtype=int) - 1
        for i, s in enumerate(sites):
            mod = (
                mods[i]
                if mods[i] not in mod_label_mapping
                else mod_label_mapping[mods[i]]
            )
            seq[s] = f"{seq[s]}{mod}"
        seq = "".join(seq)
    else:
        seq = prec.sequence
    return f"{seq} (+{prec.charge})"


def add_corr(df_obs, df_library):
    df_merged = df_library.merge(df_obs, on="mz")
    df_merged.columns = [
        c.replace("_x", "_pred").replace("_y", "_obs") for c in df_merged.columns
    ]

    df_merged = df_merged.query("intensity_pred!=0 & intensity_obs!=0")
    df_corr = (
        df_merged.groupby("rt_idx")[["intensity_pred", "intensity_obs"]]
        .corr()
        .reset_index()
        .melt(id_vars=["rt_idx", "level_1"])
        .query("value!=1")
        .query('level_1=="intensity_obs"')
        .rename(columns={"value": "corr_coeff"})[["corr_coeff", "rt_idx"]]
    )
    df_obs_corr = df_obs.merge(df_corr, on="rt_idx")

    df_obs_corr["rt_cat"] = df_obs_corr.apply(
        lambda x: f"{x.rt_idx}, r={round(x.corr_coeff, 2)}", axis=1
    )
    return df_obs_corr


def plot_mirror_byRT(
    dense,
    mz_library,
    intensity_library,
    label_library,
    precursor_entry,
    add_corr_coeff=True,
    width=800,
    height=300,
):
    seq = format_prec_entry(precursor_entry)
    df_library = convert_library_to_df(mz_library, intensity_library)
    df_obs = split_dense_byRT(dense, mz_library, label_library)

    if add_corr_coeff:
        df_obs = add_corr(df_obs, df_library)

    return plot_mirror_byRT_from_dfs(df_obs, df_library, seq, width, height)


def plot_mirror_byRT_from_dfs(df_obs, df_library, title="", width=800, height=300):
    theo = plot_theo(df_library)
    obs = plot_obs_byRT(df_obs)

    # middle line to separate obs vs theo
    middle_line = (
        altair.Chart(pd.DataFrame({"sep": [0]}))
        .mark_rule(size=3)
        .encode(y="sep", color=altair.value("lightGray"))
    )

    ## title
    title = altair.TitleParams(
        text=title,
        # fontSize=16,
        fontWeight="bold",
        # anchor="start",  # Align title to the left
        color="black",
    )

    return (obs + theo + middle_line).properties(
        width=width, height=height, title=title
    )


def plot_obs_byRT(df_plot):
    annotation_kws = {"align": "left", "angle": 270, "baseline": "middle"}

    anno = [
        altair.Tooltip("mz", format=".3f", title="m/z"),
        altair.Tooltip(
            "norm_intensity", format=".2f", title="Intensity"
        ),  # format='.%'
    ]

    df_plot["label_color"] = df_plot["frag_label"].apply(
        lambda x: fragment_colos[x[0]]
        if (x != "" and x[0]) in fragment_colos
        else "#0000FF"
    )

    num_rts = df_plot[["rt_cat"]].drop_duplicates().shape[0]
    rt_colors = _get_colors(n_colors=num_rts)  # get_palette_for_RT(num_rts)
    color = altair.Color(
        "rt_cat", scale=altair.Scale(scheme="set2", range=rt_colors), title="RT Scan"
    )
    x = altair.X(
        "mz",
        axis=altair.Axis(title="m/z", titleFontStyle="italic", grid=True),
        scale=altair.Scale(
            nice=False,
            padding=5,
            zero=False,
            domain=[50, (max(df_plot["mz"]) // 10) * 10 + 50],
        ),
    )

    y = altair.Y(
        "sum(norm_intensity):Q",
        axis=altair.Axis(title=["Intensity"], format=".2f", grid=True),  # format='%'
        # scale=altair.Scale(nice=True, padding=0)
    )
    color2 = altair.Color(
        "label_color", scale=None
    )  # altair.Scale(scheme='set2'), title='RT Scan')
    frag_anno = (
        altair.Chart(df_plot)
        .mark_text(dx=10, **annotation_kws)
        .encode(x=x, y=y, text="frag_label", color=color2)
    )

    return frag_anno + (
        altair.Chart(df_plot)
        .mark_bar(size=3)
        .encode(x=x, y=y, color=color, tooltip=anno)
    )


def plot_theo(df_plot):
    anno = [
        altair.Tooltip("mz", format=".3f", title="m/z"),
        altair.Tooltip("intensity", format=".2f", title="Intensity"),  # format='.%'
    ]

    df_plot["minus_intensity"] = -df_plot["norm_intensity"]
    x = altair.X(
        "mz",
        axis=altair.Axis(title="m/z", titleFontStyle="italic", grid=True),
        scale=altair.Scale(
            nice=False,
            padding=5,
            zero=False,
            domain=[50, (max(df_plot["mz"]) // 10) * 10 + 50],
        ),
    )

    y = altair.Y(
        "minus_intensity",
        axis=altair.Axis(title=["Intensity"], format="", grid=True),  # format='%'
        scale=altair.Scale(nice=True, padding=0),
    )
    return altair.Chart(df_plot).mark_rule(size=3).encode(x=x, y=y, tooltip=anno)


def plot_xic_w_background(spectrum_slice, palette_name=None, hex_colors=None):
    xic_observed = spectrum_slice[0].sum(
        axis=(1, 2)
    )  ## this gives us mz info and rt info
    ## this gives intensities for every m/z and RT scan
    df_xic = (
        pd.DataFrame(xic_observed).reset_index().rename(columns={"index": "mz_scan"})
    )
    df_xic = df_xic.melt(
        id_vars=["mz_scan"], var_name="RT_scan", value_name="intensity"
    )
    n_colors = df_xic[["RT_scan"]].drop_duplicates().shape[0]

    # Compute mean/median intensity across mz scans
    median_data = df_xic.groupby("RT_scan", as_index=False)["intensity"].mean()

    try:  # if hex_colors is not None:
        background_colors = hex_colors[:n_colors]
    except:
        background_colors = _get_colors(palette_name=palette_name, n_colors=n_colors)

    # Create background bars to distinguish RT scans
    background = (
        altair.Chart(df_xic)
        .mark_bar(opacity=0.05)
        .encode(
            x=altair.X("RT_scan:O", axis=altair.Axis(title="RT scan")),
            # x=altair.X('RT_scan:O', axis=None),
            y=altair.value(1),  # Dummy value for full-height bars
            color=altair.Color(
                "RT_scan:O",
                scale=altair.Scale(
                    domain=list(df_xic["RT_scan"].unique()), range=background_colors
                ),
                legend=None,
            ),
        )
        .properties(width=400, height=300, title="")
    )

    # Create the line plot with a continuous color legend
    base_chart = (
        altair.Chart(df_xic)
        .mark_line()
        .encode(
            x=altair.X(
                "RT_scan:Q", axis=altair.Axis(title=None, labels=False, ticks=False)
            ),  #'RT_scan:Q',
            y="intensity:Q",
            color=altair.Color(
                "mz_scan:Q",
                scale=altair.Scale(scheme="greys"),
                legend=altair.Legend(title="Fragment No."),
            ),
        )
    )

    # Add median intensity line in red with legend
    median_chart = (
        altair.Chart(median_data)
        .mark_line(color="red")
        .encode(
            x=altair.X(
                "RT_scan:Q", axis=altair.Axis(title=None, labels=False, ticks=False)
            ),  #'RT_scan:Q',
            # x='RT_scan:Q',
            y="intensity:Q",
        )
    )

    # Combine charts
    chart = background + base_chart + median_chart

    return chart


def plot_correlations(spectrum_slice, mz_library):
    intensity_slice = spectrum_slice[0].sum(axis=1).sum(axis=1)
    normalized_intensity_slice = normalize_profiles(intensity_slice)
    median_profile = median_axis(normalized_intensity_slice, axis=0)
    corr_list = correlation_coefficient(median_profile, intensity_slice)

    df_corrs = pd.DataFrame({"mz": mz_library, "corr_coeff": corr_list})

    # Create the bar plot
    corr_p = (
        altair.Chart(df_corrs)
        .mark_circle()
        .encode(
            x=altair.X("mz:Q", axis=altair.Axis(title="m/z")),
            y=altair.X("corr_coeff:Q", axis=altair.Axis(title="Correlation")),
            # size='corr_coeff:Q'  # Circle size based on the correlation coefficient
            size=altair.Size(
                "corr_coeff:Q",  # scale=altair.Scale(scheme='greys'),
                legend=altair.Legend(title="Correlation"),
            ),
        )
    )
    return corr_p


def mirror_w_xic_w_corrs(
    spectrum_slice,
    mz_library,
    intensity_library,
    fragment_library,
    precursor_df,
    selected_hash,
    width=600,
    height=300,
):
    precursor_entry = precursor_df[
        precursor_df["mod_seq_charge_hash"] == selected_hash
    ].iloc[0]

    mirror = plot_mirror_byRT(
        spectrum_slice,
        mz_library,
        intensity_library,
        fragment_library,
        precursor_entry,
        width=width * 0.8,
        height=height,
    )
    xic = plot_xic_w_background(
        spectrum_slice
    )  # .properties(width=width*0.2, height=height*0.3)

    corrs = plot_correlations(spectrum_slice, mz_library)

    mirror_corrs = altair.vconcat(
        mirror.properties(width=width * 0.75, height=height * 0.8),
        corrs.properties(width=width * 0.75, height=height * 0.2),
    ).resolve_scale(x="shared", color="independent", size="independent")

    p = (
        (mirror_corrs | xic.properties(width=width * 0.25, height=height * 0.5))
        .resolve_scale(color="independent")
        .configure_view(stroke=None)
    )

    return p
