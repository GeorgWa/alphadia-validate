from xic_utils import normalize_profiles, median_axis, correlation_coefficient
import altair as alt
import pandas as pd
import numpy as np

def plot_advanced_xic(spectrum_slice, dia_data):

    


    intensity_slice = spectrum_slice[0].sum(axis=1).sum(axis=1)
    normalized_intensity_slice = normalize_profiles(intensity_slice)
    median_profile = median_axis(normalized_intensity_slice, axis=0)
    corr_list = correlation_coefficient(median_profile, intensity_slice)
  
  
    # Create a DataFrame for plotting
    df_fragments = pd.DataFrame()
    for fragment in range(normalized_intensity_slice.shape[0]):
        temp_df = pd.DataFrame({
            'Retention Time': range(len(normalized_intensity_slice[fragment])),
            'Intensity': normalized_intensity_slice[fragment],
            'Fragment': f'Fragment {fragment} (r={corr_list[fragment]:.2f})',
            'Type': 'Fragment Profile'
        })
        df_fragments = pd.concat([df_fragments, temp_df])
    
    df_median = pd.DataFrame({
        'Retention Time': range(len(median_profile)),
        'Intensity': median_profile,
        'Type': 'Median Profile'
    })

    # Create fragment profiles layer with different colors
    fragment_layer = alt.Chart(df_fragments).mark_line(
        opacity=0.8
    ).encode(
        x='Retention Time:Q',
        y='Intensity:Q',
        color='Fragment:N',
        detail='Fragment:N'
    )

    # Create median profile layer
    median_layer = alt.Chart(df_median).mark_line(
        color='black',
        size=4
    ).encode(
        x='Retention Time:Q',
        y='Intensity:Q'
    )

    # Combine layers and configure chart
    chart = (fragment_layer + median_layer).properties(
        width=400,
        height=400,
        title='Fragment Intensity Profiles'
    ).configure_axis(
        titleFontSize=12
    ).configure_legend(
        titleFontSize=12,
        labelFontSize=10
    )

    return chart

