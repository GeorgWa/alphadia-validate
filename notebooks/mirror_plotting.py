import altair
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import math

fragment_colos = {'a': '#388E3C', 'b': '#1976D2', 'c': '#00796B',
          'x': '#7B1FA2', 'y': '#D32F2F', 'z': '#F57C00',
          'unknown': '#212121', None: "#808080" }

def split_dense_byRT(dense, mz_library, label_library):
    intensity_observed_full = dense[0].sum(axis=(1, 2)) ## do not sum over rt have shape (int, rt)
    
    intensity_flat = intensity_observed_full.ravel()
    mz_flat = np.tile(mz_library, dense[0].shape[3])
    #intensity_flat.shape, mz_flat.shape

    # Recover the original (mz, rt) index from the flattened array
    # Flatten mz library
    reconstructed_indices = np.unravel_index(np.arange(intensity_flat.size), intensity_observed_full.shape)
    mz_flat = mz_library[reconstructed_indices[0]]
    label_flat = label_library[reconstructed_indices[0]]

    indices = np.indices(intensity_observed_full.shape)
    rt_indices = indices[1].ravel()      # Retention time indices
    #cylce_indices = ...

    df_obs = pd.DataFrame({'mz':mz_flat, 
                           'frag_label':label_flat,
                        'intensity':intensity_flat,
                        'rt_idx': rt_indices,
                       # 'cycle_idx': cycle_indices
                       } )
    df_obs['rt_cat'] = pd.Categorical(df_obs.rt_idx)
    max_int = df_obs.groupby(['mz'])['intensity'].sum().reset_index()['intensity'].max()
    df_obs['norm_intensity'] = df_obs['intensity']/max_int

    df_obs['frag_label'] = df_obs.apply(lambda x: '' if x.intensity<=0 else x.frag_label, axis=1)
    
    return df_obs

def convert_library_to_df(mz_library, intensity_library):
    df_library = pd.DataFrame({'mz': mz_library, 'intensity': intensity_library})
    df_library['norm_intensity'] = df_library['intensity']/ max(df_library['intensity'])
    return df_library

def interpolate_colors_green_to_red(num_colors):
    colors = []
    for i in range(num_colors):
        red = int(255 * (i / (num_colors - 1))**0.7)  # Interpolating red channel
        green = int(255 * (1 - (i / (num_colors - 1))**0.3))  # Interpolating green channel
        colors.append(f"#{red:02X}{green:02X}00")  # Keep blue at 0
    return colors

# Function to interpolate colors from red to yellow (existing function)
def interpolate_colors_red_to_yellow(num_colors):
    colors = []
    for i in range(num_colors):
        green = int(255 * (i / (num_colors ))**0.3)  # Linear interpolation for the green channel
        colors.append(f"#FF{green:02X}00")  # Format as HEX, keeping red at 255 and blue at 0
    return colors

def get_palette_for_RT(num_rt_scans):
    if num_rt_scans ==1:
        return ['#FF0000']
    left_cols = interpolate_colors_green_to_red(math.floor(num_rt_scans/2.)+1)[:-1]
    right_cols = interpolate_colors_red_to_yellow(math.ceil(num_rt_scans/2.))
    return left_cols + right_cols

def format_prec_entry(prec):
    if type(prec.mods)=='str':
        mod_label_mapping = {'Carbamidomethyl':'', 'Oxidation':'(Ox)'}
        prec = precursor_df.iloc[9]
        seq = list(prec.sequence)
        mods = [m.split('@')[0] for m in prec.mods.split(';')]
        sites = np.array(prec.mod_sites.split(';'), dtype=int) -1 
        for i, s in enumerate(sites):
            mod = mods[i] if mods[i] not in mod_label_mapping else mod_label_mapping[mods[i]]
            seq[s] = f'{seq[s]}{mod}'
        seq = ''.join(seq) 
    else:
        seq = prec.sequence
    return f'{seq} (+{prec.charge})'

def plot_mirror_byRT(dense, mz_library, intensity_library, label_library, 
                     precursor_entry,
                     width=800, height=300):
    seq = format_prec_entry(precursor_entry)
    df_library = convert_library_to_df(mz_library, intensity_library)
    df_obs = split_dense_byRT(dense, mz_library, label_library)
    
    return plot_mirror_byRT_from_dfs(df_obs, df_library, seq, width, height)
    
def plot_mirror_byRT_from_dfs(df_obs, df_library, title='', width=800, height=300):
    theo = plot_theo(df_library)
    obs = plot_obs_byRT(df_obs)

    # middle line to separate obs vs theo
    middle_line = (altair.Chart(pd.DataFrame({'sep': [0]}))
        .mark_rule(size=3)
        .encode( y='sep',  color=altair.value('lightGray'))
    )

    ## title
    title = altair.TitleParams(
        text=title,
        #fontSize=16,
        fontWeight="bold",
        #anchor="start",  # Align title to the left
        color="black"
    )
    
    return (obs + theo + middle_line).properties(width=width, height=height, title=title)
    
def plot_obs_byRT(df_plot):
    annotation_kws = {'align': 'left',
                          'angle': 270, 
                          'baseline': 'middle'}
    
    anno = [ altair.Tooltip('mz', format='.3f', title='m/z'),
             altair.Tooltip('norm_intensity', format='.2f', title='Intensity') #format='.%'
           ] 
    
    df_plot['label_color'] = (df_plot['frag_label']
        .apply(lambda x: fragment_colos[x[0]] if (x!='' and x[0]) in fragment_colos else '#0000FF')
    )
    
    num_rts = df_plot[['rt_cat']].drop_duplicates().shape[0]
    rt_colors = get_palette_for_RT(num_rts)
    color = altair.Color('rt_cat', scale=altair.Scale(scheme='set2', range=rt_colors), title='RT Scan')
    x = altair.X('mz', axis=altair.Axis(title='m/z', titleFontStyle='italic', grid=True),
                              scale=altair.Scale(nice=False, padding=5, zero=False,
                                                 domain=[50, (max(df_plot['mz'])//10)*10 + 50] 
                                                ) ) 
            
    y = altair.Y('sum(norm_intensity):Q', 
                 axis=altair.Axis(title=['Intensity'], format='.2f', grid=True), #format='%'
                              #scale=altair.Scale(nice=True, padding=0)
                )
    color2 = altair.Color('label_color', scale=None)#altair.Scale(scheme='set2'), title='RT Scan')
    frag_anno = (altair.Chart(df_plot)
                      .mark_text(dx=10, **annotation_kws)
                      .encode(x=x, 
                              y=y, 
                              text='frag_label',
                              color=color2
                              )
                )
    
    return frag_anno + (altair.Chart(df_plot)
             .mark_bar(size=3)
             .encode(x=x, y=y, color=color, tooltip=anno)
    )

def plot_theo(df_plot):
    anno = [ altair.Tooltip('mz', format='.3f', title='m/z'),
             altair.Tooltip('intensity', format='.2f', title='Intensity') #format='.%'
           ] 

    df_plot['minus_intensity'] = -df_plot['norm_intensity']
    x = altair.X('mz', axis=altair.Axis(title='m/z', titleFontStyle='italic', grid=True),
                              scale=altair.Scale(nice=False, padding=5, zero=False,
                                                 domain=[50, (max(df_plot['mz'])//10)*10 + 50] 
                                                ) ) 
            
    y = altair.Y('minus_intensity', axis=altair.Axis(title=['Intensity'], format='', grid=True), #format='%'
                              scale=altair.Scale(nice=True, padding=0))
    return (altair.Chart(df_plot)
             .mark_rule(size=3)
             .encode(x=x, y=y, tooltip=anno)
    )


