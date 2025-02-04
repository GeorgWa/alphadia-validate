import altair
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def split_dense_byRT(dense, mz_library):
    intensity_observed_full = dense[0].sum(axis=(1, 2)) ## do not sum over rt have shape (int, rt)
    
    intensity_flat = intensity_observed_full.ravel()
    mz_flat = np.tile(mz_library, dense[0].shape[3])
    #intensity_flat.shape, mz_flat.shape

    # Recover the original (mz, rt) index from the flattened array
    # Flatten mz library
    reconstructed_indices = np.unravel_index(np.arange(intensity_flat.size), intensity_observed_full.shape)
    mz_flat = mz_library[reconstructed_indices[0]]

    indices = np.indices(intensity_observed_full.shape)
    rt_indices = indices[1].ravel()      # Retention time indices
    #cylce_indices = ...

    df_obs = pd.DataFrame({'mz':mz_flat, 
                        'intensity':intensity_flat,
                        'rt_idx': rt_indices,
                       # 'cycle_idx': cycle_indices
                       } )
    df_obs['rt_cat'] = pd.Categorical(df_obs.rt_idx)
    max_int = df_obs.groupby(['mz'])['intensity'].sum().reset_index()['intensity'].max()
    df_obs['norm_intensity'] = df_obs['intensity']/max_int
    return df_obs

def convert_library_to_df(mz_library, intensity_library):
    df_library = pd.DataFrame({'mz': mz_library, 'intensity': intensity_library})
    df_library['norm_intensity'] = df_library['intensity']/ max(df_library['intensity'])
    return df_library
    
def plot_mirror_byRT(dense, mz_library, intensity_library,  width=800, height=300):
    df_library = convert_library_to_df(mz_library, intensity_library)
    df_obs = split_dense_byRT(dense, mz_library)
    return plot_mirror_byRT_from_dfs(df_obs, df_library, width, height)
    
    
def plot_mirror_byRT_from_dfs(df_obs, df_library, width=800, height=300):
    theo = plot_theo(df_library)
    obs = plot_obs_byRT(df_obs)

    # middle line to separate obs vs theo
    middle_line = (altair.Chart(pd.DataFrame({'sep': [0]})).mark_rule(size=3).encode( y='sep',  color=altair.value('lightGray')))

    return (obs + theo + middle_line).properties(width=width, height=height)
    
def plot_obs_byRT(df_plot):
    anno = [ altair.Tooltip('mz', format='.3f', title='m/z'),
             altair.Tooltip('norm_intensity', format='.2f', title='Intensity') #format='.%'
           ] 
    
    color = altair.Color('rt_cat', scale=altair.Scale(scheme='set2'), title='RT Scan')
    x = altair.X('mz', axis=altair.Axis(title='m/z', titleFontStyle='italic', grid=True),
                              scale=altair.Scale(nice=False, padding=5, zero=False,
                                                 domain=[50, (max(df_plot['mz'])//10)*10 + 50] 
                                                ) ) 
            
    y = altair.Y('sum(norm_intensity):Q', 
                 axis=altair.Axis(title=['Intensity'], format='.2f', grid=True), #format='%'
                              #scale=altair.Scale(nice=True, padding=0)
                )
    return (altair.Chart(df_plot)
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


