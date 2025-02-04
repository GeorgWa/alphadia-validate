import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .slicer import SpectrumSlicer, get_flat_library_entry_by_hash

def isoform_attention_plot(spectral_library_flat, precursor_df, dia_data, selected_hash):

    slicer = SpectrumSlicer(spectral_library_flat, precursor_df, dia_data)
    
    precursor_entry = precursor_df[precursor_df['mod_seq_charge_hash'] == selected_hash].iloc[0]

    mz_library, intensity_library, data_slice = slicer.get_by_hash(selected_hash)

    library_isoforms = spectral_library_flat.precursor_df.loc[
        np.logical_and(spectral_library_flat.precursor_df['sequence'] == precursor_entry['sequence'],
                       spectral_library_flat.precursor_df['mods'] == precursor_entry['mods']), :]

    fragment_tables = []
    for i, e in library_isoforms.iterrows():
        t = spectral_library_flat.fragment_df[e['flat_frag_start_idx']:e['flat_frag_stop_idx']]
        t['isoform'] = e['mod_sites']
        fragment_tables.append(t)

    spec = pd.concat(fragment_tables).groupby('mz').agg({'isoform': lambda s: '+'.join(s)})
    spec['discrimination'] = 1 - (spec['isoform'].str.count('\+')) / library_isoforms.shape[0]

    library_spec = pd.DataFrame({'mz': mz_library, 'intensity': intensity_library})
    library_spec['isoform'] = list(spec.reindex(library_spec['mz'])['isoform'])

    intensity_observed = data_slice[0].sum(axis=(1,2,3))
    intensity_observed_normalized = intensity_observed / intensity_observed.max()
    library_spec['intensity'] = library_spec['intensity'] / library_spec['intensity'].max()

    plt.stem(mz_library, intensity_observed_normalized)
    for n, (iii, gr) in enumerate(library_spec.groupby('isoform')):
        plt.stem(gr['mz'], -1 * gr['intensity'], linefmt=plt.colormaps['tab10'](n + 1), label=iii)

    plt.legend()
    plt.title('{} {} {}'.format(precursor_entry['sequence'], precursor_entry['mods'], precursor_entry['mod_sites']))
    plt.show()
