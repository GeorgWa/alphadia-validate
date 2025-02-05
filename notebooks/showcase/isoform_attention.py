import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from alphabase.protein.fasta import SpecLibFasta
from alphabase.spectral_library.flat import SpecLibFlat
from .slicer import SpectrumSlicer, get_flat_library_entry_by_hash

def fragment_isoform_annotation(precursor_entry, fixed_mods = {'Carbamidomethyl@C'}):    
    #parsing precursor entry
    sequence = precursor_entry['sequence']
    charge = precursor_entry['charge']

    if precursor_entry['mods'] == '':
        mod_counter = Counter()
    else:
        mod_counter = Counter(precursor_entry['mods'].split(';'))
    
    var_mod_count = [mod_counter[mod] for mod in mod_counter.keys() if mod not in fixed_mods]
    if len(var_mod_count) == 0:
        max_var_mod_num = 0
    else:
        max_var_mod_num = max(var_mod_count)


    potential_fragements = [f'{type}_z{z}' for type in ['b', 'y'] for z in range(1, charge + 1)]

    #create all possible isoforms
    precursor_df = pd.DataFrame({
        'sequence': [sequence],
        'charge': [charge],
    })

    fastalib = SpecLibFasta(
        potential_fragements, 
        var_mods=[mod for mod in mod_counter.keys() if mod not in fixed_mods],
        fix_mods=[mod for mod in mod_counter.keys() if mod in fixed_mods],
        min_var_mod_num=max_var_mod_num,
        max_var_mod_num=max_var_mod_num
    )

    fastalib.precursor_df = precursor_df
    fastalib.add_modifications()
    fastalib.calc_fragment_mz_df()
    flatlib = SpecLibFlat()
    flatlib.parse_base_library(fastalib)

    #TODO look for better way to annotate isoforms
    fragment_tables = []
    for _, e in flatlib.precursor_df.iterrows():
        t = flatlib.fragment_df[e['flat_frag_start_idx']:e['flat_frag_stop_idx']]
        t['isoform'] = e['mod_sites']
        fragment_tables.append(t)

    #TODO better grouping?
    isoform_map = pd.concat(fragment_tables).groupby('mz').agg({'isoform': lambda s: '+'.join(s.unique())})

    return isoform_map

def isoform_attention_plot(spectral_library_flat, precursor_df, dia_data, selected_hash):

    slicer = SpectrumSlicer(spectral_library_flat, precursor_df, dia_data)
    
    precursor_entry = precursor_df[precursor_df['mod_seq_charge_hash'] == selected_hash].iloc[0]

    mz_library, intensity_library, data_slice = slicer.get_by_hash(selected_hash)

    isoform_map = fragment_isoform_annotation(precursor_entry)
    
    library_spec = pd.DataFrame({'mz': mz_library, 'intensity': intensity_library})
    library_spec['isoform'] = list(isoform_map.reindex(library_spec['mz'])['isoform'])

    intensity_observed = data_slice[0].sum(axis=(1,2,3))
    intensity_observed_normalized = intensity_observed / intensity_observed.max()
    library_spec['intensity'] = library_spec['intensity'] / library_spec['intensity'].max()

    plt.stem(mz_library, intensity_observed_normalized)
    for n, (iii, gr) in enumerate(library_spec.groupby('isoform')):
        plt.stem(gr['mz'], -1 * gr['intensity'], linefmt=plt.colormaps['tab10'](n + 1), label=iii)

    plt.legend()
    plt.title('{} {} {}'.format(precursor_entry['sequence'], precursor_entry['mods'], precursor_entry['mod_sites']))
    plt.show()
