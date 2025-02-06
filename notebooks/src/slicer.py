import numpy as np
import pandas as pd
from alphabase.spectral_library.flat import SpecLibFlat


def get_library_entry_by_hash(speclib, hash_, min_intensity=0.01):
    speclib_entry = speclib.precursor_df[
        speclib.precursor_df["mod_seq_charge_hash"] == hash_
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


def get_flat_library_entry_by_hash(speclib_flat, hash_, min_intensity=0.01):
    speclib_entry = speclib_flat.precursor_df[
        speclib_flat.precursor_df["mod_seq_charge_hash"] == hash_
    ].iloc[0]

    flat_frag_start_idx = speclib_entry.flat_frag_start_idx
    flat_frag_stop_idx = speclib_entry.flat_frag_stop_idx

    fragment_mz = (
        speclib_flat.fragment_df["mz"]
        .iloc[flat_frag_start_idx:flat_frag_stop_idx]
        .to_numpy()
        .flatten()
    )
    fragment_intensity = (
        speclib_flat.fragment_df["intensity"]
        .iloc[flat_frag_start_idx:flat_frag_stop_idx]
        .to_numpy()
        .flatten()
    )

    fragment_label = (
        speclib_flat.fragment_df["fragment_label"]
        .iloc[flat_frag_start_idx:flat_frag_stop_idx]
        .to_numpy()
        .flatten()
    )

    fragment_mask = fragment_intensity > min_intensity

    fragment_mz = fragment_mz[fragment_mask]
    fragment_intensity = fragment_intensity[fragment_mask]
    fragment_label = fragment_label[fragment_mask]

    # sort both by mz
    fragment_order = np.argsort(fragment_mz)
    fragment_mz = fragment_mz[fragment_order]
    fragment_intensity = fragment_intensity[fragment_order]
    fragment_label = fragment_label[fragment_order]

    return speclib_entry, fragment_mz, fragment_intensity, fragment_label


type_map = {97: "a", 98: "b", 99: "c", 120: "x", 121: "y", 12: "z"}
loss_map = {0: "", 18: "H2O", 17: "NH3"}


def get_ion_labels(row):
    ls_label = loss_map[int(row["loss_type"])]
    ls_label = "_" + ls_label if ls_label != "" else ""
    label = (
        f"{type_map[row['type']]}{int(row['number'])}{ls_label}(+{int(row['charge'])})"
    )
    return label


class SpectrumSlicer:
    def __init__(
        self,
        spectral_library_flat: SpecLibFlat,
        precursor_df: pd.DataFrame,
        dia_data,
    ):
        # Get fragment annotations
        spectral_library_flat.fragment_df["fragment_label"] = (
            spectral_library_flat.fragment_df.apply(get_ion_labels, axis=1)
        )
        self.spectral_library_flat = spectral_library_flat

        self.precursor_df = precursor_df
        self.dia_data = dia_data

    def get_by_hash(
        self, selected_hash: int
    ) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """#### 4 Visualize precursor data

        :returns: np.ndarray
        The spectrum data `dense` is a 5 dimensional numpy array with a dense slice of the spectrum space.
        The dimensions are:
        - 0: either intensity information 0 or relative mass error 1
        - 1: index of the fragment mz which was queried
        - 2: ion mobility dimension (will be zero for DIA data)
        - 3: The observations in the DIA cycle. As there might be multiple quadrupole windows where the precursor was detected, this will be a list of observations.
        - 4: Retention time datapoints.
        """

        speclib_entry, mz_library, intensity_library, fragment_library = (
            get_flat_library_entry_by_hash(self.spectral_library_flat, selected_hash)
        )
        precursor_entry = self.precursor_df[
            self.precursor_df["mod_seq_charge_hash"] == selected_hash
        ].iloc[0]

        jit_data = self.dia_data.jitclass()
        precursor_query = np.array(
            [[speclib_entry.precursor_mz, speclib_entry.precursor_mz]], dtype=np.float32
        )
        scan_limits = np.array(
            [[precursor_entry.scan_start, precursor_entry.scan_stop, 1]], dtype=np.int64
        )
        frame_limits = np.array(
            [[precursor_entry.frame_start, precursor_entry.frame_stop, 1]],
            dtype=np.int64,
        )

        spectrum_slice, _precursor_index = jit_data.get_dense(
            frame_limits,
            scan_limits,
            mz_library,
            30,
            precursor_query,
        )

        return mz_library, intensity_library, spectrum_slice, fragment_library
