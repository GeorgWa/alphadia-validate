import os
from pathlib import Path

import pandas as pd
from alphabase.spectral_library.base import SpecLibBase

from alphadia.data.alpharaw_wrapper import Thermo
from alphadia.test_data_downloader import DataShareDownloader
from alphabase.spectral_library.flat import SpecLibFlat

# Bulk injections of HeLa cell lysate acquired on the Orbitrap Astral
RAW_DATA_URL_LIST = (
    "https://datashare.biochem.mpg.de/s/339jg5HtGrwLwDN/download?files=20231017_OA2_TiHe_ADIAMA_HeLa_200ng_Evo011_21min_F-40_07.raw",
)

# results from search_1.10.0.ipynb
PRECURSORS_TSV_URL = "https://datashare.biochem.mpg.de/s/fsCqlT757ttVWI8"
SPECLIB_URL = "https://datashare.biochem.mpg.de/s/VxakpS6mhM2IwxJ"


from joblib import Memory

memory = Memory("cachedir")


@memory.cache
def prepare_data(
    main_folder: Path,
) -> tuple[pd.DataFrame, SpecLibBase, SpecLibFlat, Thermo]:
    """Prepare raw & results data.

    If the main_folder/output does not already contain results (precursors.tsv and speclib.hdf) from
    a previous search, this function will download the raw data and return the dataframes and objects

     Note: this is just required if you did not run the search yourself in search_1.10.0.ipynb

    :param main_folder: folder in which data is expected / will be downloaded to
    :return:
    """
    precursor_df, current_raw_path, speclib_path = _download_data(main_folder)

    spectral_library = SpecLibBase()
    spectral_library.load_hdf(speclib_path)

    dia_data = Thermo(current_raw_path)

    spectral_library_flat = SpecLibFlat()
    spectral_library_flat.parse_base_library(spectral_library)

    return precursor_df, spectral_library, spectral_library_flat, dia_data


def _download_data(main_folder: Path) -> tuple[pd.DataFrame, str, str]:
    """Download data if not already present."""
    data_folder = main_folder / "data"
    output_folder = main_folder / "output"
    for folder in [data_folder, output_folder]:
        if not os.path.exists(folder):
            os.makedirs(folder)
    raw_file_paths = [
        DataShareDownloader(url, str(data_folder)).download()
        for url in RAW_DATA_URL_LIST
    ]
    precursors_tsv_path = DataShareDownloader(
        PRECURSORS_TSV_URL, str(output_folder)
    ).download()
    speclib_path = DataShareDownloader(SPECLIB_URL, str(output_folder)).download()

    current_raw_path = raw_file_paths[0]
    current_raw_name = os.path.basename(current_raw_path).replace(".raw", "")
    precursor_df = pd.read_csv(precursors_tsv_path, sep="\t")
    precursor_df = precursor_df[precursor_df["run"] == current_raw_name]

    return precursor_df, current_raw_path, speclib_path


def display_spectral_library(spectral_library):
    print("precursor_df")
    display(
        spectral_library.precursor_df[
            [
                "precursor_mz",
                "sequence",
                "mods",
                "mod_sites",
                "charge",
                "mod_seq_charge_hash",
                "frag_start_idx",
                "frag_stop_idx",
            ]
        ].head()
    )

    print("fragment_mz_df")
    display(spectral_library.fragment_mz_df.head())

    print("fragment_intensity_df")
    display(spectral_library.fragment_intensity_df.head())

    from alphabase.spectral_library.flat import SpecLibFlat

    spectral_library_flat = SpecLibFlat()
    spectral_library_flat.parse_base_library(spectral_library)

    print("spectral_library_flat.fragment_df")
    display(spectral_library_flat.fragment_df)
