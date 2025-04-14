import os
from pathlib import Path

import pandas as pd
from alphabase.spectral_library.base import SpecLibBase

from alphadia.data.alpharaw_wrapper import Thermo, Sciex, MzML
from alphadia.data.bruker import TimsTOFTranspose

from alphabase.tools.data_downloader import DataShareDownloader
from alphabase.spectral_library.flat import SpecLibFlat

# Bulk injections of HeLa cell lysate acquired on the Orbitrap Astral
RAW_DATA_URL = "https://datashare.biochem.mpg.de/s/VfqtW5p9MJ0kxAC/download?files=20231017_OA2_TiHe_ADIAMA_HeLa_200ng_Evo011_21min_F-40_07.mzML"


# results from search_1.10.0.ipynb
PRECURSORS_TSV_URL = "https://datashare.biochem.mpg.de/s/VfqtW5p9MJ0kxAC/download?files=precursors.tsv"
SPECLIB_URL = "https://datashare.biochem.mpg.de/s/VfqtW5p9MJ0kxAC/download?files=speclib.hdf"


from joblib import Memory

memory = Memory(".memoize_cachedir")


@memory.cache
def prepare_data(
    main_folder: Path,
    *,
    download_data: bool = False,
    raw_file_name: str = None,
    precursors_file_name: str = None,
    speclib_file_name: str = None,
) -> tuple[pd.DataFrame, SpecLibBase, SpecLibFlat, Thermo]:
    """Prepare raw & results data and return data objects.

    If the download_data is true, this function will download the required data.

    :param main_folder: folder in which data is expected / will be downloaded to
    :param download_data: whether to download the data
    :param precursors_file_name: name of the precursors file (required if download_data is False)
    :param raw_file_name: name of the raw file (required if download_data is False)
    :param speclib_file_name: name of the speclib file (required if download_data is False)
    :return:
    """
    if download_data:
        precursors_tsv_path, raw_file_path, speclib_path = _download_data(main_folder)
    else:
        if (
            precursors_file_name is None
            or raw_file_name is None
            or speclib_file_name is None
        ):
            raise ValueError(
                "Please provide the file names for the precursors, raw file and speclib files"
            )

        precursors_tsv_path = main_folder / precursors_file_name
        raw_file_path = main_folder / raw_file_name
        speclib_path = main_folder / speclib_file_name

    current_raw_name = raw_file_path.stem
    precursor_df = pd.read_csv(precursors_tsv_path, sep="\t")
    precursor_df = precursor_df[precursor_df["run"] == current_raw_name]

    spectral_library = SpecLibBase()
    spectral_library.load_hdf(speclib_path)

    if raw_file_path.suffix.lower() == ".mzml":
        dia_data = MzML(str(raw_file_path))
    elif raw_file_path.suffix.lower() == ".raw":
        dia_data = Thermo(str(raw_file_path))
    elif raw_file_path.suffix.lower() == ".wiff":
        dia_data = Sciex(str(raw_file_path))
    elif raw_file_path.suffix.lower() == ".d":
        dia_data = TimsTOFTranspose(str(raw_file_path))
    else:
        raise ValueError(
            f"Unsupported file type: {raw_file_path.suffix}. Supported types are .mzML, .raw, .wiff, and .d"
        )

    spectral_library_flat = SpecLibFlat()
    spectral_library_flat.parse_base_library(spectral_library)

    return precursor_df, spectral_library, spectral_library_flat, dia_data


def _download_data(main_folder: Path) -> tuple[Path, Path, Path]:
    """Download data if not already present."""
    os.makedirs(main_folder, exist_ok=True)

    raw_file_path = DataShareDownloader(RAW_DATA_URL, str(main_folder)).download()

    precursors_tsv_path = DataShareDownloader(
        PRECURSORS_TSV_URL, str(main_folder)
    ).download()
    speclib_path = DataShareDownloader(SPECLIB_URL, str(main_folder)).download()

    return Path(precursors_tsv_path), Path(raw_file_path), Path(speclib_path)


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
