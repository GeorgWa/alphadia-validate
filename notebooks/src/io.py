import os
from pathlib import Path

import pandas as pd
from alphabase.spectral_library.base import SpecLibBase
from alphabase.spectral_library.flat import SpecLibFlat
from alphabase.tools.data_downloader import DataShareDownloader
from alphadia.raw_data.alpharaw_wrapper import MzML, Sciex, Thermo
from alphadia.raw_data.bruker import TimsTOFTranspose

# Bulk injections of HeLa cell lysate acquired on the Orbitrap Astral
EXAMPLE_RAW_DATA_URL = "https://datashare.biochem.mpg.de/s/VfqtW5p9MJ0kxAC/download?files=20231017_OA2_TiHe_ADIAMA_HeLa_200ng_Evo011_21min_F-40_07.mzML"

EXAMPLE_PRECURSORS_TSV_URL = (
    "https://datashare.biochem.mpg.de/s/VfqtW5p9MJ0kxAC/download?files=precursors.tsv"
)
EXAMPLE_SPECLIB_URL = (
    "https://datashare.biochem.mpg.de/s/VfqtW5p9MJ0kxAC/download?files=speclib.hdf"
)


def load_alphadia_data(
    main_folder: Path,
    *,
    download_urls: tuple[str, str, str] | None = None,
    raw_file_name: str = None,
    precursors_file_name: str = None,
    speclib_file_name: str = None,
) -> tuple[pd.DataFrame, SpecLibBase, SpecLibFlat, Thermo]:
    """Load AlphaDIA results (precursors.tsv + speclib.hdf) and raw data.

    If the download_data is true, this function will download the required data.

    :param main_folder: folder in which data is expected / will be downloaded to
    :param download_urls: a optional tuple (raw_file_url, precursors_tsv_url, speclib_url) with the URLs to download the data
    :param precursors_file_name: name of the precursors file (required if download_data is False)
    :param raw_file_name: name of the raw file (required if download_data is False)
    :param speclib_file_name: name of the speclib file (required if download_data is False)
    :return:
    """
    if download_urls:
        precursors_tsv_path, raw_file_path, speclib_path = (
            _download_alphadia_example_data(main_folder, *download_urls)
        )
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

    print("Reading raw file ...")
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

    print("Parsing spectral library ...")
    spectral_library_flat = SpecLibFlat()
    spectral_library_flat.parse_base_library(spectral_library)

    return precursor_df, spectral_library, spectral_library_flat, dia_data


def _download_alphadia_example_data(
    main_folder: Path, raw_data_url: str, precursors_tsv_url: str, speclib_url: str
) -> tuple[Path, Path, Path]:
    """Download AlphaDIA data if not already present."""
    os.makedirs(main_folder, exist_ok=True)

    raw_file_path = DataShareDownloader(raw_data_url, str(main_folder)).download()

    precursors_tsv_path = DataShareDownloader(
        precursors_tsv_url, str(main_folder)
    ).download()
    speclib_path = DataShareDownloader(speclib_url, str(main_folder)).download()

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
