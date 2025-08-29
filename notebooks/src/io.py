import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from alphabase.peptide.precursor import hash_precursor_df
from alphabase.psm_reader import psm_reader_provider
from alphabase.psm_reader.dia_psm_reader import register_readers
from alphabase.spectral_library.base import SpecLibBase
from alphabase.spectral_library.flat import SpecLibFlat
from alphabase.tools.data_downloader import DataShareDownloader
from alphadia.raw_data.alpharaw_wrapper import MzML, Sciex, Thermo
from alphadia.raw_data.bruker import TimsTOFTranspose
from peptdeep.pretrained_models import ModelManager

# Bulk injections of HeLa cell lysate acquired on the Orbitrap Astral
EXAMPLE_RAW_DATA_URL = "https://datashare.biochem.mpg.de/s/VfqtW5p9MJ0kxAC/download?files=20231017_OA2_TiHe_ADIAMA_HeLa_200ng_Evo011_21min_F-40_07.mzML"

EXAMPLE_PRECURSORS_TSV_URL = (
    "https://datashare.biochem.mpg.de/s/VfqtW5p9MJ0kxAC/download?files=precursors.tsv"
)
EXAMPLE_SPECLIB_URL = (
    "https://datashare.biochem.mpg.de/s/VfqtW5p9MJ0kxAC/download?files=speclib.hdf"
)

EXAMPLE_DIANN_REPORT_URL = (
    "https://datashare.biochem.mpg.de/s/Z8bOcBz6QViCVKN/download?"
)
EXAMPLE_DIANN_RAW_DATA_URL = (
    "https://datashare.biochem.mpg.de/s/ajTYDbBcVH0G6c2/download?"
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
        precursors_tsv_path, raw_file_path, speclib_path = _download_data(
            main_folder, *download_urls
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
    if "run" in precursor_df:
        precursor_df = precursor_df[precursor_df["run"] == current_raw_name]
    elif "raw_name" in precursor_df:
        precursor_df = precursor_df[precursor_df["raw_name"] == current_raw_name]

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


def _download_data(
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


def _download_diann_data(
    main_folder: Path, diann_report_url: str, raw_data_url: str
) -> tuple[Path, Path]:
    """Download DIA-NN data if not already present."""
    os.makedirs(main_folder, exist_ok=True)

    diann_report_path = DataShareDownloader(
        diann_report_url, str(main_folder)
    ).download()
    raw_file_path = DataShareDownloader(raw_data_url, str(main_folder)).download()

    return Path(diann_report_path), Path(raw_file_path)


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


def load_diann_data(
    diann_report_path: str | Path | None = None,
    raw_file_path: str | Path | None = None,
    *,
    main_folder: Path | None = None,
    download_urls: tuple[str, str] | None = None,
    model_device: str = "cpu",
    charged_frag_types: list[str] | None = None,
) -> tuple[pd.DataFrame, SpecLibBase, SpecLibFlat, Any]:
    """Load DIA-NN data and create spectral libraries with automatic instrument/NCE detection.

    Args:
        diann_report_path: Path to DIA-NN report file (e.g., report.parquet)
        raw_file_path: Path to raw file (.mzML, .raw, .wiff, or .d)
        main_folder: Folder for downloading data (required if download_urls provided)
        download_urls: Optional tuple (diann_report_url, raw_file_url) to download data
        model_device: Device for model computation ("cpu" or "gpu")
        charged_frag_types: Fragment ion types to include

    Returns:
        Tuple of (precursor_df, speclib_base, spectral_library_flat, dia_data)
    """
    if download_urls:
        if main_folder is None:
            raise ValueError("main_folder must be provided when using download_urls")
        diann_report_url, raw_file_url = download_urls
        diann_report_path, raw_file_path = _download_diann_data(
            main_folder, diann_report_url, raw_file_url
        )

    if diann_report_path is None or raw_file_path is None:
        raise ValueError(
            "Both diann_report_path and raw_file_path must be provided (either directly or via download_urls)"
        )
    raw_path = Path(raw_file_path)
    if raw_path.suffix.lower() == ".raw":
        dia_data = Thermo(str(raw_path))
        instrument = "Orbitrap"
    else:
        raise ValueError(
            f"Unsupported raw format: {raw_path.suffix}. " "Expecting .raw"
        )

    if "nce" in dia_data.spectrum_df.columns and not dia_data.spectrum_df.empty:
        nce = float(dia_data.spectrum_df["nce"].max())
    else:
        raise ValueError("dia_data.spectrum_df lacks an 'nce' column or is empty")

    register_readers()
    modification_mapping = {
        "Phospho@S": "S(Phospho)",
        "Phospho@T": "T(Phospho)",
        "Phospho@Y": "Y(Phospho)",
    }
    diann_reader = psm_reader_provider.get_reader(
        "diann", modification_mapping=modification_mapping
    )
    diann_df = diann_reader.import_file(str(diann_report_path))

    precursor_df = diann_df.copy()
    precursor_df["instrument"] = instrument
    precursor_df["nce"] = nce
    precursor_df = hash_precursor_df(precursor_df)
    precursor_df = precursor_df.drop_duplicates("mod_seq_charge_hash")
    precursor_df["rt"] = precursor_df["rt"] * 60
    precursor_df["rt_start"] = precursor_df["rt_start"] * 60
    precursor_df["rt_stop"] = precursor_df["rt_stop"] * 60
    assert precursor_df["mod_seq_charge_hash"].is_unique, "Hashes not unique!"

    if charged_frag_types is None:
        charged_frag_types = [
            "b_z1",
            "b_z2",
            "y_z1",
            "y_z2",
            "b_modloss_z1",
            "b_modloss_z2",
            "y_modloss_z1",
            "y_modloss_z2",
        ]

    speclib_base = SpecLibBase(charged_frag_types=charged_frag_types)
    speclib_base.precursor_df = precursor_df
    speclib_base.calc_fragment_mz_df()

    model_mgr = ModelManager(device=model_device)
    model_mgr.instrument = instrument
    model_mgr.nce = nce
    intensity_df = model_mgr.predict_ms2(speclib_base.precursor_df)
    speclib_base._fragment_intensity_df = intensity_df

    speclib_base.precursor_df["nAA"] = speclib_base.precursor_df["sequence"].str.len()
    frag_counts = (speclib_base.precursor_df["nAA"] - 1).to_numpy(dtype=np.int64)
    frag_start_idx = np.concatenate(([0], np.cumsum(frag_counts)[:-1]))
    frag_stop_idx = np.cumsum(frag_counts)
    speclib_base.precursor_df["frag_start_idx"] = frag_start_idx
    speclib_base.precursor_df["frag_stop_idx"] = frag_stop_idx
    assert frag_stop_idx[-1] == len(speclib_base.fragment_mz_df), "Row count mismatch!"

    spectral_library_flat = SpecLibFlat()
    spectral_library_flat.parse_base_library(speclib_base)

    return precursor_df, speclib_base, spectral_library_flat, dia_data
