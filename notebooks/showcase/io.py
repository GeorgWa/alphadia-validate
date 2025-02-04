import os

import pandas as pd
from alphabase.spectral_library.base import SpecLibBase

from alphadia.data.alpharaw_wrapper import Thermo
from alphadia.test_data_downloader import DataShareDownloader


def prepare_data(main_folder: str) -> tuple[pd.DataFrame, SpecLibBase, Thermo]:
    """Prepare raw & results data.

    Note: this is just required if you did not run the search yourself in search_1.10.0.ipynb
    :param main_folder: folder in which data is expected / will be downloaded to
    :return:
    """
    data_folder = main_folder / "data"
    output_folder = main_folder / "output"
    for folder in [data_folder, output_folder]:
        if not os.path.exists(folder):
            os.makedirs(folder)

    # Bulk injections of HeLa cell lysate acquired on the Orbitrap Astral
    raw_data_url_list = [
        "https://datashare.biochem.mpg.de/s/339jg5HtGrwLwDN/download?files=20231017_OA2_TiHe_ADIAMA_HeLa_200ng_Evo011_21min_F-40_07.raw",
    ]
    # results from search_1.10.0.ipynb
    precursors_tsv_url = "https://datashare.biochem.mpg.de/s/fsCqlT757ttVWI8"
    speclib_url = "https://datashare.biochem.mpg.de/s/VxakpS6mhM2IwxJ"

    raw_file_paths = [DataShareDownloader(url, data_folder).download() for url in raw_data_url_list]
    precursors_tsv_path = DataShareDownloader(precursors_tsv_url, output_folder).download()
    speclib_path = DataShareDownloader(speclib_url, output_folder).download()

    current_raw_name = os.path.basename(raw_file_paths[0]).replace('.raw', '')
    precursor_df = pd.read_csv(precursors_tsv_path, sep='\t')
    precursor_df = precursor_df[precursor_df['run'] == current_raw_name]

    spectral_library = SpecLibBase()
    spectral_library.load_hdf(speclib_path)

    current_raw_path = raw_file_paths[0]
    dia_data = Thermo(current_raw_path)

    return precursor_df, spectral_library, dia_data