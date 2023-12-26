import os
from tqdm.auto import tqdm
import requests
import warnings
import pandas as pd
_sqlit_chembl = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_33_sqlite.tar.gz"
ROOT_DATA_DIR = "F:\data"


def fetch_file(url,out_path,verbose=True):
    """
    downloads url contents to out_path
    https://stackoverflow.com/questions/37573483/progress-bar-while-download-file-over-http-with-requests

    Parameters
    ----------
    url : _type_
        _description_
    out_path : _type_
        _description_
    verbose : bool, optional
        _description_, by default True

    Returns
    -------
    _type_
        _description_
    """
    # Streaming, so we can iterate over the response.
    response = requests.get(url, stream=True)
    total_size_in_bytes= int(response.headers.get('content-length', 0))
    block_size = 1024 #1 Kibibyte
    if verbose:
        print(f'{url} -> {repr(out_path)}')

    progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
    with open(out_path, 'wb') as file:
        for data in response.iter_content(block_size):
            progress_bar.update(len(data))
            file.write(data)
    progress_bar.close()
    if total_size_in_bytes != 0 and progress_bar.n != total_size_in_bytes:
        print("ERROR, something went wrong")
    return out_path

def make_dir_path(path):
    mode = 0o666
    if not os.path.exists(path):
        os.makedirs(path, mode)

def download_iuphar():
    d = os.path.abspath(ROOT_DATA_DIR)
    make_dir_path(os.path.join(d, 'IUPHAR'))

    file_path = os.path.join(d,'IUPHAR','interactions.csv')
    fetch_file('https://www.guidetopharmacology.org/DATA/interactions.csv',file_path)


def download_drugs_and_probes():
    targets = 'https://www.probes-drugs.org/media/download/targets/pd_export_02_2023_targets_standardized.csv'
    compounds = 'https://www.probes-drugs.org/media/download/compounds/pd_export_02_2023_compounds_standardized.csv'
    
    d = os.path.abspath(ROOT_DATA_DIR)
    make_dir_path(os.path.join(d, 'PROBES_AND_DRUGS'))

    t_path = os.path.join(d,'PROBES_AND_DRUGS','targets.csv')
    c_path = os.path.join(d,'PROBES_AND_DRUGS','compounds.csv')
    fetch_file(targets,t_path)
    fetch_file(compounds,c_path)




def fetch_data():
    # Check if data exists and it's hash is right
    # if data isn't read fetch it
    # return data
    pass

class DataCatalog:

    def __init__(self) -> None:
        pass

    def catalog_as_df(self) -> pd.DataFrame:
        pass

    def catalog_as_records_dict(self) -> dict:
        pass

    def load_df(self,file) -> pd.DataFrame:
        pass

    




aliagas_benchmark = pd.read_csv(
    "https://static-content.springer.com/esm/art%3A10.1007%2Fs10822-015-9838-3/MediaObjects/10822_2015_9838_MOESM1_ESM.csv"
)

def download_biodgen_adme():
    d = os.path.abspath(ROOT_DATA_DIR)
    make_dir_path(os.path.join(d, 'Biogen_ADME'))

    file_path = os.path.join(d,'Biogen_ADME','ADME_public_set_3521.csv')
    url = "https://raw.githubusercontent.com/molecularinformatics/Computational-ADME/main/ADME_public_set_3521.csv"
    fetch_file(url,file_path)
    
def test_data_dir_write():
    """
    _summary_

    Raises
    ------
    ValueError
        _description_
    FileExistsError
        _description_
    """
    d = os.path.abspath(ROOT_DATA_DIR)
    txt = 'test file used to test wrx permissions'

    if os.path.exists(d):
        test_file_path = os.path.join(d,'test.txt')
        with open(test_file_path,mode='w') as f:
            f.write(txt)
        with open(test_file_path,mode='r') as f:
            rtxt = f.read()
        if rtxt != txt:
            raise ValueError(f"{test_file_path} contents are unexpected. {repr(rtxt)} != {repr(txt)}")
    else:
        raise FileExistsError(f"Couldn't find: {repr(d)}")
    

if __name__ == "__main__":
    test_data_dir_write()
    download_drugs_and_probes()
    #download_iuphar()
    print('all systems go')

