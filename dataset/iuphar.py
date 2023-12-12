import os
from tqdm.auto import tqdm
import requests
_sqlit_chembl = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_33_sqlite.tar.gz"
ROOT_DATA_DIR = "F:\data"


def fetch_file(url,out_path,verbose=True):
    """
    _summary_

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


def download_iuphar():
    d = os.path.abspath(ROOT_DATA_DIR)
    mode = 0o666
    path = os.path.join(d, 'IUPHAR')
    if not os.path.exists(path):
        os.makedirs(path, mode)

    file_path = os.path.join(d,'IUPHAR','interactions.csv')
    fetch_file('https://www.guidetopharmacology.org/DATA/interactions.csv',file_path)


def test_data_dir_write():
    
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
    download_iuphar()
    print('all systems go')

