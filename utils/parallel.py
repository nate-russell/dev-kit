from tqdm.auto import tqdm
from multiprocessing import Pool, cpu_count
from time import sleep



def tqdm_imap(iterable,function,n_procs=None,chunksize=100):
    """
    _summary_

    Parameters
    ----------
    iterable : _type_
        _description_
    function : _type_
        _description_
    n_procs : _type_, optional
        _description_, by default None
    chunksize : int, optional
        _description_, by default 100

    Returns
    -------
    _type_
        _description_

    Raises
    ------
    ValueError
        _description_
    """

    if n_procs is None:
        n_procs = cpu_count()
    if not (isinstance(n_procs,int) and n_procs >= 1):
        raise ValueError(f"n_procs needs to be a positive int, it is type {type(n_procs)}")

    if n_procs > 1:
        with Pool(processes=n_procs) as pool:
            p = pool.imap(function, iterable, chunksize=chunksize)
            res = [r for r in tqdm(p,
                                desc=f'IMAP f={function.__name__} N={n_procs} ChunkSize={chunksize}',
                                total=len(iterable))
                                ]
    else:
        # useful for debugging
        res = [function(i) for i in tqdm(iterable,desc=f'Serial f={function.__name__}')]

    return res

def basic(x):
    sleep(0.01)
    return x**2

if __name__ == "__main__":
    import numpy as np
    x = np.random.normal(-1,1,1000)

    tqdm_imap(x,basic)
    tqdm_imap(x,basic,n_procs=1)
    tqdm_imap(x,basic,n_procs='a')