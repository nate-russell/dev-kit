import pynndescent
import numpy as np
from devkit.utils import tqdm_imap
import numpy as np
import pandas as pd
import umap
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import matplotlib.pyplot as plt
from tqdm.auto import tqdm
from devkit.chemistry import mol_dash_scatter, smi_pair_to_dash



def smi_to_fingerprint(smiles,radius=2, nBits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
        arr = np.zeros((1,), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    else:
        return None


class ActivityCliffAnalysis(object):

    def __init__(self,smiles,y,precomputed_fp=None) -> None:
        pass
        if precomputed_fp is None:
            fingerprints = tqdm_imap(smiles,smi_to_fingerprint,n_procs=1)
        else:
            fingerprints = precomputed_fp

        X = np.vstack(fingerprints)
        index = pynndescent.NNDescent(X,metric='rogerstanimoto')
        index.prepare()

        self.neighbors,self.distances = index.neighbor_graph

        print(self.neighbors.shape)
        print(self.distances.shape)

        #plt.imshow(self.distances,interpolation='nearest', aspect='auto')
        #plt.show()

        pairs = []
        y_deltas = []
        xx_deltas = []
        smiles_pairs = []
        seen_pairs = set()

        for i, (narray,darray) in tqdm(enumerate(zip(self.neighbors,self.distances)),total=len(self.neighbors)):
            # print('mol',i)
            for j in range(len(narray)):
                n_idx = narray[j]
                n_dist = darray[j]
                pair = tuple(sorted([i,n_idx]))
                if pair not in seen_pairs:
                    pairs.append(pair)
                    seen_pairs.add(pair)
                    y_deltas.append(np.abs(y[i] - y[n_idx]))
                    xx_deltas.append(n_dist)
                    smiles_pairs.append(".".join([smiles[i],smiles[n_idx]]))   
            
        mol_dash_scatter(x=y_deltas,
                        y=xx_deltas,
                        c=y_deltas,
                        smi=smiles_pairs,
                        draw_func=smi_pair_to_dash)