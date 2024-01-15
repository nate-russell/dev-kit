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
from scipy.spatial.distance import pdist, squareform
import numpy as np
import seaborn as sns


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
            self.fingerprints = tqdm_imap(smiles,smi_to_fingerprint,n_procs=1)
        else:
            self.fingerprints = precomputed_fp

        self.X = np.vstack(self.fingerprints)
        self.smiles = smiles
        self.y = y
        self.pairs = []
        self.y_deltas = []
        self.xx_deltas = []
        self.smiles_pairs = []
        self.seen_pairs = set()

    def add_pynn_neighbors(self,n=30):

        index = pynndescent.NNDescent(self.X,n_neighbors=n,metric='rogerstanimoto')
        index.prepare()

        self.neighbors,self.distances = index.neighbor_graph

        for i, (narray,darray) in tqdm(enumerate(zip(self.neighbors,self.distances)),total=len(self.neighbors)):
            for j in range(len(narray)):
                n_idx = narray[j]
                n_dist = darray[j]
                pair = tuple(sorted([i,n_idx]))
                if pair not in self.seen_pairs:
                    self.pairs.append(pair)
                    self.seen_pairs.add(pair)
                    self.y_deltas.append(np.abs(self.y[i] - self.y[n_idx]))
                    self.xx_deltas.append(n_dist)
                    self.smiles_pairs.append(".".join([self.smiles[i],self.smiles[n_idx]])) 

    def nn_imshow(self):
        plt.imshow(self.distances,interpolation='nearest', aspect='auto')
        plt.show()
                    
    def add_dense_pairs(self):
        pairwise_distances = squareform(pdist(self.fingerprints, metric='rogerstanimoto'))
        N,_ = pairwise_distances.shape
        for i in range(N):
            for j in range(i,N):
                pair = tuple(sorted([i,j]))
                if pair not in self.seen_pairs:
                    self.pairs.append(pair)
                    self.seen_pairs.add(pair)
                    self.y_deltas.append(np.abs(self.y[i] - self.y[j]))
                    self.xx_deltas.append(pairwise_distances[i,j])
                    self.smiles_pairs.append(".".join([self.smiles[i],self.smiles[j]]))
        
        plt.imshow(pairwise_distances)
        plt.colorbar()
        plt.xlabel("X")
        plt.ylabel("X'")
        plt.show()

        sns.histplot(pairwise_distances.flatten(),bins=50)
        plt.title("Pairwise Distance Histogram\nRogersTanimoto(X,X')")
        plt.show()

        sns.displot(x=self.y_deltas, y=self.xx_deltas,bins=(20,20),
                    #binwidth=(1, 0.01),
                      cbar=True)
        plt.xlabel("Abs(y-y')")
        plt.ylabel("RogersTanimoto(X,X')")
        plt.show()



    def dash_scatter(self):
        mol_dash_scatter(x=self.y_deltas,
                        y=self.xx_deltas,
                        c=self.y_deltas,
                        smi=self.smiles_pairs,
                        draw_func=smi_pair_to_dash)