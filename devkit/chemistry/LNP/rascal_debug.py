"""
rdRascalMCES.FindMCES is great but the maxBondMatchPairs = 1000; limit is making it hard to work with lipids. 
wrapper for options: https://github.com/rdkit/rdkit/blob/73b4da2a44a297758def89dd4f9eda8cf0b4e907/Code/GraphMol/RascalMCES/Wrap/rdRascalMCES.cpp#L155
the actual options: https://github.com/rdkit/rdkit/blob/73b4da2a44a297758def89dd4f9eda8cf0b4e907/Code/GraphMol/RascalMCES/RascalOptions.h
the code https://github.com/rdkit/rdkit/blob/73b4da2a44a297758def89dd4f9eda8cf0b4e907/Code/GraphMol/RascalMCES/RascalMCES.cpp
"""

import rdkit
print('rdkit version: ',rdkit.__version__)
import sys
print(sys.version)

from rdkit import Chem
from rdkit.Chem import rdRascalMCES

def rascal_wrapper(mola,molb):
    opts = rdRascalMCES.RascalOptions()
    opts.similarityThreshold = 0.0
    opts.returnEmptyMCES = True
    opts.singleLargestFrag = True
    opts.allBestMCESs = True
    opts.maxBondMatchPairs = 10000
    opts.completeAromaticRings = False
    opts.timeout = -1
    results = rdRascalMCES.FindMCES(mola, molb,opts)
    print(f'Results {results}')
    print(f"N Results: {len(results)}")

    print("tier1sim:",results[0].tier1Sim,"tier2sim:",results[0].tier2Sim)
    
    for i,res in enumerate(results):
        print(f"Res {i}")
        print(f'MCES SMARTS : {res.smartsString}')
        print(f'Matching Bonds : {res.bondMatches()}')
        print(f'Matching Atoms : {res.atomMatches()}')


# These are the Lipids i actually wanted to compare
print('\nLipids I care about')
real_mol_1 =  Chem.MolFromSmiles('CCCCCCCC\C=C/CCCCCCCCNC(=O)C(CCCCCOC(=O)CCCCCCC\C=C/CCCCCCCC)NCCCN(C)CCCN')
real_mol_2 =  Chem.MolFromSmiles('CCCCCCCCCCC(=O)OCCCCCC(NCCCN(C)CCCN)C(=O)NCCCCCCCC\C=C/CCCCCCCC')
rascal_wrapper(real_mol_1,real_mol_2)


# Larger heteratom ring systems seem to work ok so i think my issue is related to carbon chains
# This was an example molecukle pair i made to test the limits
print('\nReturns results')
gets_results_1 =  Chem.MolFromSmiles('CCCC=CCCCC=CCCCCCCCCCCCCCCCCCCCC1CNCCC1')
gets_results_2 =  Chem.MolFromSmiles('CCCC=CCCCCC=CCCCCCCCCCCCCCCCCCCC1CNCCC1')
rascal_wrapper(gets_results_1,gets_results_2)

# I add one extra carbon and I hit the limit
print('\nNo Results: 1 carbon too long resulting in Too many matching bond pairs (1032) so can\'t continue.')
too_long_1 =  Chem.MolFromSmiles('CCCC=CCCCC=CCCCCCCCCCCCCCCCCCCCCC1CNCCC1')
too_long_2 =  Chem.MolFromSmiles('CCCC=CCCCCC=CCCCCCCCCCCCCCCCCCCCC1CNCCC1')
rascal_wrapper(too_long_1,too_long_2)