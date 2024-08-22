from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor
import seaborn as sns
import numpy as np
import rdkit
from collections import defaultdict
from IPython.display import SVG, display
from rdkit import Geometry
import matplotlib.pyplot as plt
import networkx as nx

COMMON_HETEROCYLES_AND_FUNCTIONAL_GROUPS = {
        'Benzene':'c1ccccc1',
        'Pyridine':'n1ccccc1',
        'Nitrile':'C#N',
        'Trifluoromethane':'C(F)(F)F',
        'Ether':'[OD2]([#6])[#6]',
        'Ketone': 'cc(=O)c',

        # CAME FROM GPT - DO NOT TRUST
        'Pyridine': 'n1ccccc1',
        'Furan': 'o1cccc1',
        'Thiophene': 's1cccc1',
        'Pyrrole': 'n1cccc1',
        #'Imidazole': 'n1cnc[nH]1',
        'Piperidine': 'C1CCNCC1',
        'Morpholine': 'C1COCCN1',
        'Oxazole': 'o1cncc1',
        'Thiazole': 's1cncc1',
        'Indole': 'c1ccc2[nH]ccc2c1',
        'Quinoline': 'c1ccc2ncccc2c1',
        'Isoquinoline': 'c1ccc2c(c1)cccc2',
        'Pyrimidine': 'n1ccnc1',
        'Purine': 'n1cnc2c1ncn2',
        'Oxepane': 'C1CCOCC1',
        'Aziridine': 'C1CCN1',
        'Oxirane': 'C1CO1',
        'Thiane': 'C1CS1',
        'Dioxane': 'C1COCCO1',
        'Oxathiane': 'C1COS1',
        'Dithiane': 'C1CSCS1',
        'Dioxolane': 'C1OCOC1',
        'Morpholine oxide': 'C1CO[N+]([O-])C1',
        '1,3,4-Oxadiazole': 'o1nnnco1',
        '1,2,4-Triazole': 'c1nn[nH]n1',
        'Thiirane': 'C1CS1',
        '1,3,2-Dioxaborinane': 'B1COCO1',
        'Oxathiirane': 'C1COCS1',
        '1,2,4-Oxathiazolidine': 'C1CSSO1',
        'Oxadiazole': 'n1nc(no1)',
        'Pyrrolidine': 'C1CCNC1',
        'Diaziridine': 'C1CNN1',
        '1,4-Dioxane': 'C1OCCOCC1',
        '1,4-Dithiane': 'C1SCSCSC1',
        '1,2,3-Triazole': 'c1nn[nH]n1',
        '1,3-Dioxane': 'C1OCCOC1',
        '1,2,4-Oxadiazine': 'o1nnc[nH]o1',
        '1,2,3-Oxadiazole': 'o1nnc[oH]1',
        '1,2,3-Thiadiazole': 's1nnc[nH]s1',
        '1,2,4-Thiadiazole': 's1nnc[sH]1',
        '1,2,3-Triazol-4-one': 'c1nnoc1=O',
        '1,2,3-Triazol-4-ol': 'c1nnoc1',
        '1,2,3-Triazine': 'c1nnnc1',
        '1,2,4-Triazol-5-one': 'c1n[nH]c(=O)n1',
        '1,2,4-Triazol-5-ol': 'c1n[nH]c(=O)o1',
        '1,2,4-Triazol-3-one': 'c1n[nH]c(=O)c1',
        'Alkane': '[#6]',
        'Alkene': '[#6]=[#6]',
        'Alkyne': '[#6]#[#6]',
        'Alkyl halide': '[#6][F,Cl,Br,I]',
        'Alcohol': '[#6][OH]',
        'Ether': '[#6][O][#6]',
        'Aldehyde': '[#6][C]=O',
        'Ketone': '[#6][C](=O)[#6]',
        'Carboxylic acid': '[#6][C](=O)[O][H]',
        'Ester': '[#6][C](=O)[O][#6]',
        'Amide': '[#6][C](=O)[N][#6]',
        'Amine (primary)': '[#6][N][H2]',
        'Amine (secondary)': '[#6][N][H1]([#6])[#6]',
        'Amine (tertiary)': '[#6][N]([#6])[#6]',
        'Nitro group': '[#7](=[O])=[O]',
        'Nitrile': '[#6][C]#N',
        'Sulfide': '[#6][S][#6]',
        'Sulfoxide': '[#6][S](=O)[#6]',
        'Sulfone': '[#6][S](=O)(=O)[#6]',
        'Sulfonic acid': '[#6][S](=O)(=O)[O][H]',
        'Halide': '[F,Cl,Br,I]',
        'Aryl halide': '[#6][F,Cl,Br,I]',
        'Phenol': 'c[OH]',
        'Aromatic amine': 'c[NH2]',
        'Azide': '[N]=[N]=[N]',
        'Diazo group': '[#6]N=N',
        'Epoxide': '[#6]1[#8][#6]1',
        'Thiol': '[#6][SH]',
        'Isocyanate': 'N=C=O',
        'Nitrone': '[#6][N]=[O]',
        # Add more functional groups as needed
    }



def labels_to_colors(x):
    ux = np.unique(x)
    # List of RGB triplets
    rgb_values = sns.color_palette("tab20b", len(ux))
    # Map label to RGB
    color_map = dict(zip(set(ux), rgb_values))
    colors = [color_map[e] for e in x]
    return colors

def hash_to_range(string_list, vals):
    hash_values = [hash(s) for s in string_list]
    range_size = len(vals)
    mapped_values = {s:vals[hash_val % range_size] for s,hash_val in zip(string_list,hash_values)}
    return mapped_values

def to_mol(obj):
    mol = None
    if isinstance(obj,Chem.Mol):
        return obj
    if isinstance(obj,str):
        mol = Chem.MolFromSmiles(obj)
    if mol is not None and isinstance(obj,str):
        mol = Chem.MolFromSmarts(obj)
    return mol

def smarts_to_aidx_bidx(mol,query_mol):

    #print(mol,query_mol)
    query_mol = to_mol(query_mol)
    mol = to_mol(mol)
    #print(mol,query_mol)

    Rhit_ats,Rhit_bonds = [],[]

    if query_mol is not None:
        for hit_ats in mol.GetSubstructMatches(query_mol):
            #print("Query Atoms:",hit_ats)
            hit_bonds = []
            for bond in query_mol.GetBonds():
                aid1 = hit_ats[bond.GetBeginAtomIdx()]
                aid2 = hit_ats[bond.GetEndAtomIdx()]
                hit_bonds.append(mol.GetBondBetweenAtoms(aid1,aid2).GetIdx())

            #print(hit_ats,hit_bonds)
            Rhit_ats.append(hit_ats)
            Rhit_bonds.append(hit_bonds)

    return Rhit_ats,Rhit_bonds

def draw_mol_with_annotations(mol,named_smarts=COMMON_HETEROCYLES_AND_FUNCTIONAL_GROUPS,highlightRadius = 0.5,k=1,h=400,w=400):
    """
    _summary_

    Parameters
    ----------
    mol : _type_
        _description_
    named_smarts : _type_, optional
        _description_, by default COMMON_HETEROCYLES_AND_FUNCTIONAL_GROUPS
    highlightRadius : float, optional
        _description_, by default 0.5
    k : int, optional
        _description_, by default 1
    h : int, optional
        _description_, by default 400
    w : int, optional
        _description_, by default 400

    Returns
    -------
    _type_
        _description_
    """


    rgb_values = sns.color_palette("tab20b", 90)
    #rgb_values = sns.color_palette("hsv", 90)
    name_to_rgb = hash_to_range(named_smarts.keys(),rgb_values)
    #print("Named, Hashed To Colors:",name_to_rgb)

    bid_to_color = defaultdict(list)
    aid_to_color = defaultdict(list)
    
    name_to_aidx = []
    for name,smarts in named_smarts.items():
        try:
            Rhit_ats,Rhit_bonds = smarts_to_aidx_bidx(mol,smarts)
            for index,(hit_ats,hit_bonds) in enumerate(zip(Rhit_ats,Rhit_bonds)):
                aids = []
                for aid in hit_ats:
                    aid_to_color[aid].append((*name_to_rgb[name.split('_')[0]],0.99))
                    aids.append(aid)
                for bid in hit_bonds:
                    bid_to_color[bid].append((*name_to_rgb[name.split('_')[0]],0.99))
                name_to_aidx.append((f"{name}_{index}",aids))
        except IndexError:
            pass
        except TypeError:
            pass


    h = h
    w = w

    x = rdDepictor.Compute2DCoords(mol)
    rdDepictor.StraightenDepiction(mol)

    d2d = Chem.Draw.MolDraw2DSVG(w * 1, h * 1, w, h)
    dopts = d2d.drawOptions()

    dopts.addStereoAnnotation = True
    dopts.maxFontSize = -1
    dopts.baseFontSize = 0.4
    dopts.annotationFontScale = 0.15

    dopts.additionalAtomLabelPadding = 0.15
    dopts.highlightRadius = highlightRadius
    dopts.clearBackground = False
    #dopts.highlightAtoms = hit_ats
    #dopts.highlightBonds = hit_bonds
    #dopts.setHighlightColour((0,.9,.9,.8))
    #dopts.setBackgroundColour((0,.9,.9,.3))
    #d2d.DrawMoleculeWithHighlights(mol,highlightAtoms = dict(aid_to_color),highlightBonds = dict(bid_to_color))

    d2d.DrawMoleculeWithHighlights(mol,'',dict(aid_to_color),dict(bid_to_color),{},{})

    sssr = Chem.GetSSSR(mol)
    atom_coords = mol.GetConformer().GetPositions() # Get Coords
    ring_system_names = set()
    for rs in sssr:
        rs_aidx = set(rs)
        for name,aidx in name_to_aidx:
            name_aidx = set(aidx)
            X,Y = [],[]
            ps = []
            if rs_aidx.issubset(name_aidx):
                for aid in sorted(rs_aidx):
                    for atom_idx, coord in enumerate(atom_coords):
                        if aid == atom_idx:
                            x, y, _ = coord # print(f"Atom {atom_idx + 1}: ({x:.3f}, {y:.3f})")
                            pos = Geometry.Point2D(x,y)
                            ps.append(pos)
                            X.append(x)
                            Y.append(y)
                            ring_system_names.add(name)
            # Fill Poly
            if ps:
                d2d.SetFillPolys(True)
                d2d.SetColour((*name_to_rgb[name.split('_')[0]],0.5))
                d2d.DrawPolygon(ps)

    name_to_centorid = {}
    for name,aidx in name_to_aidx:
        X,Y = [],[]
        ps = []
        for aid in aidx:
            for atom_idx, coord in enumerate(atom_coords):
                if aid == atom_idx:
                    x, y, _ = coord # print(f"Atom {atom_idx + 1}: ({x:.3f}, {y:.3f})")
                    pos = Geometry.Point2D(x,y)
                    ps.append(pos)
                    X.append(x)
                    Y.append(y)

        # Get Centroid
        x_mean = np.mean(X)
        y_mean = np.mean(Y)
        i = Geometry.Point2D(x_mean,y_mean)
        name_to_centorid[name] = i

    # Figure out positions of arrow annotations
    G = nx.Graph()
    fix_pos = {}
    pos = {}
    for atom_idx, coord in enumerate(atom_coords):
        x, y, _ = coord
        G.add_node(atom_idx)
        pos[atom_idx] = (x,y)
    
    for name,i in name_to_centorid.items():
        if name not in ring_system_names:
            label_node = f"{name}_anchor"
            G.add_node(name)
            G.add_node(label_node)
            pos[label_node] = (i.x,i.y)
            G.add_edge(name,label_node)

    fixed_nodes = list(pos.keys())  # Get a list of fixed node names
    pos = nx.spring_layout(G, k=k, pos=pos, fixed=fixed_nodes)
    name_to_arrow_pos = dict()
    for node,pos_arr in pos.items():
        if node in name_to_centorid:
            x,y = pos_arr[0],pos_arr[1]
            name_to_arrow_pos[node] = Geometry.Point2D(x,y)
    
    # Annotate the structures
    for name,aidx in name_to_aidx:
        #print(name,aidx)
        i = name_to_centorid[name]
        if name in ring_system_names:
            d2d.SetColour((*name_to_rgb[name.split('_')[0]],0.5))
            d2d.DrawString(name.split('_')[0],i,0,rawCoords=False)
        else:
            #d2d.SetOffset(0,0)
            #d2d.SetLineWidth(40)
            arrow_label = name_to_arrow_pos[name]
            ni = Geometry.Point2D(float(i.x),float(i.y))
            arrow_end = Geometry.Point2D(float(name_to_arrow_pos[name].x),float(name_to_arrow_pos[name].y))
            offset_percent = 0.15
            x_delta = abs(arrow_end.x - i.x)*offset_percent
            y_delta = abs(arrow_end.y - i.y)*offset_percent

            
            if arrow_end.x > i.x:
                arrow_end.x += -x_delta
                ni.x += x_delta
            else:
                arrow_end.x += x_delta
                ni.x += -x_delta

            if arrow_end.y > i.y:
                arrow_end.y += -y_delta
                ni.y += y_delta
            else:
                arrow_end.y += y_delta
                ni.y += -y_delta

            # if the arrow is mostly vertical then center the label
            # if the arrow is mostly horizontal then left or right justify it
            justify = 0
            deg = np.degrees(np.arctan2(y_delta,x_delta))
            #print(f"Degrees: {deg}")
            if deg < 30:
                if arrow_end.x > i.x:
                    justify = 1
                else:
                    justify = 2

            d2d.SetColour((*name_to_rgb[name.split('_')[0]],1))#(*name_to_rgb[name.split('_')[0]],1)
            d2d.DrawLine(arrow_end,i,rawCoords=False)
            d2d.DrawString(name.split('_')[0],arrow_label,justify,rawCoords=False)

    # Wrap it up
    d2d.FinishDrawing()
    svg_str = d2d.GetDrawingText()
    return svg_str