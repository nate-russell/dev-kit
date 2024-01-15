import dash
from dash import dcc, html
from dash.dependencies import Input, Output
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from base64 import b64encode
from rdkit import Chem
from rdkit import Geometry
import io
import random, string
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D


# Example data
data = [
    {"smiles": "CCO", "label": 1},
    {"smiles": "CCN", "label": 2},
    {"smiles": "CCOc1ccccc1O", "label": 1},
    {"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O", "label": 2}
]

# Convert RDKit molecule to SVG
def svg_2_dash_img(svg_str, dash_id=None):

    try:
        b64 = b64encode(svg_str.encode('utf-8')).decode('utf-8')
    except AttributeError:
        b64 = b64encode(svg_str).decode('utf-8')

    src_ready_svg = 'data:image/svg+xml;base64, {}'.format(b64)
    if dash_id is not None:
        drawing = html.Img(className='dash_mol_drawing', src=src_ready_svg)
    else:
        drawing = html.Img(id=dash_id, className='dash_mol_drawing', src=src_ready_svg)
    return drawing

def mol_to_svg(mol):
    h = 200
    w = 300
    d2d = Chem.Draw.MolDraw2DSVG(w * 1, h * 1, w, h)
    d2d.drawOptions().addStereoAnnotation = True
    d2d.drawOptions().maxFontSize = -1
    d2d.drawOptions().annotationFontScale = 0.5
    d2d.drawOptions().additionalAtomLabelPadding = 0.15
    d2d.SetFontSize(6)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    svg_str = d2d.GetDrawingText()
    return svg_str


def smi_to_dash(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol:
        svg_data = mol_to_svg(mol)
        return svg_2_dash_img(svg_data,dash_id='hover_mol_img')
    else:
        return html.Div(f'{smi} is not valid molecule')
    
def smi_pair_to_dash(smi):

    two_smi = smi.split('.')
    assert len(two_smi) == 2
    smiles1,smiles2 = two_smi
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    ps = rdFMCS.MCSParameters()
    ps.AtomCompareParameters.RingMatchesRingOnly = True
    ps.AtomTyper = rdFMCS.AtomCompare.CompareAny
    mcs = rdFMCS.FindMCS([mol1,mol2],ps)
    core=Chem.MolFromSmarts(mcs.smartsString)

    try:
        AllChem.Compute2DCoords(mol1)
        match1=mol1.GetSubstructMatch(core)
        match2=mol2.GetSubstructMatch(core)
        coords=[mol1.GetConformer().GetAtomPosition(x) for x in match1]
        coords2D = [Geometry.Point2D(pt.x,pt.y) for pt in coords]
        coordDict={}
        for i,coord in enumerate(coords2D):
            coordDict[match2[i]] = coord
        AllChem.Compute2DCoords(mol2,coordMap=coordDict)
    except Exception:
        pass

    patt = core

    hit_ats = list(mol1.GetSubstructMatch(patt))
    hit_bonds = []
    for bond in patt.GetBonds():
        aid1 = hit_ats[bond.GetBeginAtomIdx()]
        aid2 = hit_ats[bond.GetEndAtomIdx()]
        hit_bonds.append(mol1.GetBondBetweenAtoms(aid1,aid2).GetIdx())


    h = 400
    w = 1500
    d2d = Chem.Draw.MolDraw2DSVG(w * 1, h * 1, w, h)
    dopts = d2d.drawOptions()
    dopts.addStereoAnnotation = True
    dopts.maxFontSize = -1
    dopts.annotationFontScale = 0.5
    dopts.additionalAtomLabelPadding = 0.15
    dopts.highlightAtoms = hit_ats
    dopts.highlightBonds = hit_bonds
    dopts.setHighlightColour((0,.9,.9,.8))
    #dopts.setBackgroundColour((0,.9,.9,.3))
    d2d.DrawMolecule(mol1,highlightAtoms = hit_ats,highlightBonds = hit_bonds)
    d2d.FinishDrawing()
    svg_str1 = d2d.GetDrawingText()

    #    d1 =  Chem.Draw.MolDraw2DSVG(1500, 500) # or MolDraw2DCairo to get PNGs
    #rdMolDraw2D.PrepareAndDrawMolecule(d1, mol1, highlightAtoms=hit_ats,highlightBonds=hit_bonds)
    #svg_str1 = d1.GetDrawingText()

    hit_ats = list(mol2.GetSubstructMatch(patt))
    hit_bonds = []
    for bond in patt.GetBonds():
        aid1 = hit_ats[bond.GetBeginAtomIdx()]
        aid2 = hit_ats[bond.GetEndAtomIdx()]
        hit_bonds.append(mol2.GetBondBetweenAtoms(aid1,aid2).GetIdx())

    h = 400
    w = 1500
    d2d = Chem.Draw.MolDraw2DSVG(w * 1, h * 1, w, h)
    dopts = d2d.drawOptions()
    dopts.addStereoAnnotation = True
    dopts.maxFontSize = -1
    dopts.annotationFontScale = 0.5
    dopts.additionalAtomLabelPadding = 0.15
    dopts.highlightAtoms = hit_ats
    dopts.highlightBonds = hit_bonds
    dopts.setHighlightColour((0,.9,.9,.8))
    #dopts.setBackgroundColour((0,.9,.9,.3))
    d2d.DrawMolecule(mol2,highlightAtoms = hit_ats,highlightBonds = hit_bonds)
    d2d.FinishDrawing()
    svg_str2 = d2d.GetDrawingText()

    img1 = svg_2_dash_img(svg_str1, dash_id="pair mol 1")
    img2 = svg_2_dash_img(svg_str2, dash_id="pair mol 2")


    return html.Div([img1,img2])


import dash_bootstrap_components as dbc

def mol_dash_scatter(x,y,c,smi,draw_func=smi_to_dash):
    # Create Dash app layout
    rstr = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))
    app = dash.Dash(f"{__name__}_{rstr}",external_stylesheets=[dbc.themes.ZEPHYR])

    import plotly.graph_objects as go # or plotly.express as px
    fig = go.Figure()

    fig.add_trace(
        go.Scattergl(
            x = x,
            y = y,
            mode = 'markers',
            text=smi,
            marker={'size': 12,'color': c,'colorscale': 'Viridis'}
        )
    )
    fig.update_layout(
        hovermode = 'closest',
        title = 'Scatter Plot of Molecules',
        height = 800,
        width = 800
    )
    graph = dcc.Graph(figure=fig,id='scatter-plot')

    """
    graph = dcc.Graph(
            id='scatter-plot',
            figure={
                'data': [
                    {
                        'x': x,
                        'y': y,
                        'text': smi,
                        'mode': 'markers',
                        'marker': {
                            'size': 12,
                            'color': c,
                            'colorscale': 'Viridis'
                        }
                    }
                ],
                'layout': {
                    'hovermode': 'closest',
                    'title': 'Scatter Plot of Molecules',
                    'height':'800px',
                    'width':'800px'
                }
            }
        )
    """
    row = dbc.Row(
            [
                dbc.Col(html.Div([graph])),
                dbc.Col(html.Div(id='molecule-svg-div')),
            ]
        )
    
    app.layout = html.Div([row])

    # Callback to update the SVG on hover
    @app.callback(
        Output('molecule-svg-div', 'children'),
        [Input('scatter-plot', 'hoverData')]
    )
    def update_svg_on_hover(hover_data):
        if hover_data:
            selected_smiles = hover_data['points'][0]['text']
    
            dash_img = draw_func(selected_smiles)
            return html.Div([
                html.H5(f'Molecule: {selected_smiles}'),
                dash_img
            ])
            
    app.run_server(debug=True)
