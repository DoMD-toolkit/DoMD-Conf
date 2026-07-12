import numpy as np
from rdkit import Chem

from embed_molecule import embed_molecule
from misc.logger import logger
from misc.pipeline import topology_builder

logger.setLevel('INFO')
logger.propagate = True

smiC = 'C=CC(=C)C'

reaction_template = {
    'C-C': {
        'cg_reactant_list': [('C', 'C')],
        'smarts': '[C:1]=[C:2][C:3](C)=[C:4].[C:5]=[CH1:6]>>[C:1]/[C:2]=[C:3](C)/[C:4][C:5]=[C:6]',
        'prod_idx': [0]
    },
}
mol_config = {
    'C': {'smiles': smiC, 'file': None},
}

config = {
    'reactions_template': {
        'C-C': {
            'cg_reactant_list': [('C', 'C')],
            'smarts': '[C:1]=[C:2][C:3](C)=[C:4].[C:5]=[CH1:6]>>[C:1]/[C:2]=[C:3](C)/[C:4][C:5]=[C:6]',
            'prod_idx': [0]
        },
    },
    'reactants': {
        'C': {'smiles': 'C=CC(=C)C', 'file': None},
    }
}

xmlfile = 'cg.xml'

from misc.io.xml import XmlParser
from misc.parser import read_cg_topology

xml = XmlParser(xmlfile)
cg_sys, cg_mols = read_cg_topology(xml, mol_config)
box = np.array([xml.box.lx, xml.box.ly, xml.box.lz]).astype(float) * 10
cg_mol = cg_mols[0]
i = 0
# for i, cg_mol in enumerate(cg_mols):
rdmol, mol_graph = topology_builder(mol_config, reaction_template, cg_mol, mol_idx=i)
print(
    f"Embedding molecule {i} with {rdmol.GetNumAtoms()} atoms (molecule graph with {mol_graph.number_of_nodes()} nodes) and {cg_mol.number_of_nodes()} CG beads...")
import networkx as nx

mol_graph = nx.Graph()
for atom in rdmol.GetAtoms():
    atom_id = atom.GetIdx()
    mol_graph.add_node(atom_id, res_id=atom.GetIntProp('res_id'), res_name=atom.GetProp('res_name'),
                       global_res_id=atom.GetIntProp('global_res_id'))
# raise
conf = embed_molecule(rdmol, cg_mol, mol_graph, box=box)
if rdmol.GetNumConformers() == 0 and conf is not None:
    rdmol.AddConformer(conf, assignId=True)
Chem.MolToPDBFile(rdmol, f"mol_{i}.pdb", flavor=4)
