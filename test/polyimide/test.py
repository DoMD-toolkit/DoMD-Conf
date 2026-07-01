import os

from misc.pipeline import build_aa_topology, topology_builder
from misc.logger import logger
from rdkit import Chem
logger.setLevel('INFO')
logger.propagate = True
from misc.io.sdf import write_mols_to_sdf

smiA = 'Nc1ccc(Oc2ccc(N)cc2)cc1'
smiB = 'O1C(=O)c2cc(Oc3cc4C(=O)OC(=O)c4cc3)ccc2C1=O'

reaction_template = {
        'B-A': {
            'cg_reactant_list': [('A', 'B')],
            'smarts': '[#7H2:1].[#6:3](=[#8:4])[#8:2][#6:5]=[#8:6]>>[#6:3](=[#8:4])[#7:1][#6:5]=[#8:6].[#8:2]',
            'prod_idx': [0]
        }
    }
mol_config = {
    'A': {'smiles': smiA, 'file': None},
    'B': {'smiles': smiB, 'file': None},
}

config  = {
'reactions_template' : {
        'B-A': {
            'cg_reactant_list': [('A', 'B')],
            'smarts': '[#7H2:1].[#6:3](=[#8:4])[#8:2][#6:5]=[#8:6]>>[#6:3](=[#8:4])[#7:1][#6:5]=[#8:6].[#8:2]',
            'prod_idx': [0]
        }
    },
'reactants' : {
    'A': {'smiles': 'Nc1ccc(Oc2ccc(N)cc2)cc1', 'file': None},
    'B': {'smiles': 'O1C(=O)c2cc(Oc3cc4C(=O)OC(=O)c4cc3)ccc2C1=O', 'file': None},
}
}
xmlfile = 'cg.xml'

rdmols = build_aa_topology(mol_config, reaction_template, xmlfile, large=200, chunks_per_d=1)
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)
for i, mol in enumerate(rdmols):
    fname = f"mol_{i:0>3d}.pdb"
    fpath = os.path.join(output_dir, fname)
    Chem.MolToPDBFile(mol, fpath, flavor=4)
#from misc.io.xml import XmlParser
#from misc.parser import read_cg_topology
#xml = XmlParser(xmlfile)
#cg_sys, cg_mols = read_cg_topology(xml, mol_config)
#cg_mol = cg_mols[0]
#i = 0
##for i, cg_mol in enumerate(cg_mols):
#rdmol, mol_graph = topology_builder(mol_config, reaction_template, cg_mol, mol_idx=i)
#for g in mol_graph.nodes(data=True):
#    print(g)
