import os

from misc.pipeline import run_sdf_mode, topology_builder
from domd_xyz.embed_molecule import embed_molecule
from misc.logger import logger
from rdkit import Chem
from misc.parser import parse_cg_topology
logger.setLevel('INFO')
logger.propagate = True
from misc.io.sdf import write_mols_to_sdf
import numpy as np

reaction_template = {
        'P-P': {
            'cg_reactant_list': [('P', 'P')],
            'smarts': '[C:1][CH1:2].[CH3:3]>>[C:1][C:2][C:3]',
            'prod_idx': [0]
        },
        'CN-A': {
            'cg_reactant_list': [('CN', 'A')],
            'smarts': '[N:1][H:5].[OH1:2][C:3]=[O:4]>>[N:1][C:3]=[O:4].[O:2][H:5]',
            'prod_idx': [0]
        },
        'A-P': {
            'cg_reactant_list': [('A', 'P')],
            'smarts': '[CH2:1][CH1:2].[CH3:3]>>[C:1][C:2][C:3]',
            'prod_idx': [0]
        },
}
mol_config = {
    'P': {'smiles': 'CC(C(=O)OC)C', 'file': None},
    'A': {'smiles': 'OC(=O)CC(C(=O)OC)C', 'file': None},
    'CN': {'smiles': '[N:1]([H:2])[H:3]', 'file': None},

    #'R': {'smiles': None, 'file': 'siRNA.pdb', 'is_rigid':True, 'body_id':0},
}
body_configs = {
    # key is body_id, value is a dict with 'file', 'mapping', and 'reacting_sites'
    0 :{
        'file': 'POSS.pdb',
        'mapping': {
            5:  {'atom_index': [ 0, 32, 33],  'smarts': '[N:1]([H:2])[H:3]'},  # the key is local CG bead index, the value is local atom index
            6:  {'atom_index': [23, 38, 39],  'smarts': '[N:1]([H:2])[H:3]'},
            7:  {'atom_index': [28, 51, 52],  'smarts': '[N:1]([H:2])[H:3]'},
            8:  {'atom_index': [31, 58, 59],  'smarts': '[N:1]([H:2])[H:3]'},
        },
    },
}
# --------------------- CG pre-equilibrium ------------------- #
## please run CG polymerization and pre-equilibrium

#xmlfile = 'POSS-graft-4-PMMA10.xml'
xmlfile = 'final_cg.xml'
cg_mols, box_tensor = parse_cg_topology(xmlfile, mol_config, rigid_configs=body_configs)
cg_mol = cg_mols[0]
aa_mol, aa_graph = topology_builder(mol_config, reaction_template, rigid_configs=body_configs, cg_graph=cg_mol)
aa_mol, aa_graph = embed_molecule(aa_mol, cg_mol, aa_graph, box_tensor, large=200)
aa_mols = run_sdf_mode(mol_config, reaction_template, xmlfile, rigid_configs=body_configs, reactions=None, large=200, chunks_per_d=1, output_sdf_path='POSS-g-PMMA_aa.sdf')
#from rdkit.Chem import AllChem
#AllChem.EmbedMolecule(aa_mol, useRandomCoords=True)
Chem.MolToPDBFile(aa_mol, "POSS-g-PMMA_aa.pdb", flavor=4)


