from rdkit import Chem

from embed_molecule import embed_molecule
from misc.logger import logger
from misc.parser import parse_cg_topology
from misc.pipeline import topology_builder

logger.setLevel('INFO')
logger.propagate = True

reaction_template = {
    'a': {
        'cg_reactant_list': [('PEO', 'PEO')],
        'smarts': '[C:1][O:2].[O:3][C:4]>>[C:1][O:2][C:4].[O:3]',
        'prod_idx': [0]
    },
    'b': {
        'cg_reactant_list': [('Au_G', 'PEO')],
        'smarts': '[Au:1].[O:2][C:3]>>[Au:1][O:2][C:3]',
        'prod_idx': [0]
    },
}
mol_config = {
    'PEO': {'smiles': 'OCCO', 'file': None},
    'Au_G': {'smiles': '[Au]', 'file': None},
    # 'R': {'smiles': None, 'file': 'siRNA.pdb', 'is_rigid':True, 'body_id':0},
}
body_configs = {
    # key is body_id, value is a dict with 'file', 'mapping', and 'reacting_sites'
    0: {
        'file': 'au.pdb',
        'mapping': {
            0: {'atom_index': [0], 'smarts': '[Au]'},  # the key is local CG bead index, the value is local atom index
            18: {'atom_index': [18], 'smarts': '[Au]'},
            49: {'atom_index': [49], 'smarts': '[Au]'},
            50: {'atom_index': [50], 'smarts': '[Au]'},
            146: {'atom_index': [146], 'smarts': '[Au]'},
        },
    },
    1: {
        'file': 'au.pdb',
        'mapping': {
            49: {'atom_index': [49], 'smarts': '[Au]'},
        },
    }
}
# --------------------- CG pre-equilibrium ------------------- #
## please run CG polymerization and pre-equilibrium

xmlfile = 'out_au_peo_au_large.xml'
cg_mols, box_tensor = parse_cg_topology(xmlfile, mol_config, rigid_configs=body_configs)
cg_mol = cg_mols[0]
aa_mol, aa_graph = topology_builder(mol_config, reaction_template, rigid_configs=body_configs, cg_graph=cg_mol)
aa_mol, aa_graph = embed_molecule(aa_mol, cg_mol, aa_graph, box_tensor, large=200)

# from rdkit.Chem import AllChem
# AllChem.EmbedMolecule(aa_mol, useRandomCoords=True)
Chem.MolToPDBFile(aa_mol, "au_aa_embed.pdb", flavor=4)
