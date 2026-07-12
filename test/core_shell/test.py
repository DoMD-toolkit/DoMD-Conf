from misc.logger import logger
from misc.pipeline import topology_builder

logger.setLevel('INFO')
logger.propagate = True

smiO = 'OCC'
smiN = 'OC(=O)C(CC(=O)NCCNCCN)N'
smiR = 'CC'
reaction_template = {
    'O-O': {
        'cg_reactant_list': [('O', 'O')],
        'smarts': '[CH3:1].[OH1:2]>>[C:1][O:2]',
        'prod_idx': [0]
    },
    'N-N': {
        'cg_reactant_list': [('N', 'N'), ],
        'smarts': '[C:1](=[O:2])[OH1:3].[NH2:4][CH1:5]>>[C:1](=[O:2])[N:4][C:5].[O:3]',
        'prod_idx': [0]
    },
    'N-O': {
        'cg_reactant_list': [('O', 'N')],
        'smarts': '[CH3:1].[NH2:2][CH1:3]>>[C:1][N:2][C:3]',
        'prod_idx': [0]
    },
}
mol_config = {
    'O': {'smiles': smiO, 'file': None},
    'N': {'smiles': smiN, 'file': None},
    # 'R': {'smiles': None, 'file': 'siRNA.pdb', 'is_rigid':True, 'body_id':0},
}
body_configs = {
    # key is body_id, value is a dict with 'file', 'mapping', and 'reacting_sites'
    0: {
        'file': 'siRNA.pdb',
        'mapping': {
            2: [120],  # the key is local CG bead index, the value is local atom index
            3: [121],
        },
        'reacting_sites': {
            2: {'smiles': 'O', 'type': 'R1'},  # the key is local CG index
            3: {'smiles': 'N', 'type': 'R2'},
        },
    }
}
# --------------------- CG pre-equilibrium ------------------- #
## please run CG polymerization and pre-equilibrium

xmlfile = 'cg.xml'
from misc.parser import parse_cg_topology

cg_mols, box_tensor = parse_cg_topology(xmlfile, mol_config, rigid_configs=body_configs)
print(f"Parsed {len(cg_mols)} coarse-grained molecules from {xmlfile}. Box dimensions: {box_tensor}")
for cg_mol in cg_mols:
    print(
        f"CG molecule has {cg_mol.number_of_nodes()} nodes and {cg_mol.number_of_edges()} edges. Is rigid: {cg_mol.graph['rigidity']}. Type: {cg_mol.graph.get('type', 'N/A')}")
# rdmols = build_aa_topology(mol_config, reaction_template, xmlfile, large=600, chunks_per_d=3)
# output_dir = "output"
# os.makedirs(output_dir, exist_ok=True)
# write_mols_to_sdf(rdmols, os.path.join("core_shell.sdf"),fragments=True)
# from rdkit.Chem import rdMolTransforms
#
## 构建 4x4 的缩放对角矩阵，最后一位固定为 1.0（齐次坐标约束）
# scale_matrix = np.diag([0.1, 0.1, 0.1, 1.0])
#
# for i, mol in enumerate(rdmols):
#    # 1. 结构清洗
#    Chem.SanitizeMol(mol)
#
#    # 2. 提取 Conformer 并通过矩阵操作直接将坐标缩小 10 倍
#    if mol.GetNumConformers() > 0:
#        conf = mol.GetConformer()
#        rdMolTransforms.TransformConformer(conf, scale_matrix)
#
#    # 3. 输出 PDB 文件（此时坐标已是原本的 0.1 倍）
#    output_path = os.path.join(output_dir, f"mol_{i:0>3d}.pdb")
#    Chem.MolToPDBFile(mol, output_path, flavor=4)


# from misc.io.xml import XmlParser
# from misc.parser import read_cg_topology
# import networkx as nx
# xml = XmlParser(xmlfile)
# cg_sys, cg_mols = read_cg_topology(xml, mol_config)
# cg_mol = cg_mols[0]
# i = 0
##for i, cg_mol in enumerate(cg_mols):
# rdmol, mol_graph = topology_builder(mol_config, reaction_template, cg_mol, mol_idx=i)
# rdmol, mol_graph = embed_molecule(rdmol, cg_mol, mol_graph, box=cg_mol.graph['box'], large=600, chunk_per_d=1)
# Chem.MolToPDBFile(rdmol, f"mol_{i:0>3d}.pdb", flavor=4)
