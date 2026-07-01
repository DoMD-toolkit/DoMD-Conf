import numpy as np
from rdkit import Chem
import logging
from domd_xyz.embed_molecule import embed_molecule
from domd_topology.topology_builder import topology_builder
from misc.parser import post_process_aa_mol, parse_cg_xml_topology
logger = logging.getLogger(__name__)


def build_aa_topology(mols_config, reaction_template, xml_path, reactions=None, rigid_meta=None, large=500, chunks_per_d=1):
    """
    Lightweight wrapper to reconstruct AA topology from CG system,
    generate 3D conformers, and inject essential metadata.

    Args:
        mols_config (dict): Monomer configuration dictionary, e.g., {'C': {'smiles': smiC, 'file': None}}
        reaction_template (dict): Reaction SMARTS patterns and topology rules.
        xml_path (str): Path to the GALAMOST CG XML configuration file.
        reactions (list/tuple, optional): Explicit sequence of reactions.
            If None, infers from the XML bond section.
        rigid_meta (dict, optional): Metadata for rigid molecules, e.g., {'C': {'rigid_aidxs_map': {...}}}
        large (int, optional): Threshold for large systems to adjust embedding parameters.
        chunks_per_d (int, optional): Number of chunks per dimension for embedding large systems.

    Returns:
        list[Chem.Mol]: A list of RDKit molecules with conformers and injected metadata.
    """
    logger.info('Parsing CG XML and extracting box dimensions...')

    cg_mols, box_tensor, reactions = parse_cg_xml_topology(xml_path, mols_config, reactions=reactions, rigid_meta=rigid_meta)
    # 5. Serial 3D coordinate embedding and metadata injection
    final_rdmols = []

    for i, cg_mol in enumerate(cg_mols):
        # Topology building
        aa_mol_h, aa_graph = topology_builder(mols_config, reaction_template, cg_mol, mol_idx=i)
        # Molecule embedding (using default safe thresholds for large systems)
        aa_mol_h, aa_graph = embed_molecule(aa_mol_h, cg_mol, aa_graph, box=box_tensor, large=large, chunk_per_d=chunks_per_d)
        # SDF format post processing, inject metadata: resname, res_id, box_tensor
        aa_mol_h = post_process_aa_mol(aa_mol_h, box_tensor)
        final_rdmols.append(aa_mol_h)
        logger.info(f'Successfully processed molecule {i+1} / {len(cg_mols)}')

    return final_rdmols

