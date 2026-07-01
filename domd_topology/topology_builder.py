from .reactor import Reactor

def topology_builder(reactants_config, reaction_template, cg_graph=None, reactions=None, mol_idx=0):
    """Builds the all-atom topology from coarse-grained input using the Reactor class.

        Args:
            reactants_config (dict): Configuration for reactant molecules, including SMILES and file paths.
            reaction_template (dict): Reaction SMARTS patterns and topology rules.
            cg_graph (networkx.Graph, optional): Coarse-grained graph representation of the system.
            reactions (list/tuple, optional): Explicit sequence of reactions. If None, inferred from cg_graph.
            mol_idx (int, optional): Index of the molecule being processed (for logging purposes).

        Returns:
            tuple: A tuple containing:
                - list[Chem.Mol]: List of reconstructed all-atom molecules.
                - dict: Metadata associated with the reconstructed molecules.
    """
    if cg_graph is None and reactions is None:
        raise ValueError("Either cg_graph or reactions must be provided for topology building.")
    if reactions is None:
        reactions = []
        for i,j, bond in cg_graph.edges(data=True):
            bondtype = bond['bond_type']
            reactions.append((bondtype,i,j))
    reactor = Reactor(reactants_config, reaction_template)
    aa_mol_h, aa_graph = reactor.process(cg_graph, reactions, mol_idx=mol_idx)

    return aa_mol_h, aa_graph