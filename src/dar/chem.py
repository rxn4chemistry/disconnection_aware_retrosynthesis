from __future__ import annotations

from collections import OrderedDict
from enum import Enum, auto, unique
from typing import List, Set

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rxn.chemutils.conversion import canonicalize_smiles


@unique
class AtomEnvironment(Enum):
    CHANGED = auto()
    SAME = auto()


def remove_mapping(smiles: str) -> str:
    """
    Removes mapping from a SMILES string.

    Args:
        smiles: Atom-Mapped SMILES string

    Returns:
        SMILES string without atom-mapping
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)

    [atom.ClearProp("molAtomMapNumber") for atom in mol.GetAtoms()]

    return Chem.MolToSmiles(mol)


def remove_rxn_mapping(rxn_smiles: str) -> str:
    """
    Removes mapping from a reaction SMILES.

    Args:
        rxn_smiles: Atom-Mapped reaction SMILES

    Returns:
        Reaction SMILES without atom-mapping
    """
    rd_rxn = rdChemReactions.ReactionFromSmarts(rxn_smiles)
    rdChemReactions.RemoveMappingNumbersFromReactions(rd_rxn)

    return rdChemReactions.ReactionToSmiles(rd_rxn)


def standardise_reaction_component(component_smiles: str) -> str:
    """
    Standardises reaction component given a SMILES string
    Reaction component example:
    'BrBr.CC(C)(c1ccc(C(=O)C2CC2)cc1)C1CC1.O'

    Args:
        component_smiles: SMILES string
    Returns:
        Standardised SMILES or None if errors
    """
    try:
        return ".".join(sorted(canonicalize_smiles(component_smiles).split(".")))
    except Exception:
        return ""


def remove_unmapped_components(component: str) -> str:
    """
    Remove non-mapped components from a reactant/reagent/product smiles string

    Args:
        component: component SMILES string potentially containing mapping information.
                    e.g. CC(C)(C)[SiH2][O:20][C:19](C)(C)[c:18]1[n:14]([CH:13]2[C:2]....

    Returns:
        component SMILES string with atom mapping information and non-mapped species
        removed or empty string if error
    """

    try:
        component_list = component.split(".")
        component_mols = [Chem.MolFromSmiles(smiles) for smiles in component_list]
        mapped_mols = []

        for mol in component_mols:
            mapping_list = [atom.GetAtomMapNum() for atom in mol.GetAtoms()]
            if mapping_list.count(0) != len(mapping_list):
                mapped_mols.append(mol)

        mapped_smiles = ".".join([Chem.MolToSmiles(mol) for mol in mapped_mols])

        return mapped_smiles

    except Exception:
        return ""


def get_atomic_neighbourhoods(smiles: str) -> OrderedDict[int, List[str]]:
    """
    Obtains a dictionary containing each atomIdx and a list of its bonding environment.

    Args:
        smiles: Atom-Mapped SMILES string

    Returns:
        A dictionary containing each atomIdx and a list of its bonding environment.
    """

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    atoms = mol.GetAtoms()

    neighbour_dict = {}
    for i, atom in enumerate(atoms):
        precursor_bonds_list = []
        if atom.GetAtomMapNum() != 0:
            for bond in atom.GetBonds():
                atom_map1 = atoms[bond.GetBeginAtomIdx()].GetAtomMapNum()
                atom_map2 = atoms[bond.GetEndAtomIdx()].GetAtomMapNum()
                bond_order = bond.GetBondType()
                sorted_atom_maps = sorted([atom_map1, atom_map2])
                precursor_bonds_list.append(
                    f"{sorted_atom_maps[0]}_{sorted_atom_maps[1]}_{bond_order}"
                )
                neighbour_dict[atom.GetAtomMapNum()] = sorted(precursor_bonds_list)
    ordered_neighbour_dict = OrderedDict(sorted(neighbour_dict.items()))

    return ordered_neighbour_dict


def get_all_atom_indices(precursor_smiles: str, product_smiles: str) -> Set[int]:
    """
    Retrieves all atomIdxs common between precursors and products for a reaction
    given the SMILES strings of precursors and products.

    Args:
        precursor_smiles: Atom-mapped SMILES string for the precursor(s)
        product_smiles: Atom-mapped SMILES string for the product(s)
    Returns:
        Set of all AtomIdxs common between precursors and products
    """

    ordered_precursor_neighbour_dict = get_atomic_neighbourhoods(precursor_smiles)
    ordered_product_neighbour_dict = get_atomic_neighbourhoods(product_smiles)

    all_indices = set(ordered_product_neighbour_dict.keys()) | set(
        ordered_precursor_neighbour_dict.keys()
    )

    return all_indices


def get_atom_list(
    precursor_smiles: str, product_smiles: str, atom_environment: AtomEnvironment
) -> List[int]:
    """
    Given two sets of SMILES strings corresponding to a set of precursors and products,
    obtains a list of atomIdxs for which the atomic environment has changed,
    as defined by a change in the bonds.

    Args:
        precursor_smiles: Atom-mapped SMILES string for the precursor(s)
        product_smiles: Atom-mapped SMILES string for the product(s)
        atom_environment: "changed" for changed atoms
                        or "same" for list of equivalent atoms

    Returns:
        List of atomIdxs for which the atomic environment has changed
    """

    ordered_precursor_neighbour_dict = get_atomic_neighbourhoods(precursor_smiles)
    ordered_product_neighbour_dict = get_atomic_neighbourhoods(product_smiles)

    all_indices = set(ordered_product_neighbour_dict.keys()) | set(
        ordered_precursor_neighbour_dict.keys()
    )

    if atom_environment == AtomEnvironment.CHANGED:
        # Checks to see equivlence of atomic enviroments (AE).
        # If not equivelent, then changed AE and added to list
        atom_list = [
            atom_map
            for atom_map in all_indices
            if ordered_precursor_neighbour_dict.get(atom_map, [])
            != ordered_product_neighbour_dict.get(atom_map, [])
        ]
    elif atom_environment == AtomEnvironment.SAME:
        # Checks to see equivlence of atomic enviroments (AE).
        # If equivelent, then unchanged AE and added to list
        atom_list = [
            atom_map
            for atom_map in all_indices
            if ordered_precursor_neighbour_dict.get(atom_map, [])
            == ordered_product_neighbour_dict.get(atom_map, [])
        ]
    else:
        raise TypeError(
            """Unrecognised type: Use 'AtomEnvironment.CHANGED' for a list of changed
            atomIdxs or 'AtomEnvironment.SAME'
            for a list of unchanged atomIdxs"""
        )

    return atom_list
