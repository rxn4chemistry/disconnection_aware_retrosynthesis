import random
import re
from itertools import chain, combinations
from typing import Dict, List, Tuple

from rdkit import Chem

from disconnection_aware_retrosynthesis.chem import AtomEnvironment, get_atom_list


def get_tagged_products(precursor_smiles: str, product_smiles: str) -> str:
    """
    Given two sets of SMILES strings corresponding to a set of precursors and products,
    tags the changed atoms in the product molecule.

    Args:
        precursor_smiles: Atom-mapped SMILES string for the precursor(s)
        product_smiles: Atom-mapped SMILES string for the product(s)

    Returns:
        SMILES of the product containing tags corresponding to atoms changed in the
        reaction using [<atom>:1]
    """

    products_mol = Chem.MolFromSmiles(product_smiles, sanitize=False)
    atom_list = get_atom_list(
        precursor_smiles, product_smiles, atom_environment=AtomEnvironment.CHANGED
    )

    # Set atoms in product with a different combing env to 1
    [
        atom.SetAtomMapNum(1)
        if atom.GetAtomMapNum() in atom_list
        else atom.SetAtomMapNum(0)
        for atom in products_mol.GetAtoms()
    ]

    return Chem.MolToSmiles(products_mol)


def find_number_tags(smiles: str) -> int:
    """
    Given a SMILES string determines how many tags are present by scanning the string.
    Only identifies tags in the form: [<atom>:1]

    Args:
        smiles: The SMILES string
    Returns:
        Number of tags (int) or 'NaN' if there was an issue processing the string.
    """
    try:
        pattern = r":1\]"
        return len(re.findall(pattern, smiles))
    except Exception:
        return 0


def get_changed_ids(mol: Chem.rdchem.Mol) -> List[int]:
    """
    Obtains a list of the changed atomIdxs from a tagged molecule

    Args:
        mol: RdKit molecule object
    Returns:
        List of changed atomIdxs
    """

    # Get tagged atoms
    changed_ids = []
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() != 0:
            changed_ids.append(atom.GetIdx())
        else:
            continue

    return changed_ids


def permute_tags(smiles: str) -> List[Tuple[int, ...]]:
    """
    Given a SMILES string determines how many tag combinations are possible.
    Only considers combinations upto 4 tags.

    Args:
        smiles: The SMILES string
    Returns:
        List of tag combinations (int) or 0 if there was an issue processing the string.
    """
    # Get tagged combinations
    tagging_combos: List[Tuple[int, ...]] = []
    try:
        mol = Chem.MolFromSmiles(smiles)

        # Get tagged atoms
        changed_ids = get_changed_ids(mol)

        if len(changed_ids) > 4:
            max_tags = 4
        else:
            max_tags = len(changed_ids)

        for n in range(1, max_tags):
            for i in combinations(changed_ids, n):
                tagging_combos.append(i)

        return tagging_combos
    except Exception:
        return tagging_combos


def get_tag_combinations_id_dict(
    mol: Chem.rdchem.Mol,
) -> Dict[
    int, List[Tuple[int, ...]]
]:  # tuple can contain variable number of ints min 1
    """
    From a molecule object with tagged atoms, obtain a dictionary containing
    possible atomIdx permutations of the tagged atoms.

    Args:
        mol: RdKit molecule object
    Returns:
        Dictionary containing possible atomIdx permutations of the tagged atoms.
    """

    # Get tagged atoms
    changed_ids = get_changed_ids(mol)

    if len(changed_ids) > 4:
        max_tags = 4
    else:
        max_tags = len(changed_ids)

    # Get tagged combinations
    tagging_combos = {}

    for n in range(1, max_tags):
        combos = []
        for i in combinations(changed_ids, n):
            combos.append(i)

        tagging_combos[n] = combos

    return tagging_combos


def sample_from_permutation_ids(
    tag_combination_id_dict: Dict[int, List[Tuple[int, ...]]],
    number_of_permutations: int = 1,
) -> List[Tuple[int, ...]]:
    """
    Given a dictionary containing the tag combinations, sample 'size' permutations

    Args:
        tag_combination_id_dict: Dictionary of {num_tags: [combinations]},
                                 from 'get_tag_combinations_id_dict'
        number_of_permutations: number of permutations to return
    Returns:
        List of Tuples of atomIdxs to permute for the given sample
    """
    permutation_ids = []

    if len(tag_combination_id_dict) != 0:
        for max_tag, id_list in tag_combination_id_dict.items():
            permutation = random.sample(id_list, number_of_permutations)
            permutation_ids.append(permutation)

    return list(chain(*permutation_ids))


def return_tag_combinations(smiles: str) -> int:
    """
    Given a SMILES string determines how many tag combinations are possible.
    Only considers combinations upto 4 tags.

    Args:
        smiles: The SMILES string
    Returns:
        Number of tag combinations (int) or 0 if there was an issue processing the
        string.
    """
    tagging_combos = permute_tags(smiles)
    if isinstance(tagging_combos, list):
        return len(tagging_combos)
    else:
        return 0


def permute_tagged_mol(
    mol: Chem.rdchem.Mol, permutation_ids: List[Tuple[int, ...]]
) -> List[str]:
    """
    From an rdkit molecule and corresponding ids to permute,
    return a permuted tagged SMILES

    Args:
        mol: RdKit molecule object
        permutation_ids: List of Tuples of atomIdxs to permute for the given sample
    Returns:
        SMILES of a permuted tagged molecule
    """
    permuted_tagged_smiles = []
    for ids in permutation_ids:
        for atom in mol.GetAtoms():
            if atom.GetIdx() in ids:
                atom.SetAtomMapNum(1)
            else:
                atom.SetAtomMapNum(0)

        permuted_tagged_smiles.append(Chem.MolToSmiles(mol))

    return permuted_tagged_smiles


def permute_tagged_smiles(smiles: str, number_of_permutations: int = 1) -> List[str]:
    """
    Permute the atom tags of a tagged molecule

    Args:
        smiles: Tagged SMILES string
        number_of_permutations: number of permutations to obtain
    Returns:
        List of tag permuted SMILES strings
    """
    mol = Chem.MolFromSmiles(smiles)
    tag_combinations = get_tag_combinations_id_dict(mol)
    permutation_ids = sample_from_permutation_ids(
        tag_combinations, number_of_permutations=number_of_permutations
    )
    permuted_smiles = permute_tagged_mol(mol, permutation_ids)

    return permuted_smiles
