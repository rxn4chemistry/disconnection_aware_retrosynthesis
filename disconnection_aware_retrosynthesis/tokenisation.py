from typing import Tuple

from rxn_biocatalysis_tools import (
    detokenize_enzymatic_reaction_smiles,
    tokenize_enzymatic_reaction_smiles,
)


def tokenize_and_split_enzymatic_reaction_smiles(rxn: str) -> Tuple[str, str]:
    """
    Given an enzymatic reaction SMILES tokenises and splits into reactants and products

    Args:
        rxn: Enzymatic reaction SMILES containg EC number
    Returns:
        Tuple of reactants, products
    """
    reactants, products = [
        x.strip() for x in tokenize_enzymatic_reaction_smiles(rxn).split(">>")
    ]
    return reactants, products


def detokenize_and_split_enzymatic_reaction_smiles(rxn: str) -> Tuple[str, str, str]:
    """
    Given an enzymatic reaction SMILES detokenises and splits into reactants
    and products

    Args:
        rxn: Enzymatic reaction SMILES containg EC number
    Returns:
        Tuple of reactants, ec, products
    """
    reactants, products = [
        x.strip() for x in detokenize_enzymatic_reaction_smiles(rxn).split(">>")
    ]
    print(reactants, flush=True)
    if "|" in reactants:
        reactants, ec = [x.strip() for x in reactants.split("|")]
    else:
        ec = None
    return reactants, ec, products
