# Imports
import click
import pandas as pd
from rxn.chemutils.reaction_smiles import (
    ReactionFormat,
    parse_any_reaction_smiles,
    to_reaction_smiles,
)
from rxn.utilities.containers import chunker
from rxnmapper import RXNMapper

from dar.aam import map_reactions_with_error_handling


@click.command()
@click.option(
    "--file",
    type=str,
    required=True,
    help="Enter the absolute path to the file containing rxns",
)
@click.option(
    "--reaction_column",
    type=str,
    default="rxn",
    help="Enter the name of the column containing reactions",
)
def map_reactions(file: str, reaction_column: str):
    """
    Maps reactions from a file containing unmapped reaction SMILES

    Args:
        file (str): Absolute path to a file containing reaction SMILES, one on each line

    Returns:
        A file with the original unmapped rection, the atom-mapped reaction
    """
    df = pd.read_csv(file)
    rxns = df[reaction_column].to_list()
    df[reaction_column] = [
        to_reaction_smiles(parse_any_reaction_smiles(rxn), ReactionFormat.STANDARD)
        for rxn in rxns
    ]

    # Remap reactions using the predicted precursors and apply retagging
    # Used to determine which bond was broken
    rxn_mapper = RXNMapper()

    results = []
    batch_size = 64
    for chunk in chunker(df[reaction_column].to_list(), batch_size):
        results.extend(map_reactions_with_error_handling(rxn_mapper, chunk))

    df["mapped_rxn"] = results

    df.to_csv(file.split(".")[0] + ".mapped.csv", index=False)


if __name__ == "__main__":
    map_reactions()
