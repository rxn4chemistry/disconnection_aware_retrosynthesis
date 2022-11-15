import click
import pandas as pd
from rxn.chemutils.utils import remove_atom_mapping

from dar.chem import remove_unmapped_components, standardise_reaction_component
from dar.tagging import find_number_tags, get_tagged_products, return_tag_combinations


@click.command()
@click.option(
    "--file_path",
    "-f",
    help="Absolute path to input file output from rxn_reaction_preprocessing",
)
@click.option(
    "--remove_unmapped",
    "-r",
    default=False,
    help="Removes unmapped species from reactants.",
)
def analyse_tags(
    file_path: str, remove_unmapped: bool, extract_templates: bool
) -> None:
    """
    Given an output file from 'rxn_reaction_preprocessing' containing atom-mapped reaction SMILES:
    - Generated Tags corresponding to changed atoms (disconnections)
    - Removes unmapped scpecies (optional)
    - Calculates the number of tagged atoms
    - Calculates the number of possible tag permutations
    - Generates a report on the distribution of atom tags in the dataset
    - Filters the dataset for:
        - Reactions with 0 tagged atoms (i.e. no disconnection, no change in atom environments)
        - Reactions with >10 tagged atoms (too many bond changes for a reaction to reasonably occur, low frequency)

    Args:
        file (str): Absolute path to the output file from 'rxn_reaction_preprocessing' containing atom-mapped reaction SMILES

    Returns:
        A filtered csv file with the tagged products, optionally removal of unmapped species
        csv report of atom tag distribution across the dataset
    """
    print("Reading Data....")
    data = pd.read_csv(file_path)
    data[["reactants", "products"]] = data["mapped_rxn"].str.split(">>", expand=True)

    print("Tagging Products....")
    data["tagged_products"] = data.apply(
        lambda x: get_tagged_products(x["reactants"], x["products"]), axis=1
    )

    if remove_unmapped:
        data["reactants"] = data["reactants"].apply(remove_unmapped_components)

    data["reactants"] = data["reactants"].apply(remove_atom_mapping)
    data["reactants"] = data["reactants"].apply(standardise_reaction_component)

    print("Calculating Num Tags....")
    data["num_tags"] = data["tagged_products"].apply(find_number_tags)
    data.to_csv(file_path.split(".")[0] + ".tagged.csv", header=True, index=False)

    print("Calculating Num Tag Combinations....")
    data["tag_combinations"] = data["tagged_products"].apply(return_tag_combinations)

    data.to_csv(file_path.split(".")[0] + ".tagged.csv", header=True, index=False)

    tag_stats = data.groupby(["num_tags"]).size()
    tag_stats_df = data.groupby(["num_tags"]).size()[:10].to_frame(name="counts")
    tag_stats_df.reset_index(inplace=True)
    counts_gte_ten = tag_stats[10:].sum()
    tag_stats_df.loc[len(tag_stats_df.index)] = ["10+", counts_gte_ten]
    tag_stats_df["dataset_fraction"] = tag_stats_df["counts"].apply(
        lambda x: round(x / tag_stats_df["counts"].sum() * 100, 2)
    )

    tag_stats_df.to_csv(
        file_path.split(".")[0] + ".tagged_stats.csv", header=True, index=False
    )

    data = data[(data.num_tags != 0) & (data.num_tags <= 10)].copy()

    data["tagged_rxn"] = data[["reactants", "tagged_products"]].agg(">>".join, axis=1)

    data.to_csv(
        file_path.split(".")[0] + ".tagged_filtered.csv", header=True, index=False
    )


if __name__ == "__main__":
    analyse_tags()
