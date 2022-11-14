import logging
from typing import List

from rxnmapper import RXNMapper  # noqa

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def map_reactions(mapper: RXNMapper, reactions: List[str]) -> List[str]:
    """Map multiple reaction SMILES.
    This function may raise exceptions, typically if the number of tokens is
    larger than 512."""
    chunk_results = mapper.get_attention_guided_atom_maps(
        reactions, canonicalize_rxns=True
    )
    resulting_smiles = [result["mapped_rxn"] for result in chunk_results]
    return resulting_smiles


def map_reactions_with_error_handling(
    mapper: RXNMapper, reactions: List[str]
) -> List[str]:
    """
    Map multiple reaction SMILES.
    When there is an error, the reactions will be mapped one by one, and the
    one causing the error will be replaced by an empty reaction, ">>".
    """
    try:
        return map_reactions(mapper, reactions)
    except Exception:
        logger.warning(
            f"Error while mapping chunk of {len(reactions)} reactions. "
            "Mapping them individually."
        )

    mapped_reactions = []
    for reaction in reactions:
        try:
            mapped_reaction = map_reactions(mapper, [reaction])[0]
        except Exception as e:
            logger.info(
                f"Reaction causing the error: {reaction}; "
                f"{e.__class__.__name__}: {e}"
            )
            mapped_reaction = ">>"
        mapped_reactions.append(mapped_reaction)
    return mapped_reactions
