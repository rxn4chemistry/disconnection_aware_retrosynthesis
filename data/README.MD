# Data Availability

The data used to train the models in this study can be found at:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7101695.svg)](https://doi.org/10.5281/zenodo.7101695)

## Usage of pre-processed data on Zenodo

- Download the data from Zenodo.
- You should have a file called 'complete_disconnection_labelled.csv' which you can use with our [notebook](../notebooks_and_scripts/basic_tag_and_tokenise.ipynb) to reprocess the data.
- Alternatively, you can retrain the models using the tokenised datasets.

## Usage from Lowe patent data

The disconnection aware retrosynthesis model may still be trained using the procedure outline in this repository starting from the US Patent Office extracts by Lowe,

- Download the [USPTO data](https://figshare.com/articles/dataset/Chemical_reactions_from_US_patents_1976-Sep2016_/5104873) from Figshare.
- Unzip '1976_Sep2016_USPTOgrants_smiles.7z' which we will use for demo purposes.

### Preprocessing uisng [rxn-reaction-preprocessing](https://github.com/rxn4chemistry/rxn-reaction-preprocessing) (recommended)

- [workflow](./notebooks_and_scripts/rxn_preprocessing_workflow.ipynb)

  For further information about reaction preprocessing refer to [rxn-reaction-preprocessing](https://github.com/rxn4chemistry/rxn-reaction-preprocessing)

### Basic Data Processing

- You should have a file called '1976_Sep2016_USPTOgrants_smiles.rsmi' which you can use with our [notebook](../notebooks_and_scripts/basic_preprocessing_example.ipynb) for basic data pre-processing.
