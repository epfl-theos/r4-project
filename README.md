# Project: Rule of Four

## About

This is code and documentation relating to the investigation of the "Rule of Four".
The large files containing SOAP (Smooth Overlap of Atomic Positions) vectors and databases'visualizations can be found on the Materials Cloud archive.

## Repository organization

- `src/` - Source code.
    - `src/0.Rule4_hist` - Generation of histograms to represent the Rule of Four.
    - `src/1.SOAP_vectors` - Demo of SOAP vectors generation. No need to run it (if you do with the complete datasets, you might fall into RAM issues), as all the files named \*SOAP\*.npz contain the averaged (per structure) SOAP vectors already (n_rows=n_frames, n_columns~10^5).
    - `src/2. Global geometric descriptors` - Calculation of inherited symmetries (point groups) and packing parameters.
    - `src/3. Energies_PcovR` - Principal Covariates Regression analysis regressed on formation energies per atom.
    - `src/4. Classification_PCovR` - Principal Covariates Regression analysis regressed on Random Forest classification probabilities.
    - `data.py` - data class to run noteboks with production or with example data.
    - `modules.ipynb` - python modules that jupyter notebooks run.
- `example_data/` - A reduced set of the [Materials Project](https://materialsproject.org/) database to import in the jupyter notebooks. The .json files can be directly run on [Chemiscope](https://chemiscope.org/) (color according to the 'magic probability' label and with the 'seismic' palette).

## How to run the code

The dependencies for the source code are speciefied in the `requirements.txt` file.
The easiest way to create a computational environment to run the provided source code is to use [repo2docker](https://repo2docker.readthedocs.io/):

To run with the reduced data from the `example_data/` repository:

```console
repo2docker https://github.com/epfl-theos/r4-project.git
```
To run with complete data from the [Materials Cloud Archive](https://archive.materialscloud.org/) repository, fork the whole current repository (r4-project/) on your local machine, download the Materials Cloud data folder and rename it 'r4-data'. Then run the following:

```console
repo2docker --volume=/abs/path/to/r4data:r4data/ -E r4-project/.
```
This should open a jupyter notebook on the web where all dependencies are already installed, and the data folders can be accessed by changing the paths in the single notebooks.

## How to commit to this repository

Before making commits, make sure to install the pre-commit hooks by running the following commands within the repository root directory:

```console
pip install pre-commit==2.15.0
pre-commit install
```

This will run the pre-commit hooks on each commit ensuring that certain code standards are upheld.

## Acknowledgments

This project is led by Elena Gazzarrini and managed by Simon Adorf.
Supported by Rose K. Cersonsky and Marnik Bercx.

 - MARVEL INSPIRE FELLOWSHIP (Elena)
 - H2020 MARKETPLACE (Simon)
