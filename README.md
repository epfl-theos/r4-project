# Project. The rule of four: anomalous stoichiometries of inorganic compounds

## About

This is code and documentation relating to the investigation:

*The rule of four: anomalous stoichiometries of inorganic compounds.*

The large files containing SOAP (Smooth Overlap of Atomic Positions) vectors and databases' visualizations can be found on the [Materials Cloud Archive](https://archive.materialscloud.org/record/2023.104).

## Repository organization

- `src/` - Source code.
    - `src/0.Rule4_hist.ipynb` - Generation of histograms to represent the Rule of Four.
    - `src/1.SOAP_vectors.ipynb` - Demo of SOAP vectors generation. No need to run it (if you do with the complete datasets, you might have RAM issues), as all the files named \*SOAP\*.npz contain the averaged (per structure) SOAP vectors already (n_rows=n_frames, n_columns~10^5).
    - `src/2. Global geometric descriptors.ipynb` - Calculation of inherited symmetries (point groups) and packing parameters.
    - `src/3. Energies_PcovR.ipynb` - Principal Covariates Regression analysis regressed on formation energies per atom.
    - `src/4. Classification_RF_species-invariant.ipynb` - Principal Covariates Regression analysis regressed on Random Forest classification probabilities.
    - `src/5. Classification-comparison.ipynb` - Comparison of efficiency of various classification algorithms.
    - `data.py` - data class to run noteboks with production or with example data.
    - `modules.ipynb` - python modules that jupyter notebooks run.
    - `contants.py`, `point_groups.json` and `utils.py` are functions to run the point group analysis in notebook 2. 
    - `mp_symms` and `exp_symms` are being used to prove the RoF existence by testing how the primitive unit cells are generated. This refers to the unpublished .xyz structure data in the case of the MC3D database. 
    - `test_indices_all` and `train_indices_all` refer to the last part of the classification analyses in notebooks 4 and 5. 
- `example_data/` - A reduced set of the [Materials Project](https://materialsproject.org/) database to import in the jupyter notebooks. The .json files can be directly run on [Chemiscope](https://chemiscope.org/) (color according to the 'magic probability' label and with the 'seismic' palette).

## How to run the code

The dependencies for the source code are specified in the `requirements.txt` and in the `environment.yml` files.

The easiest way to create a computational environment to run the provided source code is to use [repo2docker](https://repo2docker.readthedocs.io/):

To run with the reduced data from the `example_data/` repository:

```console
repo2docker https://github.com/epfl-theos/r4-project.git
```
To run with complete data from the [Materials Cloud Archive](https://archive.materialscloud.org/record/2023.104) submission, clone the whole current repository (r4-project/) on your local machine, download the Materials Cloud data folder and rename it 'r4-data'. Then run the following:

```console
repo2docker --volume=/abs/path/to/r4data:r4data/ -E r4-project/.
```
This should open a jupyter notebook on the web where all dependencies are already installed, and the data folders can be accessed by changing the paths in the single notebooks.


## Reading the pre-existing `.npz` files 

The `.npz` files contain information about each structure and can be opened with a dataclass, as specified in `src/data.py`. 

```
@dataclass
class R4Data:

    structures: Path  # cif-structure files
    data: Path  # only for the MP (contains energies)
    geo: Path  # geometric npz, distribution of radii, alpha-parameters, etc.
    soap: Path  # soap vectors
    chem: Path  # chemiscope file
    @classmethod
    def from_prefix(cls, prefix: Path):
        return cls(
            structures=prefix / "structures.xyz",
            data=prefix / "data.zip",
            geo=prefix / "geo.npz", 
            soap=prefix / "soap.npz",           
            chem=prefix / "chem.json.gz"
        )
```

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
