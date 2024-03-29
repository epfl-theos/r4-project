from dataclasses import dataclass
from pathlib import Path

import click

__all__ = [
    "DATA_MC3D",
    "DATA_MP",
    "R4Data",
]


@dataclass
class R4Data:
    structures: Path  # cif-structure files
    data: Path  # only for the MP (contains energies)
    geo: Path  # geometric npz, distribution of radii, alpha-parameters, etc.
    soap: Path  # soap vectors
    # soap_red: Path  # soap vectors (reduced feature selction, FPS)
    chem: Path  # chemiscope file
    # chem_red: Path  # chemiscope file (random sub-selection of 30,000)

    @classmethod
    def from_prefix(cls, prefix: Path):
        return cls(
            structures=prefix / "structures.xyz",
            data=prefix / "data.zip",
            geo=prefix / "geo.npz",
            # soap=prefix / "power_spectrum_pca_30-species.npy",
            soap=prefix / "soap.npz",
            # soap_red=prefix / "soap_red.npz",
            chem=prefix / "chem.json.gz",
            # chem_red=prefix / "chem_red.json.gz",
        )


if Path("../../r4data").exists():
    print(click.style("PRODUCTION DATA PRESENT", fg="green"))
    DATA_MC3D = R4Data.from_prefix(Path("../../r4data/MC3D").resolve())
    DATA_MP = R4Data.from_prefix(Path("../../r4data/MP").resolve())
    if not DATA_MC3D.structures.exists():
        print(
            click.style(
                "But cannot retrieve MC3D .xyz structural files because the "
                "experimental structures can not be released due to licensing constraints. Using "
                "example MP data instead.",
                fg="yellow",
            ),
        )
        DATA_MC3D = R4Data.from_prefix(Path("../example_data/").resolve())
else:
    print(
        click.style(
            "WARNING: PRODUCTION DATA NOT PRESENT, USING EXAMPLE DATA FROM "
            "REDUCED MP DATA SET!",
            fg="yellow",
        )
    )
    DATA_MC3D = R4Data.from_prefix(Path("../example_data/").resolve())
    DATA_MP = R4Data.from_prefix(Path("../example_data/").resolve())
    DATA_MC3D.soap_red = DATA_MC3D.soap
    DATA_MP.soap_red = DATA_MP.soap
