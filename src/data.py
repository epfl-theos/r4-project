from dataclasses import dataclass
from pathlib import Path

import click

__all__ = [
    "DATA_3DCD",
    "DATA_MP",
    "R4Data",
]


@dataclass
class R4Data:

    structures: Path
    data: Path
    geo: Path
    soap: Path
    soap_red: Path
    chem: Path
    chem_red: Path

    @classmethod
    def from_prefix(cls, prefix: Path):
        return cls(
            structures=prefix / "structures.xyz",
            data=prefix / "data.zip",
            geo=prefix / "geo.npz",
            soap=prefix / "soap.npz",
            soap_red=prefix / "soap_red.npz",
            chem=prefix / "chem.json.gz",
            chem_red=prefix / "chem_red.json.gz",
        )


if Path("../r4data").exists():
    print(click.style("PRODUCTION DATA PRESENT", fg="green"))
    DATA_3DCD = R4Data.from_prefix(Path("../r4data/3DCD").resolve())
    DATA_MP = R4Data.from_prefix(Path("../r4data/MP").resolve())
    if not DATA_3DCD.structures.exists():
        print(
            click.style(
                "But cannot retrieve 3DCD .xyz structural files because the "
                "database is not public yet. To get access to the files, "
                "please email Simon Adorf at simon.adorf@epfl.ch. Using "
                "example MP data instead.",
                fg="yellow",
            ),
        )
        DATA_3DCD = R4Data.from_prefix(Path("../example_data/").resolve())
else:
    print(
        click.style(
            "WARNING: PRODUCTION DATA NOT PRESENT, USING EXAMPLE DATA FROM "
            "REDUCED MP DATA SET!",
            fg="yellow",
        )
    )
    DATA_3DCD = R4Data.from_prefix(Path("../example_data/").resolve())
    DATA_MP = R4Data.from_prefix(Path("../example_data/").resolve())
    DATA_3DCD.soap_red = DATA_3DCD.soap
    DATA_MP.soap_red = DATA_MP.soap
