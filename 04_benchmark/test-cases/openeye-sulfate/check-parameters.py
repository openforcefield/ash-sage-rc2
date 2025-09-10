"""Check the parameters assigned to sulfate in a given force field."""
from openff.toolkit import ForceField, Molecule, unit

import click
from loguru import logger
import sys
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
IPythonConsole.ipython_useSVG=True 

logger.remove()
logger.add(sys.stdout)


@click.command(help=__doc__)
@click.option(
    "--force-field",
    "-ff",
    type=str,
    default="openff-2.2.1.offxml",
    help="Path to the force field file"
)
def main(
    force_field: str = "openff-2.2.1.offxml",
):
    ligand = Molecule.from_smiles("[O-]S(=O)(=O)[O-]")
    ff = ForceField(force_field)
    logger.info(f"Using force field {force_field}")
    angles = ff.label_molecules(ligand.to_topology())[0]["Angles"]
    for indices, parameter in angles.items():
        elements = [ligand.atoms[i].symbol for i in indices]
        label = "-".join(f"{el}{i + 1}" for el, i in zip(elements, indices))
        logger.info(f"{label}  |  {parameter.id} -- {parameter.angle}, {parameter.k}")

if __name__ == "__main__":
    main()