"""
Simulate a test case with a given force field.

This runs a short NVT simulation and saves the trajectory.
It saves:
- initial.pdb: the initial minimized structure
- minimized.pdb: the minimized structure
- trajectory.pdb: the trajectory of the simulation
"""

import pathlib
import sys
import click

import tqdm
import openmm
import openmm.app
from loguru import logger
from openff.toolkit import Molecule, ForceField

logger.remove()
logger.add(sys.stdout)

def create_simulation(
    interchange,
    pdb_stride: int = 500,
    trajectory_name: str = "trajectory.pdb",
) -> openmm.app.Simulation:
    """Copied from Interchange docs"""
    integrator = openmm.LangevinMiddleIntegrator(
        300 * openmm.unit.kelvin,
        1 / openmm.unit.picosecond,
        1 * openmm.unit.femtoseconds,
    )

    simulation = interchange.to_openmm_simulation(
        combine_nonbonded_forces=True,
        integrator=integrator,
    )

    # https://github.com/openmm/openmm/issues/3736#issuecomment-1217250635
    simulation.minimizeEnergy()

    simulation.context.setVelocitiesToTemperature(300 * openmm.unit.kelvin)
    simulation.context.computeVirtualSites()

    trajectory_name = str(trajectory_name)
    pdb_reporter = openmm.app.PDBReporter(trajectory_name, pdb_stride)
    simulation.reporters.append(pdb_reporter)

    return simulation



@click.command(help=__doc__)
@click.option(
    "--force-field",
    "-ff",
    type=str,
    default="../../forcefields/fb-fit-v3-single-mean-k100.offxml",
    help="The force field to use.",
    show_default=True,
)
@click.option(
    "--output-directory",
    "-o",
    "output_directory",
    type=str,
    default="output/v3-k100",
    help="The output directory to write files to.", 
    show_default=True,
)
def main(
    force_field: str = "../../forcefields/fb-fit-v3-single-mean-k100.offxml",
    output_directory: str = "output/v3-k100",
    initial_file: str = "initial.pdb",
    minimized_file: str = "minimized.pdb",
    trajectory_file: str = "trajectory.pdb",
):
    forcefield = ForceField(force_field)

    mol = Molecule.from_smiles(
        "C3=NC1=C(NC2=C1NC=C2)S3"
    )
    mol.generate_conformers(n_conformers=1)
    
    # set up output options
    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    initial_file = output_directory / initial_file
    minimized_file = output_directory / minimized_file
    trajectory_file = output_directory / trajectory_file

    ic = forcefield.create_interchange(mol.to_topology())
    ic.to_pdb(initial_file)
    logger.info(f"Wrote initial structure to {initial_file}")

    # minimize
    ic.minimize()
    ic.to_pdb(minimized_file)
    logger.info(f"Wrote minimized structure to {minimized_file}")

    simulation = create_simulation(
        ic,
        trajectory_name=trajectory_file,
    )
    n_steps = 10000
    logger.info(f"Running {n_steps} steps of dynamics")
    for _ in tqdm.tqdm(range(n_steps), desc="Simulating"):
        simulation.step(1)
    logger.info(f"Wrote trajectory to {output_directory / trajectory_file}")


if __name__ == "__main__":
    main()
