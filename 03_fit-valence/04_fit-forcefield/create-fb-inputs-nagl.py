"""
This scripts generates the ForceBalance inputs for fitting a force field.

It generates the following files:
- optimize.in: The main ForceBalance input file.
- targets/: the training targets for the optimization.
- forcefield/force-field.offxml: The force field to optimize.
"""

import os
import json
import pathlib
import typing
from loguru import logger

import click

from openff.toolkit import Molecule, ForceField
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper

from openff.qcsubmit.results import (
    OptimizationResultCollection,
    TorsionDriveResultCollection,
)
import numpy as np
import tqdm

# logging.basicConfig(
#     level=logging.INFO,
#     format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
# )
# logger = logging.getLogger(__name__)

def filter_for_smarts_or_smiles(
    entry,
    smarts_to_exclude: typing.Optional[typing.List[str]] = None,
    inchi_keys_to_exclude: typing.Optional[typing.List[str]] = None,
) -> bool:
    """
    Normal filtering with the to_records() call is incredibly slow;
    copy out the actual filtering function to speed it up

    Parameters
    ----------
    entry : OptimizationResult or TorsionDriveResult
        The entry to filter.
    smarts_to_exclude : list[str], optional
        A list of SMARTS patterns to exclude from the training set.
    inchi_keys_to_exclude : list[str], optional
        A list of InChI keys to exclude from the training set.
    
    Returns
    -------
    bool
        True if the entry should be included, False if it should be excluded.
    """
    from openff.toolkit import Molecule

    mol = Molecule.from_mapped_smiles(entry.cmiles, allow_undefined_stereo=True)

    if smarts_to_exclude:
        for smarts in smarts_to_exclude:
            if mol.chemical_environment_matches(smarts):
                return False
    
    if inchi_keys_to_exclude:
        if mol.to_inchikey(fixed_hydrogens=True) in inchi_keys_to_exclude:
            return False
    return True
        

def filter_dataset(
    dataset,
    smarts_to_exclude: typing.Optional[typing.List[str]] = None,
    smiles_to_exclude: typing.Optional[typing.List[str]] = None,
):
    """
    Filter the dataset to exclude certain SMILES and SMARTS patterns.
    This is a side-effecting function that modifies the dataset in place.

    Parameters
    ----------
    dataset : OptimizationResultCollection or TorsionDriveResultCollection
        The dataset to filter.
    smarts_to_exclude : list[str], optional
        A list of SMARTS patterns to exclude from the training set.
    smiles_to_exclude : list[str], optional
        A list of SMILES patterns to exclude from the training set.
        This is a list of SMILES strings, not patterns.
    
    Returns
    -------
    None
        The dataset is modified in place.
    """
    from openff.toolkit import Molecule

    inchi_keys_to_exclude = []
    if smiles_to_exclude:
        inchi_keys_to_exclude = [
            Molecule.from_smiles(smiles).to_inchikey(fixed_hydrogens=True)
            for smiles in smiles_to_exclude
        ]
    
    key = list(dataset.entries)[0]
    original_dataset = dataset.entries[key]
    filtered = [
        entry
        for entry in original_dataset
        if filter_for_smarts_or_smiles(
            entry,
            smarts_to_exclude=smarts_to_exclude,
            inchi_keys_to_exclude=inchi_keys_to_exclude,
        )
    ]
    dataset.entries[key] = filtered



def load_training_data(
    optimization_dataset: str,
    torsion_dataset: str,
    smarts_to_exclude: typing.Optional[str] = None,
    smiles_to_exclude: typing.Optional[str] = None,
    verbose: bool = False
):
    """
    Load the training data from the given datasets and filter out unwanted entries.
    Parameters
    ----------
    optimization_dataset : str
        The path to the optimization dataset to use, encoding an OptimizationResultCollection.
    torsion_dataset : str
        The path to the torsion dataset to use, encoding a TorsionDriveResultCollection
    smarts_to_exclude : str, optional
        The path to a file containing SMARTS patterns to exclude from the training set.
        The patterns should be separated by new lines.
    smiles_to_exclude : str, optional
        The path to a file containing SMILES patterns to exclude from the training set.
        The patterns should be separated by new lines.
    verbose : bool, optional
        Whether to print verbose logging messages.
    """
    if smarts_to_exclude is not None:
        exclude_smarts = pathlib.Path(smarts_to_exclude).read_text().splitlines()
    else:
        exclude_smarts = []

    if smiles_to_exclude is not None:
        exclude_smiles = pathlib.Path(smiles_to_exclude).read_text().splitlines()
    else:
        exclude_smiles = []

    torsion_training_set = TorsionDriveResultCollection.parse_file(torsion_dataset)
    if verbose:
        logger.info(f"Loaded torsion training set with {torsion_training_set.n_results} entries.")

    filter_dataset(
        torsion_training_set,
        smarts_to_exclude=exclude_smarts,
        smiles_to_exclude=exclude_smiles,
    )

    if verbose:
        logger.info(f"Filtered torsion training set to {torsion_training_set.n_results} entries.")

    optimization_training_set = OptimizationResultCollection.parse_file(optimization_dataset)
    if verbose:
        logger.info(f"Loaded optimization training set with {optimization_training_set.n_results} entries.")

    filter_dataset(
        optimization_training_set,
        smarts_to_exclude=exclude_smarts,
        smiles_to_exclude=exclude_smiles,
    )
    if verbose:
        logger.info(f"Filtered optimization training set to {optimization_training_set.n_results} entries.")


    return torsion_training_set, optimization_training_set


@click.command()
@click.option(
    "--tag",
    type=str,
    default="fb-fit",
    help="The tag to use for the fitting run.",
)
@click.option(
    "--optimization-dataset",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the optimization dataset to use. (JSON)",
)
@click.option(
    "--torsion-dataset",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the torsion dataset to use. (JSON)",
)
@click.option(
    "--forcefield",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the force field to use. (offxml)",
)
@click.option(
    "--valence-counts",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help=(
        "The path to the valence counts file. (JSON)"
        "This is used to determine which valence parameters to train."
    ),
)
@click.option(
    "--torsion-counts",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help=(
        "The path to the torsion counts file. (JSON)"
        "This is used to determine which torsion parameters to train."
    ),
)
@click.option(
    "--n-min-valence",
    type=int,
    default=5,
    show_default=True,
    help="The minimum number of optimizations with a parameter to train it.",
)
@click.option(
    "--n-min-torsion",
    type=int,
    default=1,
    show_default=True,
    help=(
        "The minimum number of torsiondrive records with a parameter to train it. "
        "Note, this only counts parameter IDs that match the full dihedral, "
        "not all parameter IDs that run through a torsion central bond."
    )
)
@click.option(
    "--output-directory",
    type=click.Path(exists=False, dir_okay=True, file_okay=False),
    required=True,
    help="The directory to write the ForceBalance inputs to.",
)
@click.option(
    "--frozen-angle-file",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    required=False,
    help="The path to a JSON file containing frozen angle smirks (e.g. for linear angles)",
)
@click.option(
    "--smarts-to-exclude",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    default=None,
    help=(
        "The path to a file containing a list of SMARTS patterns "
        "to exclude from the training set. "
        "The patterns should be separated by new lines."
    ),
)
@click.option(
    "--smiles-to-exclude",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    default=None,
    help=(
        "The path to a file containing a list of SMILES patterns "
        "to exclude from the training set. "
        "The patterns should be separated by new lines."
    ),
)
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    help="Whether to print verbose logging messages.",
)
@click.option(
    "--max-iterations",
    type=int,
    default=50,
    show_default=True,
    help="The maximum number of iterations to run the fitting for.",
)
@click.option(
    "--port",
    type=int,
    default=55387,
    show_default=True,
    help="The port to run the server on.",
)
@click.option(
    "--angle-k-prior",
    type=float,
    default=20.0,
    show_default=True,
    help="Angle k prior",
)
@click.option(
    "--bond-k-prior",
    type=float,
    default=20.0,
    show_default=True,
    help="Bond k prior",
)
@click.option(
    "--torsion-k-prior",
    type=float,
    default=5.0,
    show_default=True,
    help="Torsion k prior",
)
@click.option(
    "--angle-angle-prior",
    type=float,
    default=1.0,
    show_default=True,
    help="Angle angle prior",
)
@click.option(
    "--bond-length-prior",
    type=float,
    default=0.01,
    show_default=True,
    help="Bond length prior",
)
def generate(
    tag: str,
    optimization_dataset: str,
    torsion_dataset: str,
    forcefield: str,
    valence_counts: str,
    torsion_counts: str,
    n_min_valence: str,
    n_min_torsion: str,
    output_directory: str,
    frozen_angle_file: typing.Optional[str] = None,
    smarts_to_exclude: typing.Optional[str] = None,
    smiles_to_exclude: typing.Optional[str] = None,
    verbose: bool = False,
    max_iterations: int = 50,
    port: int = 55387,
    angle_angle_prior: float = 1.0,
    bond_length_prior: float = 0.01,
    angle_k_prior: float = 20.0,
    bond_k_prior: float = 20.0,
    torsion_k_prior: float = 5.0

):
    from openff.toolkit import ForceField
    from openff.bespokefit.optimizers.forcebalance import ForceBalanceInputFactory
    from openff.bespokefit.schema.fitting import OptimizationSchema, OptimizationStageSchema
    from openff.bespokefit.schema.optimizers import ForceBalanceSchema
    from openff.bespokefit.schema.targets import (
        OptGeoTargetSchema,
        TorsionProfileTargetSchema,
    )
    from openff.bespokefit.schema.smirnoff import (
        AngleHyperparameters,
        AngleSMIRKS,
        BondHyperparameters,
        ImproperTorsionHyperparameters,
        ProperTorsionHyperparameters,
        BondSMIRKS,
        ProperTorsionSMIRKS,
        ImproperTorsionSMIRKS,
    )

    torsion_training_set, optimization_training_set = load_training_data(
        optimization_dataset=optimization_dataset,
        torsion_dataset=torsion_dataset,
        smarts_to_exclude=smarts_to_exclude,
        smiles_to_exclude=smiles_to_exclude,
        verbose=verbose
    )

    ff = ForceField(forcefield)

    # check correct charge model
    assert "ToolkitAM1BCC" not in ff._parameter_handlers
    assert "ChargeIncrementModel" in ff._parameter_handlers
    charge_handler = ff.get_parameter_handler("ChargeIncrementModel")
    assert charge_handler.partial_charge_method == "openff-gnn-am1bcc-0.1.0-rc.3.pt"
    
    # remove constraints
    constraint_handler = ff.get_parameter_handler("Constraints")
    if len(constraint_handler.parameters) == 3:
        constraint_handler._parameters.pop(0) # this is the h-bond constraint

    # set up optimizer schema
    optimizer = ForceBalanceSchema(
        max_iterations=max_iterations,
        step_convergence_threshold=0.01,
        objective_convergence_threshold=0.1,
        gradient_convergence_threshold=0.1,
        n_criteria=2,
        initial_trust_radius=-1.0,
        finite_difference_h=0.01,
        extras={
            "wq_port": str(port),
            "asynchronous": "True",
            "search_tolerance": "0.1",
            "backup": "0",
            "retain_micro_outputs": "0",
        },
    )

    # set up target schema
    targets = [
        TorsionProfileTargetSchema(
            reference_data=torsion_training_set,
            energy_denominator=1.0,
            energy_cutoff=8.0,
            extras={"remote": "1"},
        ),
        OptGeoTargetSchema(
            reference_data=optimization_training_set,
            weight=0.01,
            extras={"batch_size": 30, "remote": "1"},
            bond_denominator=0.05,
            angle_denominator=5.0,
            dihedral_denominator=10.0,
            improper_denominator=10.0,
        ),
    ]

    linear_angle_smirks = []
    if frozen_angle_file:
        with open(frozen_angle_file, "r") as f:
            linear_angle_smirks = json.load(f)["smirks"]
    logger.info(f"Frozen angle SMIRKS: {linear_angle_smirks}")
        
    # load counts to determine which parameters to train
    with open(valence_counts, "r") as f:
        valence_counts = json.load(f)
    with open(torsion_counts, "r") as f:
        torsion_counts = json.load(f)

    target_parameters = []

    # select angles to parameterize
    angle_handler = ff.get_parameter_handler("Angles")
    for parameter in angle_handler.parameters:
        attrs = {"k", "angle"}
        count = valence_counts.get(parameter.id, 0)
        if count >= n_min_valence:
            if parameter.smirks in linear_angle_smirks:
                attrs = {"k"}
            target_parameters.append(
                AngleSMIRKS(
                    smirks=parameter.smirks,
                    attributes=attrs,
                )
            )
        else:
            logger.info(
                f"Not training {parameter.id} with SMIRKS {parameter.smirks} "
                f"because it has {count} records < {n_min_valence}."
            )

    # select bonds to parameterize
    bond_handler = ff.get_parameter_handler("Bonds")
    for parameter in bond_handler.parameters:
        count = valence_counts.get(parameter.id, 0)
        if count >= n_min_valence:
            target_parameters.append(
                BondSMIRKS(
                    smirks=parameter.smirks,
                    attributes={"k", "length"},
                )
            )
        else:
            logger.info(
                f"Not training {parameter.id} with SMIRKS {parameter.smirks} "
                f"because it has {count} records < {n_min_valence}."
            )

    # torsions
    torsion_handler = ff.get_parameter_handler("ProperTorsions")

    for parameter in torsion_handler.parameters:
        count = torsion_counts.get(parameter.id, 0)
        if count >= n_min_torsion:
            original_k = parameter.k
            attributes = {f"k{i + 1}" for i in range(len(original_k))}
            target_parameters.append(
                ProperTorsionSMIRKS(
                    smirks=parameter.smirks,
                    attributes=attributes,
                )
            )
        else:
            logger.info(
                f"Not training {parameter.id} with SMIRKS {parameter.smirks} "
                f"because it has {count} records < {n_min_torsion}."
            )

    # re-fit impropers directly from force field
    improper_torsion_handler = ff.get_parameter_handler("ImproperTorsions")
    for improper in improper_torsion_handler.parameters:
        original_k = improper.k
        attributes = {f"k{i + 1}" for i in range(len(original_k))}
        target_parameters.append(
            ImproperTorsionSMIRKS(
                smirks=improper.smirks,
                attributes=attributes
            )
        )

    # set up the full optimization schema
    optimization_schema = OptimizationSchema(
        id=tag,
        initial_force_field=os.path.abspath(forcefield),
        stages=[
            OptimizationStageSchema(
                optimizer=optimizer,
                targets=targets,
                parameters=target_parameters,
                parameter_hyperparameters=[
                    AngleHyperparameters(priors={"k": angle_k_prior, "angle": angle_angle_prior}),
                    BondHyperparameters(priors={"k": bond_k_prior, "length": bond_length_prior}),
                    ProperTorsionHyperparameters(priors={"k": torsion_k_prior}),
                    ImproperTorsionHyperparameters(priors={"k": torsion_k_prior}),
                ],
            )
        ]
    )

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    # ensure we can re-create this
    optfile = f"{optimization_schema.id}.json"
    with open(optfile, "w") as f:
        f.write(optimization_schema.json(indent=2))
    

    # Generate the ForceBalance inputs
    ForceBalanceInputFactory.generate(
        os.path.join(optimization_schema.id),
        optimization_schema.stages[0],
        ff,
    )



if __name__ == "__main__":
    generate()
