import datetime
import os
import re

from amici.petab_simulate import PetabSimulator
# # For newer versions of AMICI, use the following import instead:
# from amici.petab.simulator import PetabSimulator
import libsbml
import pandas as pd
import petab
import pypesto.petab
import simplesbml

from sbml_model import PARAMETERS, SPECIES, LABELS, REACTIONS, T_PULSE
from sbml_model_util import (
    create_labeled_species,
    pulse_chase_w_labels,
)


# a list of all states that are not observables
NON_OBSERVABLES = [
    "CDP_C",
    "CDP_EA",
    "Choline",
    "Ethanolamine",
    "G3P",
    "GLY",
    "ph_C",
    "ph_EA",
    "Serine",
]


def create_synthetic_data(
        petab_problem, filename: str, model=None, noise=True
):
    """
    Create synthetic data for a given PEtab problem and save it to a file.

    This function simulates data using the provided PEtab problem and saves it to a file.
    If noise is True, random noise will be added to the simulated data. The function
    ensures that no negative values are present in the synthetic data by re-simulating
    if necessary.

    Parameters
    ----------
    petab_problem : petab.Problem
        The PEtab problem to simulate data for.
    filename : str
        The base filename to save the data to. The function will append '.tsv' to this name.
    model : amici.Model, optional
        The AMICI model to use for simulation. If provided, it can speed up compilation.
        Default is None.
    noise : bool, optional
        Whether to add noise to the simulated data. Default is True.

    Returns
    -------
    None
        The function saves the data to a file but does not return any value.
    """
    solver = None
    if model is not None:
        # Create solver instance
        solver = model.getSolver()
        solver.setRelativeTolerance(1e-12)
        solver.setAbsoluteTolerance(1e-15)

    simulator = PetabSimulator(
        petab_problem, amici_model=model
    )

    any_negative = True
    while any_negative:
        petab_problem.measurement_df = simulator.simulate(
            noise=noise,
            model_name=None,
            force_compile=False,
            solver=solver,
            as_measurement=True,
        )
        any_negative = any(petab_problem.measurement_df["measurement"] < 0)
        if any_negative:
            print("Negative values in synthetic data. Re-simulating.")
    if not noise:
        petab_problem.measurement_df["noiseParameters"] = (
            petab_problem.measurement_df["measurement"] / 10
        )
    if not noise:
        petab_problem.measurement_df.to_csv(
            f"{filename}_wo_noise.tsv", sep="\t"
        )
    petab_problem.measurement_df.to_csv(f"{filename}.tsv", sep="\t")


def create_sbml_model(parameter_dict: dict = PARAMETERS, index: int = None):
    """
    Create the SBML model of the lipidomics model.

    This function creates an SBML model with the specified parameters. It adds all
    parameters, species, and reactions to the model, and sets up pulse-chase labeling.
    The model is then saved to an XML file.

    Parameters
    ----------
    parameter_dict : dict, optional
        A dictionary of parameter names and values to use in the model.
        Default is the PARAMETERS constant from sbml_model.py.
    index : int, optional
        If provided, the model will be saved with this index in its filename and directory.
        Default is None.

    Returns
    -------
    str
        The base filename of the saved SBML model.
    """
    # MODEL
    model = simplesbml.SbmlModel(level=2, version=4)

    for parameter, value in parameter_dict.items():
        model.addParameter(parameter, value)
    for species, n_labels in SPECIES.items():
        create_labeled_species(
            species_name=species, n_labels=n_labels, model=model, labels=LABELS
        )
    for reaction in REACTIONS:
        reaction.create_reaction(labels=LABELS, model=model)
    pulse_chase_w_labels(
        labels=LABELS,
        model=model,
        t_pulse=T_PULSE,
        # t_between=T_BETWEEN
    )

    # write model to file
    # time = datetime.datetime.now().strftime("%Y_%m_%d")
    sbml_file = f"lipidomics_2023_08_30"
    if index is not None:
        dir = f"../Petab_models_230829/3_labels/{sbml_file}_{index}"
        sbml_file = f"{sbml_file}_{index}"
    else:
        dir = f"../../Petab_models/{sbml_file}"
    if not os.path.isdir(dir):
        os.makedirs(dir)

    libsbml.writeSBMLToFile(model.document, f"{dir}/{sbml_file}.xml")

    # Read in the file
    with open(f"{dir}/{sbml_file}.xml", "r") as file:
        filedata = file.read()

    # Replace the target string
    filedata = re.sub("<delay>.*?</delay>\n", "", filedata, flags=re.DOTALL)

    # Write the file out again
    with open(f"{dir}/{sbml_file}.xml", "w") as file:
        file.write(filedata)

    return f"{sbml_file}"


def create_petab_model(
    index: int = None,
    parameter_dict: dict = PARAMETERS,
    noise: bool = False,
    timepoints: list = [90],
):
    """
    Create a PEtab model from an SBML model.

    This function creates a PEtab model by first creating an SBML model using
    create_sbml_model, and then generating the necessary PEtab files (observable table,
    condition table, measurement table, and parameter table). It also creates synthetic
    data for the model if noise is True.

    Parameters
    ----------
    index : int, optional
        If provided, the model will be saved with this index in its filename and directory.
        Default is None.
    parameter_dict : dict, optional
        A dictionary of parameter names and values to use in the model.
        Default is the PARAMETERS constant from sbml_model.py.
    noise : bool, optional
        Whether to add noise to the synthetic data. Default is False.
    timepoints : list, optional
        The timepoints at which to simulate measurements. Default is [90].

    Returns
    -------
    None
        The function creates PEtab files but does not return any value.
    """
    sbml_file = create_sbml_model(parameter_dict=parameter_dict, index=index)

    dir = f"../Petab_models_230829/3_labels/{sbml_file}"
    if not os.path.isdir(dir):
        os.makedirs(dir)
    sbml_reader = libsbml.SBMLReader()
    sbml_doc = sbml_reader.readSBML(f"{dir}/{sbml_file}.xml")
    sbml_model = sbml_doc.getModel()

    # create observable table
    observables = sbml_model.getListOfSpecies()
    observable_df = petab.create_observable_df()
    for observable in observables:
        if observable.getId() in NON_OBSERVABLES:
            continue
        id_obs = f"observable_{observable.getId()}"
        new_entry = pd.DataFrame.from_dict(
            {
                "observableId": [id_obs],
                "observableName": [id_obs],
                "observableFormula": [f"{observable.getId()}"],
                "noiseFormula": [f"noiseParameter1_{id_obs}"],
                "noiseDistribution": ["normal"],
            }
        )
        observable_df = pd.concat([observable_df, new_entry])
    observable_df = observable_df.set_index("observableId")
    # create condition table
    condition_df = petab.create_condition_df(
        parameter_ids=["conditionName", "preeq_switch"],
        condition_ids=["synthetic", "preeq"],
    )
    condition_df.loc["synthetic"] = pd.Series(
        {"conditionName": "synthetic", "preeq_switch": 1}
    )
    condition_df.loc["preeq"] = pd.Series(
        {"conditionName": "preeq", "preeq_switch": 0}
    )
    # create measurement table
    measurement_df = petab.create_measurement_df()
    for observable in observable_df.index:
        preeq_ID = "preeq"
        cartesian_product = pd.core.reshape.util.cartesian_product
        # create multi timepoint measurement table
        obs, preeq, simCond, time, meas, sigma = cartesian_product(
            [[observable], [preeq_ID], ["synthetic"], timepoints, [0], [0.01]]
        )
        new_entry = pd.DataFrame(
            data={
                "observableId": obs,
                "preequilibrationConditionId": preeq,
                "simulationConditionId": simCond,
                "time": time,
                "measurement": meas,
                "datasetId": simCond,
                "noiseParameters": sigma,
            }
        )
        measurement_df = pd.concat([measurement_df, new_entry])
    # create parameter table from sbml model
    parameter_df = petab.create_parameter_df(
        sbml_model=sbml_model,
        observable_df=observable_df,
        condition_df=condition_df,
        measurement_df=measurement_df,
        include_optional=True,
        lower_bound=1e-04,
        upper_bound=1e2,
    )

    petab.write_observable_df(
        df=observable_df, filename=f"{dir}/observables.tsv"
    )
    petab.write_condition_df(
        df=condition_df, filename=f"{dir}/conditions.tsv"
    )
    petab.write_measurement_df(
        df=measurement_df, filename=f"{dir}/measurements.tsv"
    )
    petab.write_parameter_df(
        df=parameter_df, filename=f"{dir}/parameters.tsv"
    )

    petab.create_problem_yaml(
        sbml_files=f"{dir}/{sbml_file}.xml",
        condition_files=f"{dir}/conditions.tsv",
        measurement_files=f"{dir}/measurements.tsv",
        parameter_file=f"{dir}/parameters.tsv",
        observable_files=f"{dir}/observables.tsv",
        yaml_file=f"{dir}/{sbml_file}.yaml",
    )

    # create synthetic data
    petab_problem = petab.Problem.from_yaml(f"{dir}/{sbml_file}.yaml")
    importer = pypesto.petab.PetabImporter(
        petab_problem=petab_problem,
        output_folder="./amici_models/lipidomics_2022_09_12",
    )
    problem = importer.create_problem(
        n_threads=8,
        force_compile=False,
    )
    create_synthetic_data(
        petab_problem=petab_problem,
        filename=f"{dir}/measurements",
        model=problem.objective.amici_model,
        noise=False,
    )
    if noise:
        create_synthetic_data(
            petab_problem=petab_problem,
            filename=f"{dir}/measurements",
            model=problem.objective.amici_model,
            noise=True,
        )


if __name__ == "__main__":
    df = pd.read_csv("parameter_set.csv", sep=",", index_col=0)
    par_set_dict = df.transpose().to_dict()
    for i_index, index in enumerate(par_set_dict):
        create_petab_model(
            index=index,
            parameter_dict=par_set_dict[index],
            noise=True,
            timepoints=[90],
        )
