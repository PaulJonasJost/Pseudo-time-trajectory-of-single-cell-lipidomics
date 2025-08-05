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

from sbml_model import PARAMETERS, SPECIES, LABELS, REACTIONS, T_PULSE, T_TOTAL
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
    """create synthetic data for a given petab problem
    and save it to filename. Can also provide the amici MODEL to speed up
    compilation."""
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


def create_sbml_model(
    parameter_dict: dict = PARAMETERS, index: int = None, labels: list = LABELS
):
    """
    Creates the SBML model of the lipidomics model.
    """
    # MODEL
    model = simplesbml.SbmlModel(level=2, version=4)

    for parameter, value in parameter_dict.items():
        # if the parameter starts with "k_fa_",
        # only include it if it continues with a valid label
        if parameter.startswith("k_fa_"):
            if parameter == "k_fa_influx" or parameter == "k_fa_deg":
                model.addParameter(parameter, value)
            if not any(label in parameter for label in labels):
                continue
        model.addParameter(parameter, value)
    for species, n_labels in SPECIES.items():
        create_labeled_species(
            species_name=species, n_labels=n_labels, model=model, labels=labels
        )
    for reaction in REACTIONS:
        reaction.create_reaction(labels=labels, model=model)
    pulse_chase_w_labels(
        labels=labels, model=model, t_pulse=T_PULSE, t_total=T_TOTAL
    )

    # write model to file
    sbml_file = f"lipidomics_2023_08_30_{len(labels)}_labels"
    if index is not None:
        dir = f"../Petab_models_230829/multilabels/{sbml_file}_{index}"
        sbml_file = f"{sbml_file}_{index}"
    else:
        dir = f"Petab_models/{sbml_file}"
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
    labels: list = LABELS,
):
    sbml_file = create_sbml_model(
        parameter_dict=parameter_dict, index=index, labels=labels
    )

    dir = f"../Petab_models_230829/multilabels/{sbml_file}"
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
        output_folder=f"./amici_models/lipidomics_2023_04_19_"
                      f"{len(labels)}_labels",
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
    n_labels_max = 5
    for i_index, index in enumerate(par_set_dict):
        if i_index > 200:
            continue
        for i_labels in range(n_labels_max + 1):
            labels = [f"{i + 1}" for i in range(i_labels)]
            labels.append("ul")

            create_petab_model(
                index=index,
                parameter_dict=par_set_dict[index],
                noise=True,
                timepoints=[120],
                labels=labels,
            )
