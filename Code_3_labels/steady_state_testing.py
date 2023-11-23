"""
Test whether the model has more than one steady state
"""
import logging
import typing

import numpy as np
import os
import pypesto
import petab
import pypesto.petab
import amici.plotting

logger = logging.getLogger(__name__)


file_dir = os.path.dirname(__file__)
model_name = 'lipidomics_2023_05_11'
petab_dir = f"{file_dir}/../Code_multilabel_comparison/Petab_models"
index = 0

# load the model
petab_problem = petab.Problem.from_yaml(
    f"{petab_dir}/{model_name}_3_labels_{index}/"
    f"{model_name}_3_labels_{index}.yaml"
)
importer = pypesto.petab.PetabImporter(
    petab_problem=petab_problem,
    output_folder=f'../Code/amici_models/{model_name}',
    model_name=model_name
)
problem = importer.create_problem(
    n_threads=8,
    guess_steadystate=True,
    force_compile=False,
)
obj = problem.objective
model = obj.amici_model
edata = obj.edatas[0]
edata.setTimepoints([0,])

# run amici
solver = model.getSolver()
solver.setNewtonMaxSteps(10)
solver.setSensitivityMethod(amici.SensitivityMethod.forward)
solver.setSensitivityOrder(amici.SensitivityOrder.first)
rdata = amici.runAmiciSimulation(model, solver, edata)

for key, value in rdata.items():
    if key == 'preeq_status':
        print('%20s: ' % key, value)
rdatas = [rdata]

# create 1000 random states
lb = np.array([0] * len(model.getInitialStates()))
ub = np.array([100] * len(model.getInitialStates()))
startpoints = pypesto.startpoint.uniform(1000, lb, ub)

for sp in startpoints:
    model.setInitialStates(sp)
    rdata = amici.runAmiciSimulation(model, solver, edata)
    for key, value in rdata.items():
        if key == 'preeq_status':
            print('%20s: ' % key, value)
    rdatas.append(rdata)

# get values of steady states except for failed ones
xss = [rdata.x for rdata in rdatas if not np.isnan(rdata.x).any()]
# get eucledean norm between points and general case
dx = [np.linalg.norm(x, rdatas[0].x) for x in xss]


# add any comparison you want, e.g. print(np.max(dx))

print("Done")