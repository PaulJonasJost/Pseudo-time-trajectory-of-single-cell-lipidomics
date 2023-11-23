"""
This file contains all function that will be used to analyze the influence
of the number of labels on the optimization/parameter estimation.
"""
import os
import sbml_model

import pandas as pd

from typing import Sequence
from sbml_model_create import create_petab_model

def create_different_labels_models(
        n_labels: Sequence[int],
        parametersets: int = 100
):
    """
    Create models with different numbers of labels.

    Parameters
    ----------
    n_labels:
        The number of labels for each species.
    model_name:
        The name of the model.
    parametersets:
        The number of parametersets to create for each model.
    """
    # load all parameters
    df = pd.read_csv('parameter_set.csv', sep=",", index_col=0)
    par_set_dict = df.transpose().to_dict()
    # for each parameter set and each number of labels, we create a petabmodel
    for i in range(parametersets):
        for n in n_labels:
            labels = [f"16_{i}" for i in range(n)]
            labels.append("ul")
            sbml_model.LABELS = labels
            for label in labels:
                sbml_model.PARAMETERS[f'k_fa_{label}_influx'] = 0.01
            create_petab_model(
                index=i,
                parameter_dict=par_set_dict[i],
                noise=True,
                timepoints=[90],
                labels=labels
            )
        print("Done with parameterset", i)

if __name__ == '__main__':
    create_different_labels_models(n_labels=[1, 2, 3, 4, 5])
