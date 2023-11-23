"""Create a set of parameters around a mean a save them as a CSV file."""
import pandas as pd
import numpy as np

from sbml_model import PARAMETERS


def create_parameters(n_samples: int = 1000, standard_dev: float = 0.1):
    """
    Creates a parameter set with n_samples samples and standard deviation
    standard_dev.

    :param standard_dev: standard deviation of the normal distribution
    :param n_samples: number of samples
    """
    parameters = pd.DataFrame()
    for parameter, value in PARAMETERS.items():
        parameters[parameter] = pd.Series(
            [np.log10(value)] * n_samples  # log parameters
        )
    for parameter in parameters.columns:
        parameters[parameter] = parameters[parameter] + np.random.normal(
            0, 0.1, n_samples  # same standard deviation for all
        )
        # print the maximum and minimum value of the parameter
        print(
            f"Parameter {parameter} has min {min(parameters[parameter])} "
            f"and max {max(parameters[parameter])}"
        )

    # parameters need to be brought back in linear scale
    parameters = 10**parameters

    parameters.to_csv("parameter_set.csv")


if __name__ == "__main__":
    create_parameters()
