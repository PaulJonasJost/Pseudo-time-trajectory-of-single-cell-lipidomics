"""Here we define functions to create figures for the paper."""

import pypesto
import pypesto.petab

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from sbml_model import PARAMETERS


def rank_cells_by_obj_val(n_cells: int = 1000):
    """
    Rank cells by objective value.

    Create a plot that shows the cells ranked by their objective value.
    x-axis: cell rank
    y-axis: objective value
    """
    objective_values = []
    indices = []
    for i in range(n_cells):
        if i == 27:
            continue
        # load results
        results = pypesto.store.read_result(
            f"../results_h5/lipidomics_2023_08_30_{i}_100starts.hdf5",
            optimize=True,
            problem=False
        )
        # get the best objective value
        objective_values.append(results.optimize_result[0].fval)
        indices.append(i)
        print(f"Finished cell {i}.")
    # sort the objective values
    objective_vals, ind = zip(*sorted(
        zip(objective_values, indices),
        reverse=True
    ))

    # create figure
    plt.figure(figsize=(6, 3))

    # plot the objective values as points
    plt.plot(sorted(ind, reverse=True), objective_vals, '-')
    # plot best point
    plt.plot(999, objective_vals[0], 'o')
    plt.xlabel('Cell rank')
    plt.ylabel('Objective value')
    plt.savefig('Model_230830/rank_cells_by_obj_val.pdf')


def plot_residues_all_cells(
        model_name,
        n_cells: int = 1000,
        parameter_ind: int = None
):
    """Plot the residuals of all cells."""
    def get_residual(index):
        """Get the residuals of a cell."""
        # create the problem
        yaml_file = f'../Petab_models_230829/3_labels/{model_name}_{index}/' \
                    f'{model_name}_{index}.yaml'
        importer = pypesto.petab.PetabImporter.from_yaml(
            yaml_file,
            output_folder=f'../Code/amici_models/{model_name}',
            model_name=model_name
        )
        problem = importer.create_problem(
            n_threads=8,
            guess_steadystate=True,
            # force_compile=True,
        )
        problem.objective.amici_solver.setRelativeTolerance(1e-12)  # 1e-8
        problem.objective.amici_solver.setAbsoluteTolerance(1e-15)  # 1e-12
        results = pypesto.store.read_result(
            f"../results_h5/lipidomics_2023_08_30_{index}_100starts.hdf5",
            optimize=True,
            problem=False
        )
        # call the objective function
        eval = problem.objective(
            x=results.optimize_result[0].x,
            return_dict=True,
            mode='mode_res',
        )
        # if objective values are not close, raise warning
        if not np.isclose(eval['fval'], results.optimize_result[0].fval):
            print(f"Warning: objective values are not close for cell {index}.")
        # return the residuals
        return eval['res']

    # get the residuals of all cells
    residuals = [get_residual(i) for i in range(n_cells)]
    residuals = np.array(residuals)
    if parameter_ind is None:
        residuals = residuals.flatten()
    else:
        residuals = residuals[:, parameter_ind]
    # plot the residuals as histogram
    sns.histplot(data=residuals, kde=True)  # TODO: residuals are
    # multidimensional
    from scipy import stats
    ks_result = stats.normaltest(residuals)
    print(ks_result)

    plt.savefig('Model_230830/all_residuals_all_cells.pdf')


def plot_par_scatter(n_cells:int=1000, par_ind: int=0):
    """Plot the histogram of the parameters of all cells."""
    # get the parameter name
    par_name = pypesto.store.read_result(
        f"../results_h5/lipidomics_2023_08_30_0_100starts.hdf5",
        optimize=True,
        problem=True,
    ).problem.x_names[par_ind]
    true_mean = PARAMETERS[par_name]
    def get_par(index, par_ind):
        """Get the parameters of a cell."""
        results = pypesto.store.read_result(
            f"../results_h5/lipidomics_2023_08_30_{index}_100starts.hdf5",
            optimize=True,
            problem=False,
            profile=False,
            sample=False,
        )
        return results.optimize_result[0].x[par_ind]

    # plt.rcParams['text.usetex'] = True

    # get the parameters of all cells
    pars = [get_par(i, par_ind) for i in range(n_cells) if i != 27]
    # plot the parameters as histogram
    sns.histplot(pars, kde=True)
    plt.axvline(
        np.log10(true_mean), color='g', linestyle='dashed', linewidth=1,
        label='True'
    )
    # also add the empirical mean
    plt.axvline(
        np.log10(np.mean(np.power(10,pars))), color='r', linestyle='dashed',
        linewidth=1,
        label='Empirical'
    )
    plt.axvline(
        np.median(pars), color='orange',
        linestyle='dashed',
        linewidth=1,
        label='Median'
    )
    plt.legend()
    plt.title(f'Parameter {par_name}')
    plt.savefig(f'Model_230830/{par_name}_hist.pdf')


if __name__ == '__main__':
    model_name = 'lipidomics_2023_08_30'
    # # create the figures
    # rank_cells_by_obj_val(n_cells=1000)
    # # clear the plots
    # plt.clf()
    for i in range(40):
        plot_par_scatter(n_cells=1000, par_ind=i)
        plt.clf()
    plot_residues_all_cells(
        model_name=model_name,
        n_cells=1000,
        # parameter_ind=0
    )