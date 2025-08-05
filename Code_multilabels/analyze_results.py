"""Code to analyze the results"""
import logging
import typing

import numpy as np
import os
import pypesto
import petab
import pypesto.petab
import pypesto.visualize as visualize
import matplotlib.pyplot as plt
import pypesto.visualize.model_fit as model_fit
import seaborn as sns
import amici.plotting

from sbml_model import PARAMETERS

from typing import Sequence, Union

logger = logging.getLogger(__name__)


class CumulativeResults:
    """
    Class to contain all results of multiple optimizations for further tests.
    """

    def __init__(
        self,
        model_name: str,
        res_dict: str,
        n_res: int = 1200,
        pypesto_problem: Union[pypesto.Problem, None] = None,
        n_starts: int = 1000,
    ):
        """
        Init function for class.

        :param model_name:
        :param res_dict:
        :param n_res:
        """
        self.model_name = model_name
        self.n_results = n_res
        self.n_starts = n_starts

        self.results = self.load_results(res_dict)
        self.pypesto_problem = pypesto_problem
        self.petab_problems = None
        self.problems_loaded = False

        self.total_runs = np.sum(
            [len(result.optimize_result) for result in self.results]
        )
        self.fvals_best = None
        self.non_existing = None

    def __getitem__(self, index):
        """Define `optimize_result[i]` to access the i-th result."""
        try:
            return self.results[index]
        except IndexError:
            raise IndexError(
                f"{index} out of range for optimize result of "
                f"length {len(self.results)}."
            )

    def load_results(self, res_dict: str):
        """Load all results."""
        # check whether file exists and if no, save the index
        non_existing = []
        for index in range(self.n_results):
            filename_hdf5 = f"{res_dict}/{self.model_name}_{index}" \
                  f"_{self.n_starts}starts.hdf5"
            if not os.path.isfile(filename_hdf5):
                non_existing.append(index)
        results = [
            pypesto.store.read_result(
                filename_hdf5,
                optimize=True,
                problem=True,
            )
            for index in range(self.n_results)
            if index not in non_existing
        ]
        print(f"Loaded {len(results)} results, {len(non_existing)} missing.")
        self.non_existing = non_existing
        return results

    def load_petab_problems(self, petab_dir: str):
        """Load all petab problems."""
        self.petab_problems = [
            petab.Problem.from_yaml(
                f"{petab_dir}/{self.model_name}_{index}/"
                f"{self.model_name}_{index}.yaml"
            )
            for index in range(self.n_results)
        ]
        if self.pypesto_problem is not None:
            self.problems_loaded = True

    def load_pypesto_problems(
        self, petab_problem: petab.Problem = None, force_compile: bool = False
    ):
        """Load the pypesto problem"""
        if petab_problem is None:
            if self.petab_problems is None:
                raise ValueError(
                    "No petab problems loaded yet, please provide one."
                )
            petab_problem = self.petab_problems[0]
            importer = pypesto.petab.PetabImporter(
                petab_problem=petab_problem,
                output_folder=f"../Code/amici_models/{self.model_name}",
                model_name=self.model_name,
            )
            problem = importer.create_problem(
                n_threads=8,
                guess_steadystate=True,
                force_compile=force_compile,
            )
            problem.objective.amici_solver.setRelativeTolerance(1e-12)  # 1e-8
            problem.objective.amici_solver.setAbsoluteTolerance(1e-15)  # 1e-12

            self.pypesto_problem = problem
            if self.petab_problems is not None:
                self.problems_loaded = True

    def fvals_overview(self, mode: Union[str, None] = None):
        """
        Print the mean of the best fvals and the sd.

        If mode is res, rather than the fvals, we will use the differences
        between the best fval and the fval with the true parameters.

        Parameters
        ----------
        mode:
            If mode=="res", will switch to differences of fval with references.
            Only works, if we already got the problems.
        """
        if self.fvals_best is None:
            self.fvals_best = [
                result.optimize_result.fval[0] for result in self.results
            ]
        if mode == "res":
            if not self.problems_loaded:
                raise ValueError(
                    "No problems loaded yet, not able to create "
                    "differences of the objective values."
                )
            fvals_res = [
                self.fvals_best[index]
                - self.pypesto_problem.objective(
                    self.petab_problems[index].x_nominal_scaled
                )
                for index in range(len(self.fvals_best))
            ]

        fval_mean = np.mean(self.fvals_best)
        fval_sd = np.std(self.fvals_best)

        output = (
            f"## Summary of the best starts per result\n"
            f"* Maximum value: {np.max(self.fvals_best)}\n"
            f"* Minimum value: {np.min(self.fvals_best)}\n"
            f"* Mean value: {fval_mean}\n"
            f"* Standard deviation: {fval_sd}\n"
        )
        if mode == "res":
            fval_res_mean = np.mean(fvals_res)
            fval_res_sd = np.std(fvals_res)
            output += (
                f"* Mean value of the differences: {fval_res_mean}\n"
                f"* Standard deviation of the differences: {fval_res_sd}\n"
            )
            output += (
                f"* Number of runs, that archieved better results than "
                f"the reference: {np.sum(np.array(fvals_res) < 0)}\n"
            )

    def overview_failed_starts(self):
        """Print the percentage of failed starts."""
        non_finite = np.sum(
            [
                np.logical_not(np.isfinite(result.optimize_result.fval)).sum()
                for result in self.results
            ]
        )
        print(
            f"Out of {self.total_runs} total starts across "
            f"{len(self.results)}, {non_finite} starts ended in non-finite "
            f"values. These are {non_finite/self.total_runs * 100}%."
        )

    def rank_cells_by_obj_value(
        self,
        figuresize: Sequence[int] = (4, 2),
        plot_points: bool = True,
        n_points: int = 3,
        color_points: Sequence[str] = None,
    ):
        """
        Rank the cells by the objective value.

        Parameters
        ----------
        figuresize:
            Size of the figure.
        plot_points:
            If True, will plot the points of the best starts.
        n_points:
            Number of points to plot.
        color_points:
            Colors of the points.
        """
        if self.fvals_best is None:
            self.fvals_best = [
                result.optimize_result.fval[0] for result in self.results
            ]
        if color_points is None:
            color_points = ["#007f5f", "#aacc00", "#fca311", "#17c3b2"]
        sorted_indices = np.round(np.argsort(self.fvals_best)).astype(int)
        fvals_sorted = np.sort(self.fvals_best)

        # create figure
        figure = plt.figure(figsize=figuresize)

        # plot the objective values as points
        plt.plot(fvals_sorted, "-", linewidth=2)
        plt.xlabel("Cell rank")
        plt.ylabel("Objective value")

        if plot_points:
            # get the points evenly spaced including first and last
            indices = np.round(
                np.linspace(0, len(fvals_sorted) - 1, n_points).astype(int)
            )
            for i in range(n_points):
                plt.plot(
                    indices[i],
                    fvals_sorted[indices[i]],
                    "o",
                    color=color_points[i],
                    label=f"Cell {indices[i]}",
                )

        print(sorted_indices[list(indices)])
        return figure

    def plot_time_trajectories_fit(
        self,
        cell_index: int = 0,
        add_cosmetics: bool = False,
        obs_name: str = None,
        obs_labels: int = None,
        timepoints: np.ndarray = None,
    ):
        """
        Plot the fit of the time trajectory of a secific cell.

        Parameters
        ----------
        cell_index:
            Index of the cell to plot.
        add_cosmetics:
            If True, will add cosmetics to the plot.
        obs_name:
            Name of the observable to plot.
        obs_labels:
            Number of labels of the observable to plot.
        timepoints:
            Timepoints to plot.
        """
        if obs_name is None or obs_labels is None:
            obs_name = "DAG"
            obs_labels = 2
        obs_names = [
            "observable_" + obs_name + obs_labels * "_16_1",
            "observable_" + obs_name + obs_labels * "_16_2",
            "observable_" + obs_name + obs_labels * "_17_1",
        ]
        print(obs_names)
        if timepoints is None:
            timepoints = np.linspace(0, 90, 1201)

        ax = model_fit.time_trajectory_model(
            result=self.results[cell_index],
            problem=self.pypesto_problem,
            timepoints=timepoints,
            state_names=[],
            state_ids=[],
            observable_ids=obs_names,
        )
        # petab parameters
        obj = self.pypesto_problem.objective.set_custom_timepoints(
            timepoints_global=timepoints
        )
        ret = obj(
            self.petab_problems[cell_index].x_nominal_free_scaled,
            mode="mode_fun",
            sensi_orders=(0,),
            return_dict=True,
        )
        obs_ind = [
            obj.amici_model.getObservableIds().index(obs_name) for
            obs_name in obs_names
        ]
        amici.plotting.plotObservableTrajectories(
            rdata=ret["rdatas"][0],
            observable_indices=obs_ind,  # corresponding to dags
            ax=ax,
            model=obj.amici_model,
        )

        label_colors = [
            (0.502, 0.733, 0.509),
            (0.944, 0.730, 0.418),
            (0.944, 0.558, 0.744),
        ]
        color_t0 = (0.647, 0.745, 0.858)

        # working on the legend
        handles, labels = ax.get_legend_handles_labels()
        for i_han, handle in enumerate(handles):
            if i_han < 3:
                labels[i_han] = f"${obs_name}_{(i_han + 1) % 4}$ fitted"
            else:
                labels[i_han] = f"${obs_name}_{(i_han + 1) % 4}$ syn. true"
            handle.set_color(label_colors[i_han % 3])
        for i_line, line in enumerate(ax.lines):
            if i_line < 3:  # skip the original lines
                line.set_ls("--")
                line.set_markevery([-1])
                line.set_marker("o")

        if add_cosmetics:
            # extend lines to -50 for T0
            for i_line, line in enumerate(ax.lines):
                xx = line.get_xdata()
                yy = line.get_ydata()
                xx = np.concatenate(([-50], xx))
                yy = np.concatenate(([yy[0]], yy))
                line.set_data(xx, yy)
                if i_line < 3:  # skip the original lines
                    line.set_ls("--")

            # cosmetics #
            #############
            xticks = [-40, 0, 30, 60, 90]
            xlabels = ["$T_0$", "$T_1$", "$T_2$", "$T_3$", "Measurement"]
            ax.set_xticks(xticks, labels=xlabels)
            # measurement
            ax.axvline(x=90, color="r", linestyle="dashed", linewidth=1)
            # T0
            ax.axvline(x=-40, color=color_t0, linestyle="dashed", linewidth=1)
            ax.get_xticklabels()[0].set_color(color_t0)
            for i_timepoint, color in enumerate(label_colors):
                ax.axvline(
                    x=i_timepoint * (120 / (self.n_labels - 1)),
                    color=color,
                    linestyle="dashed",
                    linewidth=1,
                )
                ax.get_xticklabels()[i_timepoint + 1].set_color(color)

            ax.get_xticklabels()[-1].set_color("red")

            ax.set_xlim((-50, 95))

        labeling = [item.get_text() for item in ax.get_yticklabels()]
        labeling[1] = "0"
        ax.set_yticklabels(labeling)
        plt.title("")
        plt.ylabel("")
        plt.xlabel("")
        ax.legend(
            handles=handles,
            labels=labels,
            loc="center left",
            bbox_to_anchor=(1, 0.5)
        )
        ax.figure.set_size_inches(5, 5)

        return ax

    def plot_waterfall(self, cell_index: int = 0, include_ref: bool = False):
        """
        Plot waterfall plot for a specific cell.

        Parameters
        ----------
        cell_index:
            Index of the cell to plot.
        """
        if include_ref:
            petab_problem = self.petab_problems[cell_index]
            importer = pypesto.petab.PetabImporter(
                petab_problem=petab_problem,
                output_folder=f"../Code/amici_models/{self.model_name}",
                model_name=self.model_name,
            )
            problem = importer.create_problem(
                n_threads=8,
                guess_steadystate=True,
            )
            problem.objective.amici_solver.setRelativeTolerance(1e-12)  # 1e-8
            problem.objective.amici_solver.setAbsoluteTolerance(1e-15)

        ax = visualize.waterfall(results=self.results[cell_index])
        return ax

    def scatter_plot_parameters(
        self,
        parameter_index: int,
        cell_indices: int = None,
        alpha: float = 0.5,
        return_ax: bool = True,
        return_pars: bool = False,
    ):
        """
        Create a scatter plot of the parameters. For all starts.

        Parameters
        ----------
        parameter_index:
            Indices of the parameters to plot. If None, all parameters are
            plotted. If it is an integer, the first `parameter_indices`
            parameters are plotted.
        cell_indices:
            Indices of the cells to plot.
        alpha:
            Alpha value for the scatter plot.
        return_ax:
            If True, the axis object is returned.
        return_pars:
            If True, the parameters are returned.
            Is ignored if return_ax == True
        """
        if cell_indices is None:
            cell_indices = len(self.petab_problems)

        parameters_est = [
            result.optimize_result[0].x for result in self.results
        ]
        parameters_true = [
            petab_problem.x_nominal_free_scaled
            for petab_problem in self.petab_problems
        ]

        # turn both into numpy arrays
        parameters_est = np.array(parameters_est)
        parameters_true = np.array(parameters_true)

        # select the parameters to plot
        parameters_est = parameters_est[:cell_indices, parameter_index]
        parameters_true = parameters_true[:cell_indices, parameter_index]

        fig, ax = plt.subplots()
        for i in range(parameters_est.shape[1]):
            ax.scatter(
                parameters_true[:, i],
                parameters_est[:, i],
                label=self.pypesto_problem.x_names[i],
                alpha=alpha,
            )
        ax.set_xlabel("True value")
        ax.set_ylabel("Estimated value")

        if return_ax:
            return ax
        if return_pars:
            return parameters_est, parameters_true

    def parameter_histogram(self, n_cells: int = 1000, par_ind: int = 0):
        """Plot the histogram of the parameters of all cells."""
        # get the parameter name
        fig, ax = plt.subplots()

        par_name = self.pypesto_problem.x_names[par_ind]
        true_mean = PARAMETERS[par_name]

        # get the parameters of all cells
        pars = [result.optimize_result[0].x[par_ind]
                for result in self.results]
        # plot the parameters as histogram
        sns.histplot(pars, kde=True)
        plt.axvline(
            np.log10(true_mean),
            color="g",
            linestyle="dashed",
            linewidth=1,
            label="True",
        )
        # also add the empirical mean
        plt.axvline(
            np.log10(np.mean(np.power(10, pars))),
            color="r",
            linestyle="dashed",
            linewidth=1,
            label="Empirical",
        )
        plt.axvline(
            np.median(pars),
            color="orange",
            linestyle="dashed",
            linewidth=1,
            label="Median",
        )

        plt.title(f"Parameter {par_name}")
        ax.figure.set_size_inches(4.5, 3)

        return ax

    def box_plot_parameters(self, parameter_indices=None):
        """
        Plot a Boxplots of all parameters.
        """
        if parameter_indices is None:
            parameter_indices = range(
                len(self.petab_problems[0].x_nominal_free_scaled)
            )

        parameters_est = [
            result.optimize_result[0].x for result in self.results
        ]
        # turn into numpy arrays
        parameters_est = np.array(parameters_est)
        parameters_est = [parameters_est[:, i] for i in parameter_indices]

        positions_est = [1 + i * 3 for i in range(len(parameter_indices))]

        parameters_true = [
            petab_problem.x_nominal_free_scaled
            for petab_problem in self.petab_problems
        ]
        parameters_true = np.array(parameters_true)
        positions_true = [2 + i * 3 for i in range(len(parameter_indices))]
        parameters_true = [parameters_true[:, i] for i in parameter_indices]

        fig, ax = plt.subplots(figsize=(50, 10))
        basewidth = 0.9
        patches = []
        labels = []
        for percentiles, width, alpha in zip(
            [(5, 95), (15, 85), (25, 75)],
            [basewidth * 0.5, basewidth * 0.75, basewidth],
            ["88", "BB", "FF"],
        ):
            boxes_true = [
                compute_percentiles(
                    parameter_true, percentiles=percentiles, label="True"
                )
                for parameter_true in parameters_true
            ]
            boxes_est = [
                compute_percentiles(
                    parameter_est,
                    percentiles=percentiles,
                    label="Est"
                )
                for parameter_est in parameters_est
            ]
            patch, label = draw_boxplot(
                boxes_true,
                positions=positions_true,
                width=width,
                ax=ax,
                color=f"#89B110{alpha}",
            )
            patches.append(patch)
            labels.append(label)
            patch, label = draw_boxplot(
                boxes_est,
                positions=positions_est,
                width=width,
                ax=ax,
                color=f"#118ab2{alpha}",
            )
            patches.append(patch)
            labels.append(label)
        ax.yaxis.grid(
            True,
            linestyle="-",
            which="major",
            color="lightgrey",
            alpha=0.5
        )
        ax.xaxis.grid(False)
        patches.reverse()
        labels.reverse()

        num_boxes = len(parameter_indices) * 2

        ax.set_xlim(-1, num_boxes * 1.5 + 1)

        return ax


class CumulativeResultsLabels(CumulativeResults):
    """
    An addition to the CumulativeResults function to also encompass the
    multiple label structure.
    """

    def __init__(
        self,
        model_name: str,
        res_dict: str,
        n_res: int = 1200,
        n_labels: int = 6,
        offset: int = 2,
        pypesto_problems: Union[
            Sequence[pypesto.Problem], pypesto.Problem, None
        ] = None,
    ):
        self.n_labels = n_labels
        self.n_res_ind = int(n_res / self.n_labels)
        self.offset = offset
        self.model_name = model_name
        self.n_results = n_res
        self.res_id, self.res_labels = self.assign_result_infos()
        self.filenames = None

        self.results = self.load_results(res_dict)
        self.pypesto_problems = pypesto_problems
        self.petab_problems = None
        self.problems_loaded = False

        self.total_runs = np.sum(
            [len(result.optimize_result) for result in self.results]
        )
        self.fvals_best = None

    def load_results(self, res_dict: str):
        """Load all results."""
        filenames = []
        to_delete = []
        # check first that all files exist
        for index in range(self.n_results):
            filename = (
                f"{self.model_name}_"
                f"{int(index) // self.n_res_ind + self.offset}"
                f"_labels_{int(index) % 200}"
            )
            if not os.path.exists(f"{res_dict}/{filename}_1000starts.hdf5"):
                logger.warning(f"{filename} was not found.")
                to_delete.append(index)
                continue
            # append to filenames
            filenames.append(filename)
        self.filenames = filenames
        self.res_id = np.delete(self.res_id, to_delete)
        self.res_labels = np.delete(self.res_labels, to_delete)

        results = [
            pypesto.store.read_result(
                f"{res_dict}/{fn}_1000starts.hdf5", optimize=True, problem=True
            )
            for fn in filenames
        ]
        return results

    def load_petab_problems(self, petab_dir: str):
        """Load all petab problems."""
        self.petab_problems = [
            petab.Problem.from_yaml(f"{petab_dir}/{filename}/{filename}.yaml")
            for filename in self.filenames
        ]
        if self.pypesto_problems is not None:
            self.problems_loaded = True

    def load_pypesto_problems(
        self,
        petab_problems: Union[Sequence[petab.Problem], None] = None,
        force_compile: bool = False,
    ):
        """
        Load the pypesto problems for each label
        :param petab_problems:
            The petab problems from which to make the pypesto problems. If
            none, corresponding ones in the class are used.
        :param force_compile:
            Whether to force the recompilation of the problems.
        """
        if petab_problems is None:
            if self.petab_problems is None:
                raise ValueError(
                    "No petab problems loaded yet, please provide one."
                )
            petab_problems = [
                self.petab_problems[i * self.n_res_ind]
                for i in range(self.n_labels)
            ]
            index = [i * self.n_res_ind for i in range(self.n_labels)]
        pypesto_problems = []
        for petab_problem, inde in zip(petab_problems, index):
            importer = pypesto.petab.PetabImporter(petab_problem=petab_problem)
            print(f"CREATING MODEL FROM INDEX {inde}")
            problem = importer.create_problem(
                # n_threads=8,
                guess_steadystate=True,
                force_compile=True,
            )
            problem.objective.amici_solver.setRelativeTolerance(1e-12)  # 1e-8
            problem.objective.amici_solver.setAbsoluteTolerance(1e-15)  # 1e-12
            pypesto_problems.append(problem)

        self.pypesto_problems = pypesto_problems
        if self.petab_problems is not None:
            self.problems_loaded = True

    def assign_result_infos(self):
        """
        Assign ids and number of labels to each result
        """
        res_id = [i for i in range(self.n_res_ind)] * self.n_labels
        res_labels = [
            [i + self.offset] * self.n_res_ind for i in range(self.n_labels)
        ]
        res_labels_flat = [item for sublist in res_labels for item in sublist]
        return np.array(res_id), np.array(res_labels_flat)

    def lineplot_error_label(self, parameter_indices):
        """
        Plot the mean error for the given parameters per label as a line.
        :return: axis object
        """
        parameters_est = [
            result.optimize_result[0].x for result in self.results
        ]
        parameters_true = [
            petab_problem.x_nominal_free_scaled
            for petab_problem in self.petab_problems
        ]

        # turn both into numpy arrays
        parameters_est = np.array(parameters_est)
        parameters_true = np.array(parameters_true)

        # put them into categories
        par_err = np.abs(parameters_true - parameters_est)
        # get the errors for each label
        par_errors = [
            par_err[i * self.n_res_ind: (i + 1) * self.n_res_ind, :]
            for i in range(self.n_labels)
        ]
        par_means = np.array(
            [np.mean(par_error, axis=0) for par_error in par_errors]  # axis=1?
        )
        fig, ax = plt.subplots()
        for parameters_index in parameter_indices:
            ax.plot(
                par_means[:, parameters_index],
                "--o",
                label=f"Parameter {parameters_index}",
            )
        # change y axis to log scale
        ax.set_yscale("log")
        ax.legend()

    def box_plot_parameter(self, par_ind: int):
        """
        Draw a boxplot for the given parameter for each label.

        :param par_ind:
        :return: Axis Object
        """

        parameters_est = [
            result.optimize_result[0].x for result in self.results
        ]
        parameters_true = [
            petab_problem.x_nominal_free_scaled
            for petab_problem in self.petab_problems
        ]

        # turn both into numpy arrays
        parameters_est = np.array(parameters_est)
        parameters_true = np.array(parameters_true)

        # put them into categories
        par_err = np.abs(parameters_true - parameters_est)
        # get the errors for each label
        par_errors = [
            par_err[i * self.n_res_ind: (i + 1) * self.n_res_ind, :]
            for i in range(self.n_labels)
        ]
        par_i = np.array([par_err[:, par_ind] for par_err in par_errors])
        fig, ax = plt.subplots()
        ax.boxplot(par_i)
        return ax

    def box_plot_parameters(
        self, parameter_indices=None, n_labels: int = 0, whiskers: bool = False
    ):
        """

        :param parameter_indices:
        :param n_labels:
        :param whiskers:
        :return:
        """
        if parameter_indices is None:
            parameter_indices = range(
                len(self.petab_problems[0].x_nominal_free_scaled)
            )

        parameters_est = [
            result.optimize_result[0].x for result in self.results
        ]
        # turn into numpy arrays
        parameters_est = np.array(parameters_est)
        # cut down to only include label
        parameters_est = np.array([
            parameters_est[
                n_labels * self.n_res_ind: (n_labels + 1) * self.n_res_ind, :
            ]
        ])
        parameters_est = [parameters_est[:, i] for i in parameter_indices]

        positions_est = [1 + i * 3 for i in range(len(parameter_indices))]

        parameters_true = [
            petab_problem.x_nominal_free_scaled
            for petab_problem in self.petab_problems
        ]
        parameters_true = np.array(parameters_true)
        # cut down to only include label
        parameters_true = np.array([
            parameters_true[
                n_labels * self.n_res_ind: (n_labels + 1) * self.n_res_ind, :
            ]
        ])
        positions_true = [2 + i * 3 for i in range(len(parameter_indices))]
        parameters_true = [parameters_true[:, i] for i in parameter_indices]

        fig, ax = plt.subplots(figsize=(50, 10))
        basewidth = 0.9
        patches = []
        labels = []
        if not whiskers:
            for percentiles, width, alpha in zip(
                [(5, 95), (15, 85), (25, 75)],
                [basewidth * 0.5, basewidth * 0.75, basewidth],
                ["88", "BB", "FF"],
            ):
                boxes_true = [
                    compute_percentiles(
                        parameter_true, percentiles=percentiles, label="True"
                    )
                    for parameter_true in parameters_true
                ]
                boxes_est = [
                    compute_percentiles(
                        parameter_est, percentiles=percentiles, label="Est"
                    )
                    for parameter_est in parameters_est
                ]
                patch, label = draw_boxplot(
                    boxes_true,
                    positions=positions_true,
                    width=width,
                    ax=ax,
                    color=f"#89B110{alpha}",
                )
                patches.append(patch)
                labels.append(label)
                patch, label = draw_boxplot(
                    boxes_est,
                    positions=positions_est,
                    width=width,
                    ax=ax,
                    color=f"#118ab2{alpha}",
                )
                patches.append(patch)
                labels.append(label)
            ax.yaxis.grid(
                True,
                linestyle="-",
                which="major",
                color="lightgrey",
                alpha=0.5
            )
            ax.xaxis.grid(False)
            patches.reverse()
            labels.reverse()
        else:
            ax.boxplot(parameters_true, positions=positions_true)
            ax.boxplot(parameters_est, positions=positions_est)
        ax.legend(patches, labels, loc=9)

        # Now fill the boxes with desired colors
        num_boxes = len(parameter_indices) * 2

        ax.set_xlim(-1, num_boxes * 1.5 + 1)

        return ax

    def draw_error_boxplots(
        self, parameter_indices=None, n_labels: int = 0, absolute: bool = False
    ):
        """

        :param parameter_indices:
        :param n_labels:
        :param absolute:
        :return:
        """
        if parameter_indices is None:
            parameter_indices = range(
                len(self.petab_problems[0].x_nominal_free_scaled)
            )

        parameters_est = np.array(
            [result.optimize_result[0].x for result in self.results]
        )
        parameters_true = np.array(
            [
                petab_problem.x_nominal_free_scaled
                for petab_problem in self.petab_problems
            ]
        )
        parameters_err = parameters_true - parameters_est
        if absolute:
            parameters_err = np.abs(parameters_err)
        # cut down to only include label
        parameters_err = np.array([
            parameters_err[
                n_labels * self.n_res_ind: (n_labels + 1) * self.n_res_ind, :
            ]
        ])
        parameters_err = [parameters_err[:, i] for i in parameter_indices]

        positions_err = [1 + i * 2 for i in range(len(parameter_indices))]

        fig, ax = plt.subplots(figsize=(50, 10))
        basewidth = 0.9
        patches = []
        labels = []
        for percentiles, width, alpha in zip(
            [(5, 95), (15, 85), (25, 75)],
            [basewidth * 0.5, basewidth * 0.75, basewidth],
            ["88", "BB", "FF"],
        ):
            boxes_true = [
                compute_percentiles(
                    parameter,
                    percentiles=percentiles,
                    label="Error"
                )
                for parameter in parameters_err
            ]
            patch, label = draw_boxplot(
                boxes_true,
                positions=positions_err,
                width=width,
                ax=ax,
                color=f"#89B110{alpha}",
            )
            patches.append(patch)
            labels.append(label)
        ax.yaxis.grid(
            True,
            linestyle="-",
            which="major",
            color="lightgrey",
            alpha=0.5
        )
        ax.xaxis.grid(False)
        patches.reverse()
        labels.reverse()
        num_boxes = len(parameter_indices) * 2

        ax.set_xlim(-1, num_boxes * 1.5 + 1)

        return ax

    def calculate_fim_eigenvals(self):
        """Calculate the FIM eigenvalues for each best start."""
        eigenvals = []
        errors = []
        for i_label in range(self.n_labels):
            # get the results that have i_label + 1 labels
            indices = np.where(self.res_labels == i_label + self.offset)[0]
            eigenval_label = []
            for ind in indices:
                try:
                    eigenval_label.append(
                        calculate_FIM_eigenvalues(
                            self.results[ind],
                            self.pypesto_problems[i_label].objective
                        )
                    )
                except:
                    errors.append(f"{i_label}\t\t\t{ind}... failed")
            eigenvals.append(np.array([eigenval_label]).flatten())
        return eigenvals

    def draw_fim_0_values(
        self,
        cutoff=1e-10,
        visualization: str = "boxplot",
        scale: str = "log",
        eigenvals: list = None,
    ):
        """
        Draw Boxplots of number of 0 eigen_vals next to each other for each
        label-
        :param cutoff:
            determines which parameters are considered "0".
        :param visualization:
            "boxplot" or "hist"
        :param scale:
            "log" or "linear"
        :param eigenvals:
            if None, calculate eigenvals

        :return:
        """
        fig, ax = plt.subplots(figsize=(50, 10))
        colors = ["r", "b", "g", "o", "y"]
        if eigenvals is None:
            eigenvals = self.calculate_fim_eigenvals()
        if visualization == "boxplot":
            ax.boxplot(eigenvals)
            ax.set_yscale(scale)
        if visualization == "hist":
            for i_eig, eig in enumerate(eigenvals):
                ax.hist(eig, alpha=0.3, color=colors[i_eig])
            ax.set_yscale(scale)
        return ax


def compute_percentiles(
    data, label: str, percentiles: typing.Tuple[int] = (25, 75)
) -> dict:
    """
    Computes the percentiles of the data and returns a matplotlib box.

    :param data:
    :param label:
    :param percentiles:
    :return: dict
    """
    q1 = np.percentile(data, percentiles[0])
    q3 = np.percentile(data, percentiles[1])
    med = np.percentile(data, 50)

    return {
        "med": med,
        "q1": q1,
        "q3": q3,
        "whislo": q1,
        "whishi": q3,
        "label": label,
        "percentiles": percentiles,
    }


def draw_boxplot(
    boxes: typing.Sequence[dict],
    positions: int,
    width: float,
    ax: plt.Axes,
    color: str
):
    box = ax.bxp(
        boxes,
        positions=positions,
        showfliers=False,
        widths=width,
        boxprops=dict(color=color, facecolor=color),
        medianprops=dict(linewidth=2),
        capprops=dict(color=color),
        patch_artist=True,
    )
    label = f"Estimated {boxes[0]['percentiles']} percentiles"
    if boxes[0]["label"] == "True":
        label = f"True {boxes[0]['percentiles']} percentiles"
    return box["boxes"][0], label


def calculate_FIM_eigenvalues(result, objective):
    """
    Calculate the eigenvalue spectrum of the FIM at the best point
    :param result:
    :return:
    """
    # get the FIM
    fim_matrix = objective(result.optimize_result[0].x, sensi_orders=(2,))

    return np.linalg.eig(fim_matrix)[0]


def experimental_setup(result, problem, n_labels):
    def _assign_label(i):
        val = "ul" if i == 0 else f"{i +1}"
        return val
    obs_names = [
        f"observable_DAG_{_assign_label(i_label)}_{_assign_label(i_label)}"
        for i_label in range(n_labels)
    ]
    ax = model_fit.time_trajectory_model(
        result=result,
        problem=problem,
        timepoints=np.linspace(0, 130, 1301),
        state_names=[],
        state_ids=[],
        observable_ids=obs_names,
    )
    label_colors = [
        (0.502, 0.733, 0.509),
        (0.944, 0.730, 0.418),
        (0.944, 0.558, 0.744),
        (0.733, 0.502, 0.509),
        (0.509, 0.502, 0.733),
    ]
    color_t0 = (0.647, 0.745, 0.858)

    for i_line, line in enumerate(ax.lines):
        if i_line == 0:  # skip first line as it is the original
            continue
        xx = line.get_xdata()
        yy = line.get_ydata()
        xx = np.concatenate(([-40], xx))
        yy = np.concatenate(([yy[0]], yy))
        line.set_data(xx, yy)
        line.set_markevery([1201])
        line.set_marker("o")
        line.set_c(label_colors[i_line])
        ax.annotate(
            # Label and coordinate
            "",
            xy=(
                line.get_xdata()[1201] - i_line * (120 / (n_labels - 1)),
                line.get_ydata()[1201],
            ),
            xytext=(line.get_xdata()[1201], line.get_ydata()[1201]),
            color=line.get_color(),
            # Custom arrow
            arrowprops=dict(
                arrowstyle="->",
                color=line.get_color(),
                lw=1, ls="--"
            ),
        )
    for i in range(n_labels):  # replicate points on the original line
        ax.plot(
            [ax.lines[i].get_xdata()[1201] - (i) * (120 / (n_labels - 1))],
            [ax.lines[i].get_ydata()[1201]],
            color=label_colors[i],
            marker="o",
        )

    xticks = [(i - 1) * (120 / (n_labels - 1)) for i in range(n_labels + 1)]
    xticks[0] = -30
    xlabels = [f"$T_{i-1}" for i in range(n_labels)]
    ax.set_xticks(xticks, labels=xlabels)
    ax.axvline(x=120, color="#A7A7A7", linestyle="dashed", linewidth=1)  #
    # measurement
    ax.axvline(x=-30, color=color_t0, linestyle="dashed", linewidth=1)  # T0
    ax.get_xticklabels()[0].set_color(color_t0)
    for i_timepoint, color in enumerate(label_colors):
        ax.axvline(
            x=i_timepoint * (120 / (n_labels - 1)),
            color=color,
            linestyle="dashed",
            linewidth=1,
        )
        ax.get_xticklabels()[i_timepoint + 1].set_color(color)

    labeling = [item.get_text() for item in ax.get_yticklabels()]
    labeling[1] = "0"
    ax.set_yticklabels(labeling)
    ax.set_xlim((-40, 125))
    plt.title("")
    plt.ylabel("")
    plt.xlabel("")
    ax.figure.set_size_inches(5, 2.5)
    ax.set_yticks([])

    return ax
