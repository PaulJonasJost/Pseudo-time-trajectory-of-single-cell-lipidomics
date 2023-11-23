"""
A file to run visualizations with an existing result h5 file.
"""
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib
import numpy as np
import os
import pypesto.visualize as visualize
import pypesto.petab
import petab
import pypesto.store
import amici.plotting
import pypesto.visualize.model_fit as model_fit
import petab.visualize as petab_vis

from typing import List, Sequence, Optional


def plot_time_trajectories(results_dir, res_filename, result, problem):
    obs_names = [
        # 'observable_DAG_ul_ul',
        'observable_DAG_16_1_16_1',
        'observable_DAG_16_2_16_2',
        'observable_DAG_17_1_17_1'
    ]
    ax = model_fit.time_trajectory_model(
        result=result,
        problem=problem,
        timepoints=np.linspace(0, 110, 1101),
        state_names=[],
        state_ids=[],
        observable_ids=obs_names)
    label_colors = [
        (0.502, 0.733, 0.509),
        (0.944, 0.730, 0.418),
        (0.944, 0.558, 0.744)
    ]
    color_t0 = (0.647, 0.745, 0.858)

    # working on the legend
    handles, labels = ax.get_legend_handles_labels()
    for i_handle, handle in enumerate(handles):
        labels[i_handle] = f'$DAG_{i_handle+1}$'
        handle.set_color(label_colors[i_handle])


    for i_line, line in enumerate(ax.lines):
        if i_line == 0:  # skip first line as it is the original
            continue
        xx = line.get_xdata()
        yy = line.get_ydata()
        xx = np.concatenate(([-40], xx))
        yy = np.concatenate(([yy[0]], yy))
        line.set_data(xx, yy)
        line.set_markevery([901])
        line.set_marker('o')
        line.set_c(label_colors[i_line])
        ax.annotate(
            # Label and coordinate
            f'',
            xy=(line.get_xdata()[901]-i_line*30, line.get_ydata()[901]),
            xytext=(line.get_xdata()[901], line.get_ydata()[901]),
            color=line.get_color(),
            # Custom arrow
            arrowprops=dict(arrowstyle='->', color=line.get_color(),
                            lw=1, ls='--')
        )
        ax.annotate(
            f'$T_{i_line+1} - T_1$',
            xy=((line.get_xdata()[901] - 30 + line.get_xdata()[901])/2,
                line.get_ydata()[901]),
            xytext=(0, 10), textcoords='offset pixels',
            color=line.get_color(),
            ha='center'
        )
    for i in range(3): # replicate points on the original line
        ax.plot(
            [ax.lines[i].get_xdata()[901] - (i) * 30],
            [ax.lines[i].get_ydata()[901]],
            color=label_colors[i],
            marker='o'
        )

    xticks = [-30, 0, 30, 60, 90]
    xlabels = ['$T_0$', '$T_1$', '$T_2$', '$T_3$', 'Measurement']
    ax.set_xticks(xticks, labels=xlabels)
    ax.axvline(x=90, color='#A7A7A7', linestyle='dashed', linewidth=1)  # measurement
    ax.axvline(x=-30, color=color_t0, linestyle='dashed', linewidth=1)  # T0
    ax.get_xticklabels()[0].set_color(color_t0)
    for i_timepoint, color in enumerate(label_colors):
        ax.axvline(
            x=i_timepoint * 30,
            color=color,
            linestyle='dashed',
            linewidth=1,
        )
        ax.get_xticklabels()[i_timepoint + 1].set_color(color)
        # add item to legend
        handles.append(
            mlines.Line2D(
                [], [],
                color=color,
                marker=f'$T_{i_timepoint + 1}$',
                linestyle='None',
                markersize=10,
                markeredgewidth=0.1,
            )
        )
        labels.append(f'Start Pulse {i_timepoint + 1}')


    ax.get_xticklabels()[-1].set_color('#A7A7A7')
    for t in ax.xaxis.get_ticklines(): t.set_color('green')

    labeling = [item.get_text() for item in ax.get_yticklabels()]
    labeling[1] = '0'
    ax.set_yticklabels(labeling)
    ax.set_xlim((-40, 95))
    plt.title('')
    plt.ylabel('')
    plt.xlabel('')
    ax.legend(handles=handles, labels=labels)

    ax.figure.set_size_inches(5, 2.5)
    ax.set_yticks([])

    plt.savefig(
        results_dir + "/" + res_filename+f'_timeTrajectory_DAG.svg'
    )


MODEL_FILE = "lipidomics_2023_08_30_1"
N_OPT = 1000

if __name__ == '__main__':
    petab_problem = petab.Problem.from_yaml(
        f'../Petab_models_230829/3_labels/{MODEL_FILE}'
        f'/{MODEL_FILE}.yaml')

    importer = pypesto.petab.PetabImporter(
        petab_problem=petab_problem,
        output_folder=f'./amici_models/{MODEL_FILE}',
        model_name=MODEL_FILE
    )
    problem = importer.create_problem()
    problem.objective.amici_solver.setRelativeTolerance(1e-8)
    problem.objective.amici_solver.setAbsoluteTolerance(1e-12)

    res_filename = f'{MODEL_FILE}_{N_OPT}starts'  # _{OPTIMIZATION_DAY}'
    result_file = '../results_h5/'+res_filename+'.hdf5'

    result = pypesto.store.read_result(result_file, optimize=True,
                                       problem=False)
    result.problem = problem

    # use the folder this file is in as results directory
    results_dir = os.path.dirname(os.path.abspath(__file__))



    #
    # Figure 2 F schematic
    plot_time_trajectories(
        results_dir=results_dir,
        result=result,
        res_filename=res_filename,
        problem=problem
    )

    print('Done')
