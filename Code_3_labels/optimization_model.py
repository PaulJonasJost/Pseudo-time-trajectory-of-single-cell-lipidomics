import datetime
import logging
import sys
import os

import petab
import pypesto.petab
import pypesto.optimize as optimize


if __name__ == '__main__':
    model_name = 'lipidomics_2023_08_29'
    yaml_file = f'{os.path.dirname(__file__)}/../Petab_models_230829' \
                f'/{model_name}_{sys.argv[1]}/' \
                f'{model_name}_{sys.argv[1]}.yaml'
    petab_problem = petab.Problem.from_yaml(yaml_file)
    importer = pypesto.petab.PetabImporter(
        petab_problem=petab_problem,
        output_folder=f'../Code/amici_models/{model_name}',
        model_name=model_name
    )
    problem = importer.create_problem(n_threads=8,
                                      guess_steadystate=True,
                                      # force_compile=True,
                                      )

    problem.objective.amici_solver.setRelativeTolerance(1e-12)  # 1e-8
    problem.objective.amici_solver.setAbsoluteTolerance(1e-15)  # 1e-12

    optimizer_options = {
        'maxiter': 1e4,
        'fatol': 1e-12,
        'frtol': 1e-12
    }


    optimizer = optimize.FidesOptimizer(
        options=optimizer_options,
        verbose=logging.WARN
    )
    engine = pypesto.engine.MultiProcessEngine()
    # engine = pypesto.engine.SingleCoreEngine()
    n_starts = 1000

    # creation of savefilename
    time = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M")
    save_file = f"../results_h5/{model_name}_" \
                f"{sys.argv[1]}_{n_starts}starts.hdf5"

    history_options = pypesto.HistoryOptions(
        trace_record=True,
        storage_file=save_file
    )


    ref = problem.objective(
        petab_problem.x_nominal_free_scaled,
        return_dict=True
    )
    print(ref)

    result = optimize.minimize(
        problem=problem,
        optimizer=optimizer,
        n_starts=n_starts,
        engine=engine,
        # filename=None,
        filename=save_file,
        # history_options=history_options
    )

    print(f'Done, the results summary is: \n '
          f'{result.optimize_result.summary()}')
    print(f"The value of the petab problem parameters (true ones) is: "
          f"{ref['fval']}")
