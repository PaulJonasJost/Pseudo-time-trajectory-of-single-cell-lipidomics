# Readme for the repository

This is a repository for the corresponding paper 
"Pseudo-time trajectory of single cell lipidomics".

To reproduce the workflows and figures of the paper, please follow the instructions below.

## Requirements

The code provided here was originally run with python 3.11.1, though -- just as for the 
software packages --  we estimate that newer version should also work. The code should be 
compatible with Windows, MacOS and Linux. They Key software packages used are:
* [pyPESTO](https://github.com/ICB-DCM/pyPESTO) v0.3.1 (pip installable)
* [AMICI](https://github.com/AMICI-dev/AMICI) v0.18.0 (pip installable, but requires
   a C++ compiler, see 
   [installations instructions](https://amici.readthedocs.io/en/v0.18.0/python_installation.html))
* [Fides](https://github.com/fides-dev/fides) v0.7.8 (pip installable)
* [simplesbml](https://github.com/sys-bio/simplesbml) v2.3.0 (pip installable)

All other packages should be automatically installed with the above packages.

## Repository structure

The repository contains the following structure:

* Code_3_labels: Code for the analysis of the 3 labels dataset
* Code_multilabels: Code for the analysis of the multilabels dataset
* Results_h5: Results of the optimization, which can be used to create the figures 
  without running the optimization again. The structure of a single `.hdf5` File can 
  be found [here](https://pypesto.readthedocs.io/en/latest/storage_format.html).

We also provide the base model (See **Figure 2A**) and an example model with 3 labels 
as SBML files.

## Running Code

We split the Code into two workflows:
1. 3 Label Study Workflow: This workflow creates the models and data and runs the
   optimization of the main dataset of the paper. Visualizations that recreate parts 
   of the figures are included in a separate Jupyter notebook.
2. Multilabel Study Workflow: This workflow creates the models and data and runs the
   optimization of the multilabels dataset. Visualizations that recreate parts of the 
   figures are included in a separate Jupyter notebook.

Each workflow can be run independently and the folders are set up such that only that 
specific folder is needed to run the code (even though this creates some redundancy). 
Additionally, as the optimization takes a substantial amount of time, the creation of 
data, models and optimization can be skipped. The original results of the 
optimization are provided in the `Results_h5` folder, which can be used to create the 
figures without running the optimization again. A separate folder "Figure_1" contains
a very small model, that has been used to create **Figure 1** of the paper.

### 3 Label Study Workflow

1. Run `create_parameters.py` to create the parameter sets. You can skip this step 
   if you want to use the precomputed parameter sets in the `Code_3_labels` folder.
```shell
python create_parameters.py
```
2. Run `sbml_model_create.py` to create the sbml and petab models. They will be saved 
   the folder `Petab_models_230829/3_labels`.
```shell
python sbml_model_create.py
```
3. Run `optimization_model.py` to run the optimization for the models. The results will 
   be saved in the folder `Results_h5`. You either specifiy an index of the model to 
   run (counting from zero, corresponding to the parameter set) in the command line, 
   or run the command with "ALL" instead to run all models. This will be computationally 
   intensive!!
```shell
python optimization_model.py $INDEX
```
or
```shell
python optimization_model.py ALL
```
4. Visualizations
The notebooks are set up to use the results from the 
`Results_h5` folder, so you don't need to run the optimization again. Though the 
`sbml_model_create.py` should be run beforehand. Specifically, the visualizations are:
* `Visualizations_3_labels.ipynb`: Creates **Figure 4B, Figure 5B,C** of the paper. 
  This visualization uses the results from the `Results_h5` folder.
* `visualizations.py`: Creates the basis for **Figure 2A**. Also uses the results from 
  the `Results_h5` folder.

For a more detailed description of what each files does see [Code description](#code-description).

### Multilabel Study Workflow

1. Run `create_parameters.py` to create the parameter sets. You can skip this step 
   if you want to use the precomputed parameter sets in the `Code_multilabels` folder.
```shell
python create_parameters.py
```
2. Run `sbml_model_create.py` to create the sbml and petab models. They will be saved 
   the folder `Petab_models_230829/multilabels`.
```shell
python sbml_model_create.py
```
3. Run `optimization_model.py` to run the optimization for the models. The results will 
   be saved in the folder `Results_h5`. While you can specify the index of the model, 
   the index gets converted to the corresponding parameter set and label. In practice 
   this means that supplying an index $ind$, we get the parameter index $i_{param} = ind 
   \mod 200$ and the number of labels $n_{labels} = ind // 200$. It might be preferred 
   to run the command with "ALL" instead to run all models. This will be computationally 
   intensive!!
```shell
python optimization_model.py $INDEX
```
or
```shell
python optimization_model.py ALL
```
4. Visualizations
The notebooks are set up to use the results from the 
`Results_h5` folder, so you don't need to run the optimization again. Though the 
`sbml_model_create.py` should be run beforehand. Specifically, the visualizations are:
* `Visualizations_multilabels.ipynb`: Creates **Figure 6** of the paper. 
  This visualization uses the results from the `Results_h5` folder.

### Visualization of Figure 1

The folder `Folder1_creation` serves to create **Figure 1B**. The notebook `Figure1.
ipynb` can be run without any dependencies on other folders. This visualization does 
not use any model or result from before, but rather a **Toy Model**.

## Code description

Here is a short description of the essential code in the repository:

* Code_3_labels and Code_multilabels: 
    * create_parameters.py: Function to create the parameter sets.
    * sbml_model.py: Variables such as Parameter Values and Reactions for the base model.
    * sbml_model_utils.py: Helper functions and rules to create labeled models.
    * sbml_model_create.py: Running this file creates the models.
    * parameter_set.csv: Parameter sets for the models.
    * reaction_table.tsv: Reactions of the models.
    * optimization_model.py: Run the optimization for the models.
* Visualizations:
    * analyze_results.py: Defines the CumulativeResult class used in the visualizations.
    * visualizations.py: Creates the basis for **Figure 2A**.
    * Visualizations_3_labels.ipynb: Creates **Figure 4B, Figure 5B,C** of the paper.
    * Visualizations_multilabels.ipynb: Creates **Figure 6** of the paper.
    * Figure_1_creation:
        * Figure1.ipynb: Creates **Figure 1B**.
* Results_h5:
    * Results of the optimization, which can be used to create the figures without running 
      the optimization again. For an overview of the structure of each `.hdf5` file, 
      we refer to the
      [pyPESTO documentation](https://pypesto.readthedocs.io/en/latest/storage_format.html).