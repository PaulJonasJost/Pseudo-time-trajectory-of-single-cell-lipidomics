# Readme for the repository

This is a repository for the corresponding Paper 
"Pseudo-time trajectory of single cell lipidomics".

The repository contains the following structure:

* Code_3_labels: Code for the analysis of the 3 labels dataset
* Code_multilabels: Code for the analysis of the multilabels dataset
* Toy_models: Code for the analysis of an examplatory model
* results_h5: Optimization results saved in HDF5 format

Here is a short description of the essential code in the repository, per folder:

* Code_3_labels: 
    * analyze_results.py: Helper functions to pool the optimization results.
    * create_parameters.py: Function to create the parameter sets.
    * sbml_model.py: Variables such as Parameter Values and Reactions for the base model.
    * sbml_model_utils.py: Helper functions and rules to create labeled models.
    * sbml_model_create.py: Running this file creates the models.
    * visualizations.py: Creates **Figure 2 F**
    * parameter_set.csv: Parameter sets for the models.
    * reaction_table.tsv: Reactions of the models.
    * optimization_model.py: Run the optimization for the models.
    * Figure4.ipynb: Creates parts of **Figures 3 and 4**
    * Figure5.ipynb: Creates parts of **Figure 5**
* Code_multilabels
  * The files here have the same function as the files in Code_3_labels, but for the multilabels dataset.
* Toy_models
  * Figure1.ipynb: Creates parts of **Figure 1**