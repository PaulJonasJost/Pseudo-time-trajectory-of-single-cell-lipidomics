"""
This module defines the constants and variables used in the SBML model.

It includes definitions for:
- LABELS: The labels used in the model
- SPECIES: The species in the model and their number of fatty acids
- PARAMETERS: The parameters of the model and their values
- INITIAL_VALUES: The initial values of the species
- REACTIONS: The reactions in the model, loaded from a TSV file
"""
import pandas as pd

from sbml_model_util import Reaction

DIFF_VALUES = True

# LABELS
LABELS = ["ul", "16_1", "16_2", "17_1"]
# SPECIES - format {name: #fas}
SPECIES = {
    "FA": 1,
    "GLY": 0,
    "G3P": 0,
    "LPA": 1,
    "PA": 2,
    "DAG": 2,
    "TAG": 3,
    "MAG": 1,
    "Choline": 0,
    "ph_C": 0,
    "CDP_C": 0,
    "Ethanolamine": 0,
    "ph_EA": 0,
    "CDP_EA": 0,
    "Serine": 0,
    "PC": 2,
    "PE": 2,
    "PS": 2,
    "LPC": 1,
    "G3PC": 0,
    "LPE": 1,
    "G3PE": 0,
}
# PARAMETERS - format: {Name: value}
PARAMETERS = {
    "k_glycerol_to_g3p": 0.01,
    "k_g3p_to_lpa": 0.01,
    "k_lpa_to_pa": 0.01,
    "k_pa_to_dag": 0.01,
    "k_dag_to_tag": 0.01,
    "k_choline_ph": 0.01,
    "k_choline_cdp": 0.01,
    "k_ea_ph": 0.01,
    "k_ea_cdp": 0.01,
    "k_pe_syn": 0.01,
    "k_pc_syn": 0.01,
    "k_pc_to_ps": 0.01,
    "k_pe_to_ps": 0.01,
    "k_ps_to_pe": 0.01,
    "k_pe_to_pc": 0.01,
    "k_pc_to_lpc": 0.01,
    "k_lpc_to_g3pc": 0.01,
    "k_g3pc_to_g3p": 0.01,
    "k_pe_to_lpe": 0.01,
    "k_lpe_to_g3pe": 0.01,
    "k_g3pe_to_g3p": 0.01,
    "k_dag_to_mag": 0.01,
    "k_mag_to_dag": 0.01,
    "k_tag_to_dag": 0.01,
    "k_mag_to_gly": 0.01,
    "k_serine_syn": 0.01,
    "k_gly_syn": 0.01,
    "k_fa_influx": 0.01,
    "k_serine_deg": 0.01,
    "k_pc_trans": 0.01,
    "k_fa_deg": 0.01,
    "k_gly_deg": 0.01,
    "k_choline_deph": 0.01,
    "k_ea_deph": 0.01,
    "k_gly_deph": 0.01,
    "k_choline_syn": 0.01,
    "k_choline_deg": 0.01,
    "k_choline_cdp_back": 0.01,
    "k_ea_cdp_back": 0.01,
    "k_ea_deg": 0.01,
    # preequilibration parameter
    "preeq_switch": 1,
}
if DIFF_VALUES:
    for i_par, par in enumerate(PARAMETERS):
        if par == "preeq_switch":
            continue
        PARAMETERS[par] += i_par * 0.01


for label in LABELS:
    PARAMETERS[f"k_fa_{label}_influx"] = 0.01

# INITIAL VALUES - format: {Name: value}
INITIAL_VALUES = {
    "FA_ul_initial": 1,
    "GLY_ul_initial": 1,
    "G3P_ul_initial": 1,
    "LPA_ul_initial": 1,
    "PA_ul_initial": 1,
    "DAG_ul_initial": 1,
    "TAG_ul_initial": 1,
    "Choline_ul_initial": 1,
    "ph_C_ul_initial": 1,
    "CDP_C_ul_initial": 1,
    "Ethanolamine_ul_initial": 1,
    "ph_EA_ul_initial": 1,
    "CDP_EA_ul_initial": 1,
    "Serine_ul_initial": 1,
    "PC_ul_initial": 1,
    "PE_ul_initial": 1,
    "PS_ul_initial": 1,
    "LPC_ul_initial": 1,
    "G3PC_ul_initial": 1,
    "LPE_ul_initial": 1,
    "G3PE_ul_initial": 1,
    "MAG_ul_initial": 1,
}

# load reaction from table
df_reactions = pd.read_csv(
    "reaction_table.tsv", sep="\t", na_filter=False
)
REACTIONS = [
    Reaction(**reaction) for reaction in df_reactions.to_dict("records")
]

T_PULSE = 30  # time of each pulse, it is for now assumed that this is constant
T_BETWEEN = None  # time in between pulses
