import logging
import simplesbml
import itertools
import numpy as np

from typing import Iterable, List
from collections import Counter
from copy import deepcopy

from itertools import combinations_with_replacement

logger = logging.getLogger(__name__)


def uniqueCombinations(lstList, n_labels):
    """
    Get unqiue list of combinations from the original list.
    :param lstList: original list
    :param n_labels: number of labels in new combination
    :return: list of unique lits
    """
    all_combos = itertools.combinations(lstList, n_labels)
    unique_label_combination = set()
    for combo_candidate in all_combos:
        if combo_candidate in unique_label_combination:
            continue
        unique_label_combination.add(combo_candidate)
    return [list(lab) for lab in unique_label_combination]


def compare_lists(list1, list2):
    """
    Compare whether two lists are identical up to permutation level.
    """
    return Counter(list1) == Counter(list2)


def combine_species_label(base_name: str, label: List[str]):
    """
    Combine the base name of the species with its labels.

    :param base_name: str, base name of the species
    :param label: list of labels to add to the species.
    :return: str, name of the labeled species
    """
    if len(label) == 0:
        return base_name
    elif len(label) == 1:
        label_complete = "".join(label)
    else:
        label_complete = "_".join(label)
    return f"{base_name}_{label_complete}"


def create_labeled_species(
    species_name: str,
    labels: Iterable[str],
    model: simplesbml.SbmlModel,
    n_labels: int = 1,
):
    """
    Create species with different label combinations.

    :param species_name: Base name of the species.
    :param model: The sbml Model which to change.
    :param n_labels: Number of labels the species can have (e.g. DG=2, TG=3)
    :param labels: Labels that are used.

    :return: List of names of the labeled species.
    """
    # create combinations of labels
    labels_per = [
        list(lab)
        for lab in combinations_with_replacement(labels, n_labels)
    ]
    species = []
    for lab in labels_per:
        labeled_species = combine_species_label(species_name, lab)
        model.addSpecies(labeled_species, 0)
        species.append(labeled_species)
    return species


def create_initial_Assignment(
    species_name: str,
    model: simplesbml.SbmlModel,
    n_labels: int = 1,
):
    """
    Create an initial assignment for a species with unlabeled fatty acids.

    This function creates an initial assignment for a species with unlabeled fatty acids.
    It constructs the species name with the appropriate number of "ul" (unlabeled) labels
    and assigns the initial value parameter to it.

    Parameters
    ----------
    species_name : str
        The base name of the species.
    model : simplesbml.SbmlModel
        The SBML model to add the initial assignment to.
    n_labels : int, optional
        The number of labels the species can have. Default is 1.

    Returns
    -------
    None
        The function modifies the model in-place but does not return any value.
    """
    initial_par_name = f"{species_name}_ul_initial"
    species_name_ul = species_name
    for i in range(n_labels):
        species_name_ul += "_ul"
    model.addInitialAssignment(symbol=species_name_ul, math=initial_par_name)


def reactions_add_label(
    react_w_label: str,
    prod_w_label: str,
    n_label_react: int,
    n_label_prod: int,
    k_reaction: str,
    labels: Iterable[str],
    model: simplesbml.SbmlModel,
    byproduct: Iterable[str] = None,
    reactant_2: str = None,
):
    """
    Create a list of reactions based on the labels of the reactants and
    products. Used if a label is added (e.g. dag -> tag). If all labels are
    kept, use reactions_keep_label instead. If labels are removed,
    use reactions_rem_label.

    :param react_w_label: Labeled reactant, should provide only the base
    name, e.g. 'tg' for 'tg_l1_l2_l3'.
    :param prod_w_label: Labeled products, should provide only the base
    name, e.g. 'tg' for 'tg_l1_l2_l3'.
    :param n_label_react: number of labels on reactant
    :param n_label_prod: number of labels on product
    :param k_reaction: reaction rate parameter
    :param model: The SbmlModel of concern.
    :param byproduct: Byproduct without a label
    :param labels: Possible labels.
    :param reactant_2: Addition reactants without a label.
    """
    # create combinations of labels
    labels_react = [
        list(lab)
        for lab in combinations_with_replacement(labels, n_label_react)
    ]
    labels_prod = [
        list(lab)
        for lab in combinations_with_replacement(labels, n_label_prod)
    ]
    i_react = 0
    for label_curr in labels_react:
        for label_add in labels:
            label = deepcopy(label_curr)
            label.append(label_add)
            # convert to label ordering as it is in the products
            label = next(
                (lab for lab in labels_prod if compare_lists(label, lab)), None
            )
            # combine labeled fatty acid with labeled species to reactants
            reactants = [combine_species_label(react_w_label, label_curr)]
            products = [combine_species_label(prod_w_label, label)]

            if reactant_2 is None:
                reactants += [combine_species_label("FA", [label_add])]
            else:
                reactants += [combine_species_label(reactant_2, [label_add])]

            if byproduct is not None:
                products += byproduct
            expression = f'{k_reaction} * {" * ".join(reactants)}'
            model.addReaction(
                reactants=reactants,
                products=products,
                expression=expression,
                rxn_id=f"{react_w_label}_to_{prod_w_label}_{i_react}",
            )
            i_react += 1


def reactions_rem_label(
    react_w_label: str,
    prod_w_label: str,
    n_label_react: int,
    n_label_prod: int,
    k_reaction: str,
    labels: Iterable[str],
    model: simplesbml.SbmlModel,
    byproduct: str = None,
    reactant_2: str = None,
):
    """
    Create a list of reactions based on the labels of the reactants and
    products. Used if a label is removed (e.g. tag -> dag). If all labels are
    kept, use reactions_keep_label instead. If a label is added,
    use reactions_add_label.

    :param react_w_label: Labeled reactants, should provide only the base
    name, e.g. 'tg' for 'tg_l1_l2_l3'.
    :param prod_w_label: Labeled products, should provide only the base
    name, e.g. 'tg' for 'tg_l1_l2_l3'.
    :param n_label_react: number of labels on reactants
    :param n_label_prod: number of labels on product
    :param k_reaction: reaction rate parameter
    :param model: The SbmlModel of concern.
    :param byproduct: Byprodukt of reaction that receives the label.
    Assumed to be FA if not set otherwise. Is not allowed to have a label
    already.
    :param labels: Possible labels.
    :param reactant_2: Addition reactants without a label.
    """
    labels_react = [
        list(lab)
        for lab in combinations_with_replacement(labels, n_label_react)
    ]
    labels_prod = [
        list(lab)
        for lab in combinations_with_replacement(labels, n_label_prod)
    ]
    i_react = 0
    for label_react in labels_react:
        label_counter = Counter(label_react)
        for key in label_counter:
            label = deepcopy(label_counter)
            # calculate propensities for losing specific label
            total = np.sum([label[key] for key in label])
            propensity = label[key] / total
            # remove specific label once
            label[key] -= 1
            # convert counter to list
            label = [x for x in label.elements()]
            # convert to label ordering as it is in the products
            label = next(
                (lab for lab in labels_prod if compare_lists(label, lab)), None
            )
            # combine labeled fatty acid with labeled species to reactants
            reactants = [combine_species_label(react_w_label, label_react)]
            products = [combine_species_label(prod_w_label, label)]
            expression = f"{k_reaction} * propensity * {reactants[0]}"
            if byproduct is not None:
                products += combine_species_label(byproduct, [key])
            else:  # default byproduct is fatty acid
                products += [combine_species_label("FA", [key])]
            if reactant_2 is not None:
                reactants += [reactant_2]
                expression += f" * {reactant_2}"
            model.addReaction(
                reactants=reactants,
                products=products,
                expression=expression,
                rxn_id=f"{react_w_label}_to_{prod_w_label}_{i_react}",
                local_params={"propensity": propensity},
            )
            i_react += 1


def reactions_keep_label(
    react_w_label: str,
    prod_w_label: str,
    n_labels: int,
    k_reaction: str,
    labels: Iterable[str],
    model: simplesbml.SbmlModel,
    byproduct: str = None,
    reactant_2: str = None,
):
    """
    Create a list of reactions based on the labels of the reactants and
    products. Used if the labels are transferred (e.g. dag -> PE).
    If labels are removed, use reactions_rem_label instead.
    If a label is added, use reactions_add_label.

    :param react_w_label: Labeled reactant, should provide only the base
    name, e.g. 'tg' for 'tg_l1_l2_l3'.
    :param prod_w_label: Labeled product, should provide only the base
    name, e.g. 'tg' for 'tg_l1_l2_l3'.
    :param n_labels: number of labels. Has to be the same on reactant and
    product.
    :param k_reaction:  reaction rate parameter
    :param model: The SbmlModel of concern.
    :param byproduct: Byproduct without a label.
    :param reactant_2: Addition reactant without a label.
    :param labels: Possible labels.
    """
    labels_react = [
        list(lab) for lab in combinations_with_replacement(labels, n_labels)
    ]
    i_react = 0
    for label_react in labels_react:
        # combine label with labeled reactant and labeled product
        reactants = [combine_species_label(react_w_label, label_react)]
        products = [combine_species_label(prod_w_label, label_react)]
        expression = f"{k_reaction} * {reactants[0]}"
        if reactant_2 is not None:
            reactants += [reactant_2]
            expression += f" * {reactant_2}"
        if byproduct is not None:
            products += [byproduct]
        model.addReaction(
            reactants=reactants,
            products=products,
            expression=expression,
            rxn_id=f"{react_w_label}_to_{prod_w_label}_{i_react}",
        )
        i_react += 1


def reactions_synthesis(
    prod_w_label: str,
    n_label_prod: int,
    k_reaction: str,
    labels: Iterable[str],
    model: simplesbml.SbmlModel,
):
    """
    Create a list of reaction, synthesizing the different labeled products.

    :param prod_w_label: Base name of the product.
    :param n_label_prod: Number of labels.
    :param k_reaction: Reaction parameter.
    :param labels: The possible labels.
    :param model: The sbml model.
    """
    expression = f"{k_reaction}"
    if n_label_prod == 0:
        model.addReaction(
            reactants=[],
            products=[prod_w_label],
            expression=expression,
            rxn_id=f"{prod_w_label}_synthesis",
        )
        return
    labels_prod = [
        list(lab)
        for lab in combinations_with_replacement(labels, n_label_prod)
    ]
    i_react = 0
    for label_prod in labels_prod:
        products = [combine_species_label(prod_w_label, label_prod)]
        model.addReaction(
            reactants=[],
            products=products,
            expression=expression,
            rxn_id=f"{prod_w_label}_synthesis_{i_react}",
        )
        i_react += 1


def reactions_degradation(
    react_w_label: str,
    n_label_react: int,
    k_reaction: str,
    labels: Iterable[str],
    model: simplesbml.SbmlModel,
    product: str = None,
    keep_fas: bool = False,
):
    """Creates a degradation reaction for the different labeled species.

    :param react_w_label: Base name of the reactant.
    :param n_label_prod: Number of labels.
    :param k_reaction: Reaction parameter.
    :param labels: The possible labels.
    :param model: The sbml model.
    :param product: The prudoct from the reaction.
    :param keep_fas: Whether to keep the FAs lost in the reaction in the system
    """
    labels_react = [
        list(lab)
        for lab in combinations_with_replacement(labels, n_label_react)
    ]
    i_react = 0
    for label_react in labels_react:
        products = []
        if product is not None:
            products = [product]
        reactants = [combine_species_label(react_w_label, label_react)]
        expression = f"{k_reaction} * {reactants[0]}"
        if keep_fas:
            label_counter = Counter(label_react)
            products.extend(
                [f"{label_counter[key]} FA_{key}" for key in label_counter]
            )
        model.addReaction(
            reactants=reactants,
            products=products,
            expression=expression,
            rxn_id=f"{react_w_label}_degradation_{i_react}",
        )
        i_react += 1


class Reaction(dict):
    def __init__(
        self,
        react_w_label: str,
        prod_w_label: str,
        n_label_react: int,
        n_label_prod: int,
        k_reaction: str,
        name: str = None,
        byproduct: str = None,
        reactant_2: str = None,
    ):
        """Constructor for Reaction."""
        super().__init__()
        self.name = None if name == "None" else name
        self.k_reaction = k_reaction
        self.react_w_label = None if react_w_label == "None" else react_w_label
        self.prod_w_label = None if prod_w_label == "None" else prod_w_label
        self.n_label_react = n_label_react
        self.n_label_prod = n_label_prod
        self.byproduct = None if byproduct == "None" else byproduct
        self.reactant_2 = None if reactant_2 == "None" else reactant_2

    def create_reaction(
            self, labels: Iterable[str], model: simplesbml.SbmlModel
    ):
        """
        Create the reaction.

        :param labels: Possible labels.
        :param model: The SbmlModel of concern.
        """
        # decide which kind of reaction it is
        if self.react_w_label is None:
            reactions_synthesis(
                prod_w_label=self.prod_w_label,
                n_label_prod=self.n_label_prod,
                k_reaction=self.k_reaction,
                labels=labels,
                model=model,
            )
            return
        if self.prod_w_label is None:
            reactions_degradation(
                react_w_label=self.react_w_label,
                n_label_react=self.n_label_react,
                k_reaction=self.k_reaction,
                labels=labels,
                model=model,
            )
            return
        label_diff = self.n_label_prod - self.n_label_react
        if label_diff == 0:
            reactions_keep_label(
                react_w_label=self.react_w_label,
                prod_w_label=self.prod_w_label,
                n_labels=self.n_label_prod,
                k_reaction=self.k_reaction,
                labels=labels,
                model=model,
                byproduct=self.byproduct,
                reactant_2=self.reactant_2,
            )
        elif label_diff == 1:
            reactions_add_label(
                react_w_label=self.react_w_label,
                prod_w_label=self.prod_w_label,
                n_label_react=self.n_label_react,
                n_label_prod=self.n_label_prod,
                k_reaction=self.k_reaction,
                labels=labels,
                model=model,
                byproduct=self.byproduct,
                reactant_2=self.reactant_2,
            )
        elif label_diff == -1:
            reactions_rem_label(
                react_w_label=self.react_w_label,
                prod_w_label=self.prod_w_label,
                n_label_react=self.n_label_react,
                n_label_prod=self.n_label_prod,
                k_reaction=self.k_reaction,
                labels=labels,
                model=model,
                byproduct=self.byproduct,
                reactant_2=self.reactant_2,
            )
        else:
            if self.n_label_prod != 0:
                raise NotImplementedError(
                    "There is not yet a function to "
                    "transfer multiple labels in one "
                    "reaction."
                )
            keep_fas = True
            if self.prod_w_label is None:
                keep_fas = False
            reactions_degradation(
                react_w_label=self.react_w_label,
                n_label_react=self.n_label_react,
                k_reaction=self.k_reaction,
                labels=labels,
                model=model,
                product=self.prod_w_label,
                keep_fas=keep_fas,
            )


def pulse_chase_w_labels(
    labels: Iterable[str],
    model: simplesbml.SbmlModel,
    t_pulse: float = 30,
    t_total: float = None,
):
    # if t_total is not None, the t_pulse will be t_total / len(labels) - 1
    if t_total is not None:
        # one labels = unlabeled
        t_pulse = 0
        if len(labels) > 1:
            t_pulse = t_total / (len(labels) - 1)
    labels_copy = deepcopy(labels)
    # create unlabeled fa influx, only before steady state or after all pulses
    model.addAssignmentRule(
        var="k_fa_ul_influx",
        math=f"(1-preeq_switch) * k_fa_influx + "
        f"piecewise(k_fa_influx, Time >="
        f" {(len(labels_copy)-1) * t_pulse}, 0)",
    )
    model.addReaction(
        reactants=[],
        products=["FA_ul"],
        expression="k_fa_ul_influx",
        rxn_id="FA_ul_synthesis",
    )
    labels_copy.remove("ul")
    for i_label, label in enumerate(labels_copy):
        k = f"k_fa_{label}_influx"
        t_start = i_label * t_pulse
        t_end = (i_label + 1) * t_pulse
        piecewise_function = (
            f"preeq_switch * piecewise(k_fa_influx,"
            f" {t_start} <= Time < {t_end}, 0)"
        )
        model.addAssignmentRule(var=k, math=piecewise_function)
        model.addReaction(
            reactants=[],
            products=[f"FA_{label}"],
            expression=f"{k}",
            rxn_id=f"FA_{label}_synthesis",
        )
