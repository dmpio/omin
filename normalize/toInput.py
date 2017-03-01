# -*- coding: utf-8 -*-

import re
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from ..utils import StringTools
from ..utils import SelectionTools

from omin.normalize.methods import *


def normFactors(peptide_data):
    """Takes peptide abundance data and returns normalization factors.

    Normalization factors are derived by the taking the sum of each column in
    DataFrame then dividing each sum by the mean of all the sums.

    Parameters
    ----------
    peptide_data : DataFrame

    Returns
    -------
    norm_factors: DataFrame

    """
    norm_factors = peptide_data.sum() / peptide_data.sum().mean()
    return norm_factors


def normalizeTo(different, normal):
    """Normalizes one DataFrame to another.

    The 'different' DataFrame is normalized to the 'normal' DataFrame using the
    normFactors function.

    Parameters
    ----------
    different : DataFrame
    normal : DataFrame

    Returns
    -------
    normalized : DataFrame
    """
    normalized = different / normFactors(normal).as_matrix()
    return normalized


class Fractions:

    def __init__(self, modifications, abundance, pepsel):

        modifications.append("Input")
        for mod in modifications:
            notlist = list(np.array(modifications)[np.array([mod != i for i in modifications])])
            self.__dict__[mod] = omin.ModificationDelinate(abundance, mod, notlist)

        modifications.remove("Input")

        for mod in modifications:
            self.__dict__[mod].addAttribute("load_normalized", normalizeTo(self.__dict__[mod].abundance, self.Input.abundance))
            self.__dict__[mod].addAttribute("load_norm_log", Logger(self.__dict__[mod].load_normalized))
            self.__dict__[mod].addAttribute("filtered",self.__dict__[mod].load_norm_log.log_div_ave.ix[pepsel.index])

    def addAttribute(self, attribute_name, attribute_data):
        self.__dict__[attribute_name] = attribute_data

    def __repr__(self):
        return "Attributes: "+", ".join(list(self.__dict__.keys()))


# === NEW MAIN CLASS ===
class NormalizedToInput(object):
    """
    Attributes
    ----------
    """
    def __init__(self, raw_peptides=None, raw_proteins=None,
                 modifications=None, genotypes=None, treatments=None):

        self.peptides_raw_abundance = SelectionTools.sep(raw_peptides,
                                                         "Abundance:")

        self.proteins_raw_abundance = SelectionTools.sep(raw_proteins,
                                                         "Abundance:")

        fraction_set = None
        try:
            proteins_input_fract = SelectionTools.sep(self.proteins_raw_abundance, "[Ii]nput")
            self.proteins_input_fract = proteins_input_fract

            #Find all "Fn:" in the raw_abundance where is the fraction number.
            flist = self.peptides_raw_abundance.columns.str.findall("F\d").tolist()
            # Reduce list of all "Fn:" to set of all "Fn:"
            fraction_set = set([i[0] for i in flist])
            self.fraction_set = fraction_set
            # For each "Fn" in the set of all "Fn:"s
            for fraction in fraction_set:
                # Separate the given fraction.
                peptides_fract = SelectionTools.sep(self.peptides_raw_abundance, fraction)

                normal_fract = normalizeTo(peptides_fract, self.proteins_input_fract)
                normal_fract.columns = "Normalized to Protein Input Fraction, " + normal_fract.columns
                # Create one large dataframe:
                self.__dict__[fraction] = peptides_fract.join(normal_fract)

        except Exception:
            print("No abundance columns contained 'F(number)'.")

    def __repr__(self):
        """Show all attributes.
        """
        return "Attributes: "+", ".join(list(self.__dict__.keys()))
