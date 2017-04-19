# -*- coding: utf-8 -*-
"""Normalization to input classes and functions

LICENSE:
Copyright 2017 James Draper, Paul Grimsrud, Deborah Muoio

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files, Omics Modeling Integrating
Normalization (OMIN), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom
the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM.
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import re
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from omin.utils import StringTools
from omin.utils import SelectionTools

from omin.normalize.methods import *


def normFactors(abundance):
    """Takes peptide abundance data and returns normalization factors.

    Normalization factors are derived by the taking the sum of each column in
    DataFrame then dividing each sum by the mean of all the sums.

    Parameters
    ----------
    abundance : DataFrame

    Returns
    -------
    norm_factors: DataFrame

    """
    norm_factors = abundance.sum() / abundance.sum().mean()
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
    normalized.columns = different.columns + ": Normalized to: " + normal.columns
    return normalized


class NormalizedToInput(object):
    """
    Attributes
    ----------
    """

    def __init__(self, parent_self):

        # self.proteins = None
        try:
            self.peptide_groups = PeptideGroups(parent_self.raw_peptides)

        except Exception:
            print("omin.normalize.toInput.Proteins FAILED")
        # try:
        #     self.proteins = Proteins(raw_proteins)
        # except Exception:
        #     print("omin.normalize.toInput.Proteins FAILED")

    def __repr__(self):
        """Show all attributes.
        """
        return "Attributes: "+", ".join(list(self.__dict__.keys()))


class PeptideGroups(object):
    def __init__(self, raw_peptides=None):
        self.input_fraction_numbers = SelectionTools.find_number_input(raw_peptides)
        self.abundance = raw_peptides.filter(regex="Abundance:")
        self.input = self.abundance.filter(regex="[Ii]nput")
        negate_term = StringTools.regexNot("[Ii]nput")
        self.other_fractions = self.abundance.filter(regex=negate_term)
        self.other_fraction_numbers = list(SelectionTools.find_fractions(self.other_fractions))

class Proteins(object):
    def __init__(self, raw_proteins=None):

        self.input_fraction_numbers = None

        # self.input_fraction_numbers = SelectionTools.find_number_input(raw_proteins)

    def __repr__(self):
        """Show all attributes.
        """
        return "Attributes: "+", ".join(list(self.__dict__.keys()))
