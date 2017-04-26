# -*- coding: utf-8 -*-
""" omin's normalization methods.

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

import difflib
import itertools
import pandas as pd
import numpy as np
from operator import itemgetter
# === Pool Methods ===

def logNormToAve(dataframe, add_to_header=None):
    """Return log2(datafame)-mean(log2(datafame)).

    Parameters
    ----------
    datafame : DataFrame

    Returns
    -------
    log2_div_ave : DataFrame

    """
    add_to_header = add_to_header or "Log2-Mean(Log2): "

    log2 = dataframe.copy().apply(np.log2)
    ave = log2.mean(axis=1)
    log2_div_ave = log2.sub(ave, axis=0)
    log2_div_ave.columns = add_to_header + log2_div_ave.columns
    return log2_div_ave


def normToPool(log2_div_ave):
    """Normalizes a DataFrame composed of a single fraction of peptide data
    that contains a single 'Control' or 'Pool' column.

    Notes
    -----
    See also: logNormToAve

    Parameters
    ----------
    log2_div_ave : DataFrame

    """
    # FIXME: This needs to be rethought.
    pool = omin.betSep(log2_div_ave, "control", "pool")[0]
    lda_div_pool = log2_div_ave.sub(pool.ix[:, 0], axis=0)
    # Rename the columns to reflect the operations on them.
    lda_div_pool.columns = [re.sub("Log2-AVE", "Log2-AVE-Pool", i) for i in lda_div_pool.columns]
    return lda_div_pool


# === Input Methods ===
def normFactors(peptide_data):
    """Take peptide abundance data and returns normalization factors.

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
    """Normalize one DataFrame to another.

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
    # Sort the different dataframe by it's columns.
    different_sort = different.columns.tolist()
    different_sort.sort()
    different = different[different_sort]
    # Sort the normal dataframe by it's columns.
    normal_sort = normal.columns.tolist()
    normal_sort.sort()
    normal = normal[normal_sort]
    # Divide the different df by the normalization factors.
    normalized = different / normFactors(normal).as_matrix()
    normalized.columns = different.columns + ": Normalized to: " + normal.columns
    return normalized


# MACHINE LEARNING LINKAGE METHODS
class MachLink(object):
    """Machine learning based DataFrame linkage methods."""

    @staticmethod
    def simillarity(string_one, string_two):
        """Return the difflib.SequnceMatcher results for two strings as ratio.

        Parameters
        ----------
        srting_one : str

        string_two : str

        Returns
        -------
        results : float
        """
        results = difflib.SequenceMatcher(None, string_one, string_two).ratio()
        return results

    @classmethod
    def column_simillarity(cls, dataframe_a, dataframe_b, term_a, term_b):
        """Return a list of alikeness coefficients.

        Parameters
        ----------
        dataframe_a : DataFrame
        dataframe_b : DataFrame
        term_a : str
        term_b : str

        Returns
        -------
        cof : list
        """
        # Get list of columns from dataframe_a filtered by term_a.
        filtered_a = dataframe_a.filter(regex=term_a).columns
        # Get list of columns from dataframe_b filtered by term_b.
        filtered_b = dataframe_b.filter(regex=term_b).columns
        # Create list of alikeness coefficients.
        cof = list(set(itertools.starmap(cls.simillarity,
                                         zip(filtered_a, filtered_b))))
        return cof

    @staticmethod
    def select_top_linked(combo_dict, numb):
        """Select a given number of keys for a combo_dict.

        FIXME : This function needs to be rethought at the moment it just
        returns all of the combinations.

        Parameters
        ----------
        combo_dict : dict
        numb : int

        Returns
        -------
        linked_fractions : list
        """
        linked = None
        try:
            numb = int(numb)
            snumb = len(combo_dict)-numb
            linked = list(dict(sorted(combo_dict.items(),
                                      key=itemgetter(1))).keys())
        except Exception:
            print("omin.normalize.methods.MachLink.select_top_linked FAILED")
        return linked


# THE LOGGER CLASS
class Logger:
    """The class Logger preforms operations on normalized peptides.

    Attributes
    ----------
    log2 : DataFrame
        Containing the log2 of the normalized data
    ave : DataFrame
        Containing the average of the aformentioned DataFrame
    log_div_ave : DataFrame
        Containing the log2-ave.
    """

    def __init__(self, normalized_data):
        """
        Parameters
        ----------
        normalized_data : DataFrame
            Takes normalized peptide or protein in a DataFrame
        """
        # Take the log2 of the data.
        logged = normalized_data.apply(np.log2)
        # Change the column names.
        logged.columns = "Log2 " + logged.columns
        # Make the logged index the same as normalized_data
        logged.index = normalized_data.index
        # Create the log2 attribute from logged.
        self.log2 = logged
        # Take the mean of the logged data.
        ave = self.log2.mean(axis=1)
        # Name the new DataFrame's column
        ave = pd.DataFrame(ave, columns=["Ave"])
        ave.index = normalized_data.index
        # Create the ave attribute from ave
        self.ave = ave
        # Divide the log2 values by their mean (subtraction is division in
        # log-space).
        log_div_ave = pd.DataFrame(self.log2.values - self.ave.values)
        log_div_ave.columns = "Log2-AVE " + normalized_data.columns
        log_div_ave.index = ave.index
        # Create the log_div_ave attribute.
        self.log_div_ave = log_div_ave

    def __repr__(self):
        return "Attributes: "+", ".join(list(self.__dict__.keys()))
