# -*- coding: utf-8 -*-
"""
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
from ..utils import StringTools
from ..utils import SelectionTools


def logNormToAve(pepdf):
    """Takes a DataFrame composed of a fraction of peptide abundances and then
    subtracts each element in each row by the sum of it's row.

    Notes
    -----
    FIXME : Make sure this method could possibly handle protein data aswell.

    Parameters
    ----------
    pepdf : DataFrame
        Limit to a single fraction of peptide abundance data.

    """
    log2 = pepdf.copy().apply(np.log2)
    ave = log2.mean(axis=1)
    log2_div_ave = log2.sub(ave, axis=0)
    log2_div_ave.columns = "Log2-AVE: " + log2_div_ave.columns
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
    # pool = omin.betSep(log2_div_ave, "control", "pool")[0]
    pool = SelectionTools.sep(log2_div_ave, "[Cc]ontrol|[Pp]ool")
    lda_div_pool = log2_div_ave.sub(pool.ix[:, 0], axis=0)
    # Rename the columns to reflect the operations on them.
    lda_div_pool.columns = [re.sub("Log2-AVE", "Log2-AVE-Pool", i) for i in lda_div_pool.columns]
    return lda_div_pool

# === CLASSES DEFINED HERE ===


class FracParse:
    """
    Attributes
    ----------
    abundance : DataFrame
    log_div_ave : DataFrame
    pool_normalized : DataFrame
    [selected conditions] : DataFrame
        These are created dynamically from the list the user selects.
    """

    def __init__(self, abundance, treatment):
        """
        Parameters
        ----------
        abundance : DataFrame
        treatment : list
        """
        self.abundance = abundance
        self.log_div_ave = logNormToAve(self.abundance)
        self.pool_normalized = normToPool(self.log_div_ave)

        for select in treatment:
            # Remove numbers and spaces from selected term
            # term = phraseWasher(select,number_separator="_",word_separator="_").lower()
            term = StringTools.phraseWasher(select, number_separator="_",
                                            word_separator="_").lower()
            # term = re.sub(" ", "_", select)
            self.__dict__[term] = omin.sep(self.pool_normalized, select)

    def addAttribute(self, attribute_name, attribute_data):
        self.__dict__[attribute_name] = attribute_data

    def __repr__(self):
        return "Attributes: "+", ".join(list(self.__dict__.keys()))


class PoolMod:
    """

    Attributes
    ----------
    abundance : DataFrame
    [genotype] : (:obj)
        FracParse object.
    """
    def __init__(self, abundance, mod, genotypes, treatment):
        self.abundance = omin.sep(abundance, mod)
        for geno in genotypes:
            self.__dict__[geno] = FracParse(omin.sep(self.abundance, geno),
                                            treatment)

    def addAttribute(self, attribute_name, attribute_data):
        self.__dict__[attribute_name] = attribute_data

    def __repr__(self):
        return "Attributes: "+", ".join(list(self.__dict__.keys()))

# === OLD MAINCLASS ===


class WithPool(object):
    """
    Attributes
    ----------
    raw : DataFrame
        The raw DataFrame with all information.
    abundance : DataFrame
        Abundance columns from raw DataFrame.
    """

    def __init__(self, raw, modifications, genotypes, treatment):

        """

        Parameters
        ----------
        raw: DataFrame

        """
        self.raw = raw
        self.abundance = omin.sep(raw, "Abundance:")

        for mod in modifications:
            self.__dict__[omin.mod_dict[mod]] = PoolMod(self.abundance, mod,
                                                        genotypes, treatment)

    def addAttribute(self, attribute_name, attribute_data):
        self.__dict__[attribute_name] = attribute_data

    def __repr__(self):
        return "Attributes: "+", ".join(list(self.__dict__.keys()))

# === NEW MAIN CLASS ===


class NormalizedToPool(object):
    """
    Attributes
    ----------
    raw_abundance: DataFrame
    """
    def __init__(self, raw_peptides=None, modifications=None,
                 genotypes=None, treatments=None):
        """
        Parameters
        ----------
        raw_peptides: DataFrame
        modifications: list
        genotypes: list
        treatments: list
        """
        self.raw_abundance = SelectionTools.sep(raw_peptides, "Abundance:")

        fraction_set = None
        try:
            # Find all "Fn:" in the raw_abundance where is the fraction number.
            flist = self.raw_abundance.columns.str.findall("F\d").tolist()
            # Reduce list of all "Fn:" to set of all "Fn:"
            fraction_set = set([i[0] for i in flist])
            self.fraction_set = fraction_set
            # For each "Fn" in the set of all "Fn:"s
            for fraction in fraction_set:
                # Separate the given fraction.
                fract_abundance = SelectionTools.sep(self.raw_abundance, fraction)
                # Log2 normalize
                fract_log = logNormToAve(fract_abundance)
                # Normalize to the pool of the given fraction.
                fract_norm = normToPool(fract_log)
                # Create one large dataframe:
                # raw_abundance + log2(raw_abundance)/AVE(raw_abundance) + log2(raw_abundance)/AVE(raw_abundance)/Pool
                self.__dict__[fraction] = fract_abundance.join(fract_log).join(fract_norm)

        except Exception:
            print("No abundance columns contained 'F(number)'.")

    def __repr__(self):
        """Show all attributes.
        """
        return "Attributes: "+", ".join(list(self.__dict__.keys()))
