# -*- coding: utf-8 -*-
"""
DEPRECATED
----------

Use built-in pandomics normalization functions.

"""

import re
import pandas as pd
import numpy as np
# from scipy.stats import ttest_ind
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
                # FIXME: use the method that has been adapted in pandas.
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
