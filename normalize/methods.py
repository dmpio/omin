# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np

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
    pool = omin.betSep(log2_div_ave, "control", "pool")[0]
    lda_div_pool = log2_div_ave.sub(pool.ix[:, 0], axis=0)
    # Rename the columns to reflect the operations on them.
    lda_div_pool.columns = [re.sub("Log2-AVE", "Log2-AVE-Pool", i) for i in lda_div_pool.columns]
    return lda_div_pool

# === Input Methods ===


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
