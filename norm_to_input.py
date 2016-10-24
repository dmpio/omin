# -*- coding: utf-8 -*-

def normFactors(peptide_data):
    """Takes peptide abundance data and returns normalization factors.

    Normalization factors are derived by the taking the sum of each column in DataFrame then dividing each sum by the
    mean of all the sums.

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

    The 'different' DataFrame is normalized to the 'normal' DataFrame using the normFactors function.

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
        logged = normalized_data.apply(np.log2)
        logged.columns = "Log2 " + logged.columns
        logged.index = normalized_data.index
        self.log2 = logged
        ave = self.log2.mean(axis=1)
        ave = pd.DataFrame(ave, columns=["Ave"])
        ave.index = normalized_data.index
        self.ave = ave
        log_div_ave = pd.DataFrame(self.log2.values - self.ave.values)
        log_div_ave.columns = "Log2-AVE " + normalized_data.columns
        log_div_ave.index = ave.index
        self.log_div_ave = log_div_ave
