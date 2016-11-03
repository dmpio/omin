# This code has been adapted from Bonferioni Calculator v1.1.py here is the original documentation:
    #Source code for “Bonferroni Calculator”
    #program requires python 3.2 or higher, available for free at http://www.python.org/getit/releases/3.2/
    #Contact Dr. Naugler at:
    #C414, Diagnostic and Scientific Centre
    #14, 3535 Research Road NW
    #Calgary AB Canada T2L 2K8
    #Ph: 403-770-3756; Email: christopher.naugler@cls.ab.ca
    #NOTE Version 1.1, released March 2013, contains an important fix for the interpretation of Bonferroni-Holm and Benjamini-Hochberg calculations
    #Specifically, in version 1.1, when Benjamini-Hochberg calculations are performed in the rank-ordered list of P-values AFTER the first significant
    #P-value has been declared, any subsequent P-values which calculate to >0.05 are defaulted to equal the largest statistically significant adjusted P-value
    #The resulting calculations are then identical to Benjamini-Hochberg corrections performed with the statistical program "R"
    #
    #For Bonferroni-Holm corrections, version 1.1 also contains an important fix.
    #Specifically, when Bonferroni-Holm corrections are performed in the rank-ordered list of P-values, calculations are stopped after the first non-significant
    #P-value has been declared, and subsequent P-values are defaulted to P=1
    #
    #Version 1.1 also adds the ability to import the list of starting P-values from a text file instead of typing them individually into the program
    #Version 1.1 also adds the interpretation "significant" and "non-significant" after all corrected p-values
    #
    #We wish to thank Nicholas Wiebe, MSc for help in coding these changes
    #
    #Please update to this latest version
    #Information about future updates can be found at: https://sites.google.com/site/christophernaugler/open-source-software-for-performing-bonferroni-and-related-corrections
    #The authors do not assume any liability for errors resulting from the use of this software.
import pandas as pd

class Comparison:
    def __init__(self, entry_number, p_value):
        self.entry_number = entry_number
        self.p_value = p_value
        self.p_corrected = 0
        # version 1.1 added new class variable do designate a non-significant value
        self.significant = 0

    def __repr__(self):
        return repr((self.entry_number, self.p_value))


def benjamini_hochberg(comparison_objects, sig_level):
    # Sorts by p-value

    sequential_comparison_objects = sorted(comparison_objects, key=lambda comparison_value: comparison_value.p_value,
                                           reverse=True)
    rank = 1
    total_measurements = len(comparison_objects)
    rank = total_measurements
    SignificantFound = 0

    lowestPvalue = float("inf")

    for x in sequential_comparison_objects:
        if not SignificantFound:

            x.p_corrected = (x.p_value * total_measurements) / rank
            if x.p_corrected < lowestPvalue:
                lowestPvalue = x.p_corrected
            else:
                x.p_corrected = lowestPvalue
            x.significant = x.p_corrected < sig_level
            if x.significant:
                SignificantFound = 1

        else:
            x.significant = 1
            x.p_corrected = (x.p_value * total_measurements) / rank
            if x.p_corrected < lowestPvalue:
                lowestPvalue = x.p_corrected
            else:
                x.p_corrected = lowestPvalue
        rank = rank - 1

    # Re-sort based on entry order
    sequential_comparison_objects = sorted(comparison_objects,
                                           key=lambda comparison_value: comparison_value.entry_number)
    return sequential_comparison_objects


def bhFDR(pval_series):
    """Returns the p-adjusted values for a series of p-values.

    Note
    ----
    Make sure your p-values Series contain NO NaNs. If you DataFrame contains NaNs the p-adjusted values will be
    incorrect.

    Parameters
    ----------
    pval_series : Series
        Series must contain floats with no NaN values.

    Returns
    -------
    p_adjust : DataFrame

    """
    comparison_objects = []
    for n, v in enumerate(pval_series):
        comparison_objects.append(Comparison(n, v))
    bhc = benjamini_hochberg(comparison_objects, .05)
    p_adjust = pd.DataFrame([bhc[i].p_corrected for i in range(len(bhc))],
                            columns=["p-adjusted"],
                            index=pval_series.index)
    return p_adjust
