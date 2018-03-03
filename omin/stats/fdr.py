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
# FIXME: ADD DEPRECATION WARNING TO ALL FUNCTIONS AND classes

# import warnings
# import pandas as pd
from pandomics import pandas as pd
# from ..utils.warnings import deprecated


class Comparison:
    """
    This code has been adapted from Bonferioni Calculator v1.1.py. For the full
    version visit the following website:
    https://sites.google.com/site/christophernaugler/open-source-software-for-performing-bonferroni-and-related-corrections
    """
    def __init__(self, entry_number, p_value):
        self.entry_number = entry_number
        self.p_value = p_value
        self.p_corrected = 0
        self.significant = 0

    def __repr__(self):
        return repr((self.entry_number, self.p_value))


def benjamini_hochberg(comparison_objects, sig_level):
    """
    This code has been adapted from Bonferioni Calculator v1.1.py. For the full
    version visit the following website:
    https://sites.google.com/site/christophernaugler/open-source-software-for-performing-bonferroni-and-related-corrections
    """

    sequential_comparison_objects = sorted(comparison_objects,
                                           key=lambda comparison_value: comparison_value.p_value,
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


# @deprecated
def bhFDR(pval_series):
    """Return the p-adjusted values for a series of p-values.

    Parameters
    ----------
    pval_series : Series

    Returns
    -------
    p_adjust : DataFrame

    """
    print("DEPRECATED")
    # Test to see if pval_series is actually a DataFrame.
    if type(pval_series) == pd.core.frame.DataFrame:
        # If it is a DataFrame test to make sure that it only has one column.
        if pval_series.shape[1] == 1:
            # If it does has just one columns slice it from the DataFrame to create a Series.
            pval_series = pval_series.ix[:,0]
        else:
            print("Please format your p-values as a pandas series and retry.")

    # Test to see if pval_series is Series.
    if type(pval_series) == pd.core.series.Series:
        # Test for NANs
        if pval_series.isnull().any():
            og_pval = pval_series
            pval_series = pval_series.dropna()

    comparison_objects = []
    for n, v in enumerate(pval_series):
        # Create and append comparison objects for each value.
        comparison_objects.append(Comparison(n, v))

    # Finally run through BH correction.
    bhc = benjamini_hochberg(comparison_objects, .05)
    p_adjust = pd.DataFrame([bhc[i].p_corrected for i in range(len(bhc))],
                            columns=["p-adjusted"],
                            index=pval_series.index)

    # Test to see if og_pval exists in local namespace
    if "og_pval" in locals():
        # If it does then reindex p_adjust to og_pval.index
        p_adjust = p_adjust.reindex(index=og_pval.index)

    return p_adjust
