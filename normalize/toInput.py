# -*- coding: utf-8 -*-
"""Normalization to input classes and functions.

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

import itertools
# from operator import itemgetter
from omin.utils import StringTools
from omin.utils import SelectionTools
from omin.normalize.methods import normFactors
from omin.normalize.methods import normalizeTo
from omin.normalize.methods import Logger
from omin.normalize.methods import MachLink
from IPython.display import display
from ipywidgets import widgets

class NormalizedToInput(object):
    """Claculate the normalized relative abundance and relative occupancy.

    Attributes
    ----------
    peptide_groups : obj:
        An instance of the peptide_groups class.
    """

    def __init__(self, parent_self):
        """Initalize NormalizedToInput class.

        Parameters
        ----------
        parent_self : obj:
            The self argument of the parent class.
        """
        try:
            self.peptide_groups = PeptideGroups(parent_self)

        except Exception:
            print("omin.normalize.toInput.PeptideGroups FAILED")

        try:
            self.proteins = Proteins(parent_self)
        except Exception:
            print("omin.normalize.toInput.Proteins FAILED")

    def __repr__(self):
        """Show all attributes.
        """
        if "_Jupyter" in globals().keys():
            view = lambda Attribute: display(self.__dict__[Attribute])
            vw = widgets.interactive(view, Attribute=list(self.__dict__.keys()))
            display(vw)
            return "Attributes: "+", ".join(list(self.__dict__.keys()))
        else:
            return "Attributes: "+", ".join(list(self.__dict__.keys()))


class PeptideGroups(object):
    """Normalize the abundance of the peptide groups.

    Attributes
    ----------
    raw_abundance : DataFrame
        Unfilter non-normalized abundance values.
    """

    def __init__(self, parent_self=None):
        """Initalize PeptideGroups class.

        Parameters
        ----------
        parent_self : obj:
            The self argument of the parent class.

        """
        # Filter for just the Abundance columns.
        self.raw_abundance = parent_self.raw_peptides.filter(regex="Abundance:")

        # Filter the abundance columns for just input columns.
        self.input_abundance = self.raw_abundance.filter(regex="[Ii]nput")

        # Filter out the input fraction(s)
        negate_term = StringTools.regexNot("[Ii]nput")
        self.ptm_abundance = self.raw_abundance.filter(regex=negate_term)

        # Create a list of the ptm containing fractions.
        ptm_fn = list(SelectionTools.find_fractions(self.ptm_abundance))
        self.ptm_fraction_numbers = ptm_fn
        # Create a list of the input containing fractions
        inp_fn = list(SelectionTools.find_fractions(self.input_abundance))
        self.input_fraction_numbers = inp_fn

        self.combos = list(itertools.product(self.ptm_fraction_numbers,
                                             self.input_fraction_numbers))
        self.combo_dict = dict()
        # SINGLE INPUT METHOD
        if parent_self._input_number == 1:
            normalized = []
            for n in self.combos:
                normalized_df = normalizeTo(self.ptm_abundance.filter(regex=n[0]),
                                            self.input_abundance.filter(regex=n[1]))
                print("Peptide Groups", n[0], "normalized to", n[1])

                normalized.append(normalized_df)
            self.normalized_abundances = normalized
        # MULTIPLE INPUT METHOD
        if parent_self._input_number > 1:
            for i in self.combos:
                linkage = MachLink.column_simillarity(self.ptm_abundance,
                                                      self.input_abundance,
                                                      term_a=i[0],
                                                      term_b=i[1])[0]
                self.combo_dict[i] = linkage

            linked_fractions = MachLink.select_top_linked(self.combo_dict, parent_self._input_number)

            normalized = []
            self.linked_fractions = linked_fractions
            for n in linked_fractions:
                fptm = self.ptm_abundance.filter(regex=n[0])
                finp = self.input_abundance.filter(regex=n[1])
                normalized_df = normalizeTo(fptm, finp)
                normalized.append(normalized_df)
                print(n[0], "normalized to", n[1])
            self.normalized_abundances = normalized

    def __repr__(self):
        """Show all attributes.
        """
        return "Attributes: "+", ".join(list(self.__dict__.keys()))


class Proteins(object):
    """Calculate the relative occupancy.

    Attributes
    ----------
    raw_abundance : DataFrame
    input_abundance : DataFrame
    """

    def __init__(self, parent_self=None):
        """Initalize PeptideGroups class.

        Parameters
        ----------
        parent_self : obj:
            The self argument of the parent class.
        """
        # Filter for just the Abundance columns.
        self.raw_abundance = parent_self.raw_proteins.filter(regex="Abundance:")

        # Filter the abundance columns for just input columns.
        self.input_abundance = self.raw_abundance.filter(regex="[Ii]nput")

    def __repr__(self):
        """Show all attributes."""
        return "Attributes: "+", ".join(list(self.__dict__.keys()))
