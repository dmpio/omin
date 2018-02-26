# -*- coding: utf-8 -*-
"""Normalization to input classes and functions.

DEPRECATED
----------
Use built-in pandomics normalization functions.

"""

import difflib
import itertools
from operator import itemgetter
from ..utils import StringTools
from ..utils import SelectionTools
from ..core import Handle


class MachLink(object):
    """Machine learning based DataFrame linkage methods.

    DEPRECATED
    ----------
    Use better idea.
    """

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


class NormalizedToInput(Handle):
    """Claculate the normalized relative abundance and relative occupancy.

    DEPRECATED
    ----------
    Use omin.core.containers

    Attributes
    ----------
    peptide_groups : obj:
        An instance of the peptide_groups class.
    """
    def __init__(self, parent_self):
        """Initalize NormalizedToInput class.
        """
        Handle.__init__(self)

    # def __init__(self, parent_self):

        try:
            self.peptide_groups = PeptideGroups(parent_self)

        except Exception:
            print("omin.normalize.toInput.PeptideGroups FAILED")

        try:
            self.proteins = Proteins(parent_self)
        except Exception:
            print("omin.normalize.toInput.Proteins FAILED")


class PeptideGroups(object):
    """Normalize the abundance of the peptide groups.

    DEPRECATED
    ----------
    Use omin.core.containers.PeptideGroups

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
                # normalized_df = normalizeTo(self.ptm_abundance.filter(regex=n[0]),
                #                             self.input_abundance.filter(regex=n[1]))
                normalized_df = self.ptm_abundance.filter(regex=n[0]).normalize_to(self.input_abundance.filter(regex=n[1]))
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
                # normalized_df = normalizeTo(fptm, finp)
                normalized_df = fptm.normalize_to(finp)
                normalized.append(normalized_df)
                print(n[0], "normalized to", n[1])
            self.normalized_abundances = normalized

    def __repr__(self):
        """Show all attributes.
        """
        return "Attributes: "+", ".join(list(self.__dict__.keys()))


class Proteins(object):
    """Calculate the relative occupancy.

    DEPRECATED
    ----------
    Use omin.core.containers.Proteins

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
