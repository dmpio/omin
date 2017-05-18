# -*- coding: utf-8 -*-
"""Tools for making selections on DataFrames"""
# LICENSE
# -------

# Copyright 2017 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach,
# Blair Chesnut, and Elizabeth Hauser.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files, (the software)), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions: The above copyright
# notice and this permission notice shall be included in all copies or
# substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS",
# WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
# TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM. OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import os
import re
import pickle
import itertools
import numpy as np
import pandas as pd
from difflib import SequenceMatcher

from .string_tools import StringTools


this_dir, _ = os.path.split(__file__)

# Load the modifications dictionary.
mod_dict_local = "/databases/mod_dict.p"

# If using windows replace "/" with "\\"
if os.name == "nt":
    mod_dict_local = mod_dict_local.replace("/", "\\")

modification_terms = pickle.load(open(this_dir+mod_dict_local, "rb"))


class SelectionTools(object):
    """
    """

    @staticmethod
    def sep(dataframe_in, search_term, strict=False,
            match=False, reverse=False):
        """DEPRICATED in favor of pandas df.filter
        Takes DataFrame and search_term and returns a new DataFrame that
        contains columns that contain that search_term.

        Parameters
        ----------
        dataframe_in : DataFrame
        search_term : str
            What kind of columns do you want to find.
        strict : bool
            Defaults to False. FIXME : Annotate this.
        match : bool
            Defaults to False. If the function will use pandas
            dataframe.columns.str.match which is more strict than
            dataframe.columns.str.search.
        Returns
        -------
        dataframe_out : DataFrame
            A DataFrame that with columns contain just search_term.

        Examples
        --------
        >>>omin.sep(mydataframe,"Search Term")
        dataframe_out
        >>>omin.sep(mydataframe,"Search Term",match=True)
        Scricter dataframe_out

        See Also
        --------
        omin.sepCon
        omin.betSep
        omin.modSel
        omin.manyModSel

        """
        dataframe_out = None
        if match:
            dataframe_out = dataframe_in[dataframe_in.columns[
                dataframe_in.columns.str.match(search_term)]].copy()
            return dataframe_out

        if dataframe_in.columns[
                    dataframe_in.columns.str.contains(search_term,
                                                      case=False)].any():
            if strict:

                dataframe_out = dataframe_in[dataframe_in.columns[
                    np.array([search_term in set(re.sub(" ", "", i).split(",")) for i in dataframe_in.columns])]]
            else:
                dataframe_out = dataframe_in[
                    dataframe_in.columns[dataframe_in.columns.str.contains(search_term, case=False)]].copy()
        else:
            print("The DataFrame has no columns that contain:", search_term)

        if reverse:
            dataframe_out = dataframe_in[
                dataframe_in.columns[~dataframe_in.columns.str.contains(search_term,)]].copy()
        return dataframe_out

    @staticmethod
    def filterRow(dataframe_in=None, on=None, term=None):
        """Return DataFrame that contains a given term on specific column.
        Parameters
        ----------
        dataframe : DataFrame
        term : str
        on : str

        Returns
        -------
        filtered : DataFrame
        """
        filtered = dataframe_in[dataframe_in[on].str.contains(term)]
        return filtered

    @staticmethod
    def modSel(dataframe_in, mod1, mod2):
        """Selects peptides with modifications

        Parameters
        ----------
        dataframe_in : DataFrame
        mod1 : str
        mod2 : str

        Returns
        -------
        dmod12 : DataFrame
            Contains peptides with both mod1 and mod2.
        dmod1 : DataFrame
            Contains peptides with mod1
        dmod2 : DataFrame
            Contains peptides with mod1

        See Also
        --------
        omin.sep
        omin.sepCon

        """
        dmod1 = dataframe_in[dataframe_in.Modifications.str.contains(mod1)]
        dmod2 = dataframe_in[dataframe_in.Modifications.str.contains(mod2)]
        dmod12 = dataframe_in[dataframe_in.Modifications.str.contains(
            mod1) ^ dataframe_in.Modifications.str.contains(mod2, case=False)]
        return dmod12, dmod1, dmod2

    @staticmethod
    def sepCon(dataframe_in, separation_term):
        """Separates dataframe_in two DataFrames one with the separation_term and
        one without the separation_term.

        Parameters
        ----------
        dataframe : DataFrame
        separation_term : str

        Returns
        -------
        with_term = DataFrame
        without_term = DataFrame

        See Also
        --------
        omin.sep
        omin.betSep
        omin.modSel
        omin.manyModSel
        omin.specSel

        """
        with_term = dataframe_in[dataframe_in.columns[
            dataframe_in.columns.str.contains(separation_term, case=False)]]
        without_term = dataframe_in[dataframe_in.columns[
            ~dataframe_in.columns.str.contains(separation_term, case=False)]]
        return with_term, without_term

    @staticmethod
    def betSep(dataframe_in=None, *args):
        """A better version of the function sepCon.

        Parameters
        ----------
        dataframe_in : DataFrame
        *args : str

        Returns
        -------
        with_terms : DataFrame
        without_terms : DataFrame

        """
        sel = np.array([dataframe_in.columns.str.contains(i, case=False)
                        for i in args])
        sel = np.any(sel, axis=0)
        with_terms = dataframe_in[dataframe_in.columns[sel]]
        without_terms = dataframe_in[dataframe_in.columns[~sel]]
        return with_terms, without_terms

    @staticmethod
    def specSel(dataframe=None, include_list=None,
                exclude_list=None, case=True):
        """Select columns whose headers contain items just the items you want.

        Parameters
        ----------
        dataframe : DataFrame
        include_list : list
        exclude_list : list

        Returns
        -------
        selected : DataFrame

        """
        allbool = np.ones(len(dataframe.columns), dtype=bool)
        for mod in include_list:
            trubool = dataframe.columns.str.contains(mod, case)
            allbool = allbool & trubool
        for ex in exclude_list:
            fakbool = dataframe.columns.str.contains(ex, case)
            allbool = allbool & ~fakbool

        if allbool.any():
            selected = dataframe[dataframe.columns[allbool]]
        else:
            print("Special select failed. Try something else.")
            selected = np.nan
        return selected

    @staticmethod
    def manyModSel(pepdf=None, terms=None, verbose=False):
        """Returns searched peptide a tuple of searched DataFrames with [a] given
        modification(s).

        Parameters
        ----------
        pepdf : DataFrame
            With peptides information.
        terms : list
            Can be any number of modifications as a string. Case does not
            matter and regex special characters can be used e.g.
            'acetyl', 'Phospho',hydroxy...methyl.glutaryl,'ect'

        Returns
        -------
        selected : tuple
            Entering more than one term last element of the tuple will contain
            all modified peptides.

        """
        selected = ()
        for term in terms:
            # term = omin.mod_dict[term]
            term = modification_terms[term]
            moddex = pepdf.Modifications.str.contains(pat=term, case=False)
            if moddex.sum() > 0:
                selected += (pepdf.ix[moddex],)
                if verbose:
                    print(moddex.sum(), "peptides with",
                          term, "modification found.")
            else:
                if verbose:
                    print("No peptides with", term, "modification were found.")
                pass
        if len(selected) > 1:
            all_select = np.bitwise_or.reduce([df.index for df in selected])
            all_selected = pepdf.ix[all_select]
            selected = selected + (all_selected,)
        return selected

    @classmethod
    def sevSel(cls, dataframe=None, term_list=None, match=False):
        """
        Parameters
        ----------
        dataframe : DataFrame
        term_list : list

        Returns
        -------
        dataframe_out : DataFrame

        """
        if type(term_list) == str:
            term_list = [term_list]
        if match:
            dataframe_out = pd.concat(
                [cls.sep(dataframe, term, match=True) for term in term_list],
                axis=1)
            return dataframe_out
        else:
            dataframe_out = pd.concat([cls.sep(dataframe, term)
                                       for term in term_list], axis=1)
            return dataframe_out

    @classmethod
    def colSelPro(cls, dataframe=None, term_list=None):
        """Returns the dataframe of the columns specified in the terms_list.

        For each term in term_list an exact match is tried before excepting
        when omin.sep method is used.

        Parameters
        ----------
        dataframe : DataFrame
        term_list : list

        Returns
        -------
        out_dataframe : DataFrame

        """
        df_list = []
        for term in term_list:
            try:
                selected_col = dataframe[term]
                df_list.append(selected_col)

            except Exception:
                selected_col = cls.sep(dataframe, term)
                df_list.append(selected_col)

        out_dataframe = pd.concat(df_list, axis=1)
        return out_dataframe

    @staticmethod
    def filter_out(dataframe, regex):
        """Return a dataframe with [a] regex(es) filtered out_dataframe.

        Parameters
        ----------
        dataframe : DataFrame

        regex : str or list

        Returns
        -------
        result : DataFrame
        """
        # Create a single negate term out of regex.
        negate_term = StringTools.regexNot(regex)
        # Use the pandas filter method to filter out columns with regex.
        result = dataframe.filter(regex=negate_term)
        return result

# === MODIFICATION ISOLATION/CLASSIFICATION TOOLS ===

    @classmethod
    def simplifyModifications(cls, dataframe):
        """Return a series of simplifed modifications.

        Parameters
        ----------
        df : DataFrame

        Returns
        -------
        simplified : Series
        """
        simplified = None
        if any(dataframe.columns == "Modifications"):
            # Compile the regular expression to be used in findall
            rx = re.compile("x(\w+)\s")
            simplified = dataframe.Modifications.apply(rx.findall)
        else:
            print("omin.utils.SelectionTools.simplifyModifications FAILED")
        return simplified

    @classmethod
    def findModifications(cls, df):
        """Return set of ALL of TYPES of modifications found.

        Parameters
        ----------
        df : Dataframe

        Returns
        -------
        found : set
        """
        found = set(
            itertools.chain.from_iterable(cls.simplifyModifications(df))
            )

        return found

    @classmethod
    def findInVivoModifications(cls, df):
        """Return list of the in vivo modifications.

        Parameters
        ----------
        df : Dataframe

        Returns
        -------
        invivo_modifications : set
        """

        present_modifications = cls.findModifications(df)
        chemical_modifications = {'Oxidation', 'Carbamidomethyl',
                                  'TMT6plex', 'TMT10plex'}
        invivo_modifications = [modification for modification in present_modifications if modification not in chemical_modifications]
        return invivo_modifications

    @staticmethod
    def find_fractions(dataframe):
        """Return list of fractions numbers from a given DataFrame.

        Parameters
        ----------
        dataframe : DataFrame

        Returns
        -------
        res : set
        """
        # Declare res as empty set.
        res = {}
        if any(dataframe.columns.str.contains("Abundance:")):
            abundance = dataframe.filter(regex="Abundance:")
            res = set([i[0] for i in abundance.columns.str.findall("F\d").tolist()])
        else:
            print("No columns in this DataFrame begin with; Abundance:")

        return res

    @classmethod
    def find_plex_number(cls, dataframe):
        """Return the TMT plex number as an int.

        Counts the number of fractions and then uses that number to devide by
        the number of abundance columns.

        Parameters
        ----------
        dataframe : DataFrame

        Returns
        -------
        plex_number : int
        """
        fraction_number = len(cls.find_fractions(dataframe))
        abundance_number = dataframe.filter(regex="Abundance:").shape[1]
        plex_number = abundance_number/fraction_number
        return plex_number

    @classmethod
    def find_number_input(cls, dataframe):
        """Return number of inputs as a float.

        Parameters
        ----------
        dataframe : DataFrame

        Returns
        -------
        number_input : float
        """
        number_input = 0.0
        try:
            plex_number = cls.find_plex_number(dataframe)
            # FIXME: Proof this against Inputase ect.
            # To do that that load the regex with ~ [Ii]nput[\w[punctuation]]
            # FIXME: Redundant filtering of raw DataFrames = memory leak
            input_abundance = dataframe.filter(regex="Abundance:").filter(regex="[Ii]nput")
            number_input = input_abundance.shape[1]/plex_number
        except Exception:
            print("utils.SelectionTools.find_number_input failed")
        return number_input
