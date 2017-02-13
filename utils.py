# -*- coding: utf-8 -*-
import re
import os
import pickle
import numpy as np
import pandas as pd

from datetime import datetime
from urllib.request import urlopen
from urllib.error import HTTPError

this_dir, _ = os.path.split(__file__)

# Load the modifications dictionary.
modification_terms = pickle.load(open(this_dir+"\databases\mod_dict.p", "rb"))


# FIXME: Create tools to search against against databases the user specifies.
# === UNIPROT TOOLS ===

def inspectObject(obj):
    """Returns a list of object attributes and their types.

    Parameters
    ----------
    obj: :obj:
        Any kind of object that is not a builtin.

    Returns
    -------
    obj_ids: list
        A list of lists to be specific.
    """
    obj_ids = []
    # Make list of builtins
    bi_list = dir(__builtins__)
    # Try to make a list of the types of things inside obj.
    try:
        # Make list all things inside of an object
        # obj_ids = [[name, type(thing).__name__] for name, thing in obj.__dict__.items()]

        for name, thing in obj.__dict__.items():
            obj_ids.append([name, type(thing).__name__, thing])

        return obj_ids

    except Exception:
        pass


def objectWalker(obj, desired_type=None, att_list=None):
    """Recursively walks through an object an genrates a list of its attributes.

    Parameters
    ----------
    obj: :obj:

    desired_type: str

    Returns
    -------
    att_list: list

    """
    # Create an attributes list if none has been asigned.
    if att_list is None:
        att_list = []

    # Create list of types that we do not want yo serch through
    no_list = dir(__builtins__)
    no_list.append("DataFrame")
    # desired_type = desired_type or "DataFrame"

    obj_inspect = inspectObject(obj)

    for i in obj_inspect:

        # Test if object attribute is desired_type.
        if re.search(desired_type, i[1], re.I) is not None:

            att_list.append(i)

        # If the object attribute is not the desired_type then...
        elif i[1] not in no_list:
            # Try to use recursively objectWalker on the object.
            try:
                objectWalker(obj.__dict__[i[0]], desired_type, att_list)
            except:
                pass
    return att_list


class FastaTools(object):
    """A collection of tools for handling FASTA formatted strings.
    """

    @classmethod
    def fasta2Seq(cls, fasta):

        """ Joins strings in list into one string. Neglecting the first which is
        generally the annotation.

        Parameters
        ----------
        fasta : str

        Returns
        -------
        one_string : str
        """
        one_string = "".join(fasta[1:])
        return one_string

    @staticmethod
    def seqNum(cls, seq):
        """Takes a sequence and returns the number for each amino acid (I think).

        Notes
        -----
            Review this functions useage and rewrite annotation.

        Parameters
        ----------
        seq : str

        Returns
        -------
        seq_num : str

        """
        seq_num = [seq[i]+str(i+1) for i in range(len(seq))]
        return seq_num

# === UNIPROT TOOLS ===


class UniprotTools(object):

    @classmethod
    def getFasta(cls, upID):
        """Takes a uniprot ID and returns a string.

        Parameters
        ----------
        upID : str

        Returns
        -------
        fasta : str

        """
        try:
            fasta = str
            link = "http://www.uniprot.org/uniprot/"+upID+".fasta"
            try:
                urlout = urlopen(link)
                # print("hey")

            except HTTPError:
                print(
                      "Try again, the formatting may be off.",
                      "Try removing dashes."
                      )

            # fasta = urlout.read().decode("utf-8")
            fasta = urlout.read().decode("utf-8").split("\n")
            return fasta
        except TypeError as err:
            print("Please enter UniProt ID as string.")

# === STRING TOOLS ====


class StringTools(object):

    @staticmethod
    def time_stamp():
        """Takes no arguments and returns datetime timestamp as string.

        Returns
        -------
        ts : str

        Examples
        --------
        >>>StingTools.time_stamp()
        '1483631005_38507'
        """
        ts = str(datetime.timestamp(datetime.now())).replace(".", "_")
        return ts

    @staticmethod
    def multiRegExOr(term_list):
        """Returns a singles string of "OR"ed regex terms from a given list.

        Parameters
        ----------
        term_list: list
            List of regular expressions.

        Returns
        -------
        all_exp: str
            A string pattern ready for use with re methods.

        Exaples
        -------
        >>>re.search(multiRegExOr(skynet.modification_terms.values()),
        >>>                       'Phospho (fraction)')

        See Also
        --------
        re.search
        re.match
        """
        all_exp = "|".join(["(%s)" % (term) for term in term_list])
        return all_exp

    @staticmethod
    def rxAndGen(ex_list):
        """Returns a regex psuedo "AND" formatted string from list of strings.

        Parameters
        ----------
        ex_list: list

        Returns
        -------
        formatted_list: str
        """

        rxAnd = lambda x: "(?=.*\\b{}\\b)".format(x)

        andList = list(map(rxAnd, ex_list))

        formatted_list = "^"+"".join(andList)+".*$"

        return formatted_list


    # @classmethod
    @staticmethod
    def int2word(num, separator="-"):
        """Transforms integers =< 999 into english words

        Parameters
        ----------
        num : int
        separator : str

        Returns
        -------
        words : str
        """
        ones_and_teens = {0: "Zero", 1: 'One', 2: 'Two', 3: 'Three',
                          4: 'Four', 5: 'Five', 6: 'Six', 7: 'Seven',
                          8: 'Eight', 9: 'Nine', 10: 'Ten', 11: 'Eleven',
                          12: 'Twelve', 13: 'Thirteen', 14: 'Fourteen',
                          15: 'Fifteen', 16: 'Sixteen', 17: 'Seventeen',
                          18: 'Eighteen', 19: 'Nineteen'}
        twenty2ninety = {2: 'Twenty', 3: 'Thirty', 4: 'Forty', 5: 'Fifty',
                         6: 'Sixty', 7: 'Seventy', 8: 'Eighty', 9: 'Ninety', 0: ""}

        if 0 <= num < 19:
            return ones_and_teens[num]
        elif 20 <= num <= 99:
            tens, below_ten = divmod(num, 10)
            if below_ten > 0:
                words = twenty2ninety[tens] + separator + \
                    ones_and_teens[below_ten].lower()
            else:
                words = twenty2ninety[tens]
            return words

        elif 100 <= num <= 999:
            hundreds, below_hundred = divmod(num, 100)
            tens, below_ten = divmod(below_hundred, 10)
            if below_hundred == 0:
                words = ones_and_teens[hundreds] + separator + "hundred"
            elif below_ten == 0:
                words = ones_and_teens[hundreds] + separator + \
                    "hundred" + separator + twenty2ninety[tens].lower()
            else:
                if tens > 0:
                    words = ones_and_teens[hundreds] + separator + "hundred" + separator + twenty2ninety[
                        tens].lower() + separator + ones_and_teens[below_ten].lower()
                else:
                    words = ones_and_teens[
                        hundreds] + separator + "hundred" + separator + ones_and_teens[below_ten].lower()
            return words

        else:
            print("num out of range")


    @classmethod
    def phraseWasher(cls, phrase, number_separator="_", word_separator=" "):
        """Replaces numerical portions of strings with words.

        Parameters
        ----------
        phrase : str
        number_separator : str
            separator to be used with int2word function.
        word_separator : str
            Defaults to " " but you can put in whatever you like.

        See Also
        --------
        int2word

        Examples
        --------
        >>>phraseWasher("10 min post",word_separator = "_")
        "Ten_min_post"
        >>>phraseWasher("65 min post",number_separator = "_", word_separator = "_")
        "Sixty_five_min_post"

        """
        new_phrase = []
        for i in phrase.split():
            if i.isnumeric():
                phrase_part = cls.int2word(int(i), number_separator)
                new_phrase.append(phrase_part)
            else:
                phrase_part = i
                new_phrase.append(phrase_part)
        washed = word_separator.join(new_phrase)
        return washed

# === Selection tools ===


class SelectionTools(object):
    """
    """
    @staticmethod
    def sep(dataframe_in, search_term, strict=False, match=False):
        """Takes DataFrame and search_term and returns a new DataFrame that
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
        return dataframe_out

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
    def specSel(dataframe=None, include_list=None, exclude_list=None, case=True):
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
            Can be any number of modifications as a string. Case does not matter
            and regex special characters can be used e.g.
            'acetyl', 'Phospho',hydroxy...methyl.glutaryl,'ect'

        Returns
        -------
        selected : tuple
            Entering more than one term last element of the tuple will contain all
            modified peptides.

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

    @staticmethod
    def sevSel(dataframe=None, term_list=None, match=False):
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
                [omin.sep(dataframe, term, match=True) for term in term_list],
                axis=1)
            return dataframe_out
        else:
            dataframe_out = pd.concat([omin.sep(dataframe, term)
                                       for term in term_list], axis=1)
            return dataframe_out

    @staticmethod
    def colSelPro(dataframe=None, term_list=None):
        """Returns the dataframe of the columns specified in the terms_list.

        For each term in term_list an exact match is tried before excepting when
        omin.sep method is used.

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

            except:
                selected_col = omin.sep(dataframe, term)
                df_list.append(selected_col)

        out_dataframe = pd.concat(df_list, axis=1)
        return out_dataframe

    @staticmethod
    def superGroup(dataframe=None, new_level=None):
        """Returns a multiindexed DataFrame with the top index named new_level.

        Parameters
        ----------
        dataframe : DataFrame
        new_level : str

        Returns
        -------
        out_df : DataFrame

        """
        if type(dataframe.columns) == pd.indexes.base.Index:
            out_df = pd.DataFrame(dataframe.values, index=dataframe.index,
                                  columns=pd.MultiIndex.from_product([[new_level],dataframe.columns]))
            return out_df
        if type(dataframe.columns) == pd.indexes.multi.MultiIndex:
            if len(dataframe.columns.levels[0]) < 1:
                levels = [list(i.values) for i in dataframe.columns.levels]
                levels = [[new_level]] + levels
                out_df = pd.DataFrame(
                    dataframe.values, index=dataframe.index,
                    columns=pd.MultiIndex.from_arrays(levels))
                return out_df
            else:
                levels = [[new_level]] + [list(i.values)
                                          for i in dataframe.columns.levels]
                labels = [list(i) for i in dataframe.columns.labels]
                new_list = list(np.linspace(0, 0, len(labels[-1]), dtype=int))
                labels = [new_list] + labels
                multi = pd.MultiIndex(levels, labels)
                out_df = pd.DataFrame(
                    dataframe.values, index=dataframe.index, columns=multi)
                return out_df

    @staticmethod
    def modKindCheck(dataframe, modification):
        """Returns a boolean dataframe True if the given mod is found.

        Parameters
        ----------
        dataframe : DataFrame
        modification : str

        Returns
        -------
        out_dataframe : DataFrame

        """
        series = dataframe.Modifications.str.contains(modification)
        out_dataframe = pd.DataFrame(series)
        out_dataframe.columns = [modification + "-peptide?"]
        return out_dataframe

    # === FILTERING FUNCTIONS ===

    @classmethod
    def masterCleanse(cls, protein_df):
        """Filters raw protein DataFrame for master proteins.

        The raw protein data from Proteome Discoverer there is a column with
        the title 'Master' this funtion scans through that column and selects
        only the proteins that end with the string "IsMasterProtein".

        Parameters
        ----------
        protein_df : DataFrame
            Raw protein DataFrame

        Returns
        -------
        clean : DataFrame
            DataFrame that contains only proteins with 'IsMasterProtein' in
            'Master' column of protein_df.
        """
        clean = protein_df.ix[protein_df.Master.str.endswith("IsMasterProtein")]
        return clean

    @classmethod
    def onePerQ(cls, protein_df):
        """Filters raw protein DataFrame for proteins that are less than 1% the
        expected q-value.

        Scans through the protein DataFrame selecting only the proteins with
        less than 1% of the expected q-value.

        Parameters
        ----------
        protein_df : DataFrame
            Raw protein DataFrame

        Returns
        -------
        clean : DataFrame
            Protein data that contains only proteins with proteins only less
            than 1% of the expected q-value.
        """
        try:
            one_per = protein_df["Exp. q-value"] < .01
            one_per = protein_df.ix[one_per]
            return one_per
        except KeyError:
            print("No 'Exp. q-value' column found. Data has not NOT been filtered at 1% FDR.")
            return protein_df

    @classmethod
    def masterOne(cls, protein_df):

        """Takes a raw protein DataFrame and filters it using first the
        'masterCleanse' function and 'onePerQ' function.

        Parameters
        ----------
        protein_df : DataFrame
            Raw proteins.

        Returns
        -------
        master_one : DataFrame
            Of master proteins with exp. q-value <1%
        """
        master = cls.masterCleanse(protein_df)
        master_one = cls.onePerQ(master)
        return master_one

    @classmethod
    def masterPep(cls, peptide_df):
        """Takes a peptide DataFrame and returns just the first master protein
        accession for each peptide.

        Notes
        -----
        Assumes the first uniprot ID list is the correct one. Peptides with no
        master protein accession will be lost however the index of peptide_df will
        be preserved.

        Parameters
        ----------
        peptide_df : DataFrame

        Returns
        -------
        master_prot_acc : DataFrame

        """
        master_prot_acc = [i.split(';')[0] for i in peptide_df['Master Protein Accessions'].dropna()]

        master_prot_acc = pd.DataFrame(master_prot_acc,
                                       index=peptide_df['Master Protein Accessions'].dropna().index, columns=['Accession'])
        return master_prot_acc

    @classmethod
    def mpaParse(cls, raw_peptides=None, master_uniprot_id="Master",
                 new_column_name="MPA"):
        """Returns a DataFrame containing only the first master protein accession.

        Parameters
        ----------
        raw_peptides : DataFrame
        master_uniprot_id : str
            The search term used in omin.sep(raw_peptides,master_prot_id).
        new_column_name : str
            Label your new dataframe

        Returns
        -------
        mpa : DataFrame

        See Also
        --------
        omin.masterPep

        """
        mpa_list = [i.split(";")[0] if type(i) == str else np.nan for i in omin.sep(raw_peptides, master_uniprot_id).ix[:, 0]]
        mpa = pd.DataFrame(mpa_list, index=raw_peptides.index, columns=[new_column_name])
        return mpa

    @classmethod
    def vLook(cls, peptides=None, proteins=None, mods=None,):
        """Returns a tuple of selected peptides and proteins.

        Takes raw peptides and protiens returns a tuple of selected peptides and
        proteins. The function can also select for a sigle modification or many
        modifications.

        Parameters
        ----------
        peptides : DataFrame
        proteins : DataFrame
        mods : list

        Returns
        -------
        peptide_select : DataFrame
        protein_select : DataFrame

        Examples
        --------
        >>>peptide_select,protein_select = vLook(raw.peptides,raw.proteins)
        >>>peptide_select,protein_select= vLook(raw.peptides,raw.proteins,["hydroxy...methyl.glutaryl"])
        >>>peptide_select,protein_select= vLook(raw.peptides,raw.proteins,["Acetyl","Phospho"])

        See Also
        --------
        manyModSel
        masterOne
        masterPep

        """

        fdr = cls.masterOne(proteins)
        if len(mods) == 0:
            mpa = cls.masterPep(peptides)
        if len(mods) == 1:
            mpa = cls.masterPep(cls.manyModSel(peptides, mods)[0])
        else:
            mpa = cls.masterPep(cls.manyModSel(peptides, mods)[-1])

        fdrdf = pd.DataFrame(fdr.Accession, index=fdr.index)

        peptide_select = mpa.merge(fdrdf, on="Accession",
                                   how="left", right_index=True)

        protein_select = mpa.merge(fdrdf, on="Accession",
                                   how="left", left_index=True)

        return peptide_select, protein_select

    # === VENN DIAGRAM FUNCTIONS ===

    @staticmethod
    def setDiff(list_A, list_B):
        """Takes difference of two lists with respect to list_A. Simillar to
        `list(set(list_A)-set(list_B))` however duplicates are not deleted and
        order is preserved.

        Parameters
        ----------
        list_A : list
        list_B : list

        Returns
        -------
        the_diff :  list
            List of the difference between sets.
        """
        the_diff = [i for i in list_A if i not in set(list_B)]

        return the_diff


    class CompOb:
        """Make a comparison object that mirrors the output of a Venn diagram.
        Attributes
        ----------
        justA : DataFrame
            The elements exclusive to list_A.
        justB : DataFrame
            The elements exclusive to list_B.
        AandB : DataFrame
            The elements found in list_A and list_B.
        """

        def __init__(self, list_A, list_B, how=None,
                     modification_labels=["Acetyl", "Phospho"]):
            """
            Parameters
            ----------
            list_A : list
            list_B : list
            how : str
                Can be compared by any column as a string. If no string is
                specified then the index is used.
            modification_labels : list
                If none is specified then it defaults to  ["Acetyl", "Phospho"]
            """
            if how is None:
                self.justA = list_A.ix[setDiff(list_A.index, list_B.index)]
                self.justB = list_B.ix[setDiff(list_B.index, list_A.index)]
                self.AandB = list_A.ix[set(list_A.index) & set(list_B.index)]
            else:
                self.justA = list_A.ix[list_A[how].isin(
                    setDiff(list_A[how], list_B[how]))]
                self.justB = list_B.ix[list_B[how].isin(
                    setDiff(list_B[how], list_A[how]))]
                combo = pd.concat([list_A, list_B], keys=modification_labels)
                self.AandB = combo.ix[combo[how].isin(
                    set(list_A[how]) & set(list_B[how]))]

    @staticmethod
    def aboveCut(trunch_object, cond, pval_kind="pval", lfc_kind="lfc", cut=.05):
        """Selects p-values above a given cutoff in a given trunch_object and
        returns a DataFrame of p-values and log fold changes.

        Parameters
        ----------
        trunch_object : (:obj)
        cond : str
        pval_kind : str
            If none is specified it defaults to "pval"
        lfc_kind : str
            If none is specified it defaults to "lfc"
        cut : float
            If none is specified it defaults to .05

        Returns
        -------
        above_cut : DataFrame
            Contains p-values and log fold change columns.
        """
        pv = omin.sep(trunch_object.__dict__[pval_kind], cond).ix[
            :, 0].sort_values()
        pv = pv[pv < cut]
        lfc = omin.sep(trunch_object.__dict__[lfc_kind], cond).ix[pv.index]
        above_cut = pd.concat([pv, lfc], axis=1)
        return above_cut

    @staticmethod
    def allComp(trunch_object, modification_object, cond, pval_kind="pval",
                lfc_kind="lfc", cut=.05):
        """
        Notes
        -----
            TODO: Include modification columns

        Parameters
        ----------
        trunch_object : (:obj)
        modification_object : (:obj)
        cond : str
        pval_kind : str
            If none is specified it defaults to "pval"
        lfc_kind : str
            If none is specified it defaults to "lfc"
        cut : float
            If none is specified it defaults to .05

        Returns
        -------
        compound : DataFrame

        """
        # Grab the P-values and LFCs
        pv_lfc = omin.aboveCut(modification_object,
                               cond, pval_kind, lfc_kind, cut)
        # Grab the master protein accession and gene id
        mpa_gn = trunch_object.mpa.ix[pv_lfc.index]
        # Grab the Modifications columns
        mods = trunch_object.peptides.Modifications[pv_lfc.index]
        # grab the MitoCarta2.0 data
        mito = trunch_object.mitopep.ix[pv_lfc.index][
            ["MitoCarta2_List", "Matrix", "IMS"]]
        # Concatenate DataFrames
        compound = pd.concat([mpa_gn, mods, pv_lfc, mito], axis=1)
        # Change NaNs for zeros
        compound = compound.fillna(0)
        return compound

    # === INTERACTIVE SELECTION ===

    @staticmethod
    def treatmentSelect(peptide_abundance):
        """Allows the user to select which columns contain treatments.

        Parameters
        ----------
        peptide_abundance : DataFrame

        Returns
        -------
        select_list : list
        """
        print(pd.DataFrame([i.split(",") for i in peptide_abundance.columns]))
        col_num = int(input("Which column has treatment data?(Enter the number)"))
        select_set = set(pd.DataFrame(
            [i.split(",") for i in peptide_abundance.columns]).ix[:, col_num])

        select_list = [re.sub(" ", "_", i.strip()) for i in select_set]
        # print(select_list)
        [print(n, i) for n, i in enumerate(select_list)]
        remove_num = int(
            input("Enter the number of any element that need to be removed."))
        select_list.remove(select_list[remove_num])
        return select_list

    @staticmethod
    def listWasher(start_list, pattern=None, replace=None):
        """Returns a list that has had the 'pattern' string replaced with the
        'replace' string.

        If no pattern or replace string are entered listWasher will ask the user if
        any of the elements in the string should be replaced and what they should
        be replaced by.

        Parameters
        ----------
        start_list : list
        pattern : str
        replace : str

        Returns
        -------
        out_list : list

        Notes
        -----
        FIXME: Add list support for this function eventually.

        Examples
        --------
        >>>conditions = ['NonEx', 'Immediate Post', '60 min Post', 'n/a']
        >>>new_list = listWasher(conditions,pattern="n/a",replace="Pool")
        >>>print(new_list)
        ['NonEx', 'Immediate Post', '60 min Post', 'Pool']

        Okay let's try with user input.
        >>>new_list = listWasher(conditions)
        Would you like to replace any terms in the list? y/n
        0 NonEx
        1 Immediate Post
        2 60 min Post
        3 n/a
        >>>y
        Which term would you like to replace?
        >>>3
        Replace 'n/a' with?
        >>>Pool
        Would you like to replace any more terms in the list? y/n
        >>>n

        """
        if replace is not None:
            out_list = [replace if x == pattern else x for x in start_list]
            return out_list
        else:
            print(start_list)
            contin = input("Would you like to replace any terms in the list? y/n ")
            out_list = start_list
            while contin == "y":
                for n, element in enumerate(out_list):
                    print(n, element)
                term_num = int(input("Which term would you like to replace? "))
                pattern = out_list[term_num]
                replace = input("Replace" + " '" + pattern + "' " + "with? ")
                out_list = [replace if x == pattern else x for x in out_list]
                contin = input(
                    "Would you like to replace any more terms in the list? y/n ")
            while contin == "n":
                return out_list

# Testing
if __name__ == "__main__":
    print("Testing utils.py ...")
    import pandas as pd
    from omin.core.handles import RawData

    try:
        data = RawData("ExampleData\crat_ex\_E749_4154_010716_PeptideGroups.txt",
                       "ExampleData\crat_ex\_E749_4154_010716_Proteins.txt")
        if data.raw_peptides.shape == (8712, 277):
            print("ExampleData has loaded correctly.")
    except:
        print("Loading ExampleData failed.")
    try:
        if type(inspectObject(data)) == list:
            print("inspectObject works!")
    except:
        print("inspectObject has failed.")

