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

import numpy as np


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

            except Exception:
                selected_col = omin.sep(dataframe, term)
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

        FIXME: TMT6 and TMT10 are used interchangeablily throughout
        Proteome Discoverer. This whole method needs to be rethought.
        Perhaps count the number of fractions and then use that number to
        devide by the number of abundance columns.

        Parameters
        ----------
        dataframe : DataFrame

        Returns
        -------
        plex_number : int
        """
        # mods = cls.findModifications(dataframe)
        # tmt_type = list(filter(lambda x: x.startswith("TMT"), mods))
        # plex_number = list(filter(lambda x: x.isdigit(), tmt_type[0]))
        # plex_number = int(plex_number[0])

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
            input_abundance = dataframe.filter(regex="Abundance:").filter(regex="[Ii]nput")
            number_input = input_abundance.shape[1]/plex_number
        except Exception:
            print("utils.SelectionTools.find_number_input failed")
        return number_input
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
            print("No 'Exp. q-value' column found.")
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

        Notes
        -----
        FIXME: Use try dropna() and reindex instead of else statement in loop.

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

    @staticmethod
    def alike(string_one, string_two):
        """DEPRECATED use omin.normalize.methods.MachLink.simillarity

        Return difflib.SequnceMatcher.ratio results for two strings.

        Parameters
        ----------
        srting_one : str

        string_two : str

        Returns
        -------
        results : float
        """
        results = SequenceMatcher(None, string_one, string_two).ratio()
        return results

    @classmethod
    def alikeness(cls, dataframe_a, dataframe_b, term_a, term_b):
        """DEPRECATED use omin.normalize.methods.MachLink.column_simillarity

        Return a list of alikeness coefficients.

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
        cof = list(set(itertools.starmap(cls.alike, zip(filtered_a,
                                                        filtered_b))))
        return cof

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

        If no pattern or replace string are entered listWasher will ask the
        user if any of the elements in the string should be replaced and what
        they should be replaced by.

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
            contin = input(
                "Would you like to replace any terms in the list? y/n "
                )
            out_list = start_list
            while contin == "y":
                for n, element in enumerate(out_list):
                    print(n, element)
                term_num = int(input("Which term would you like to replace? "))
                pattern = out_list[term_num]
                replace = input("Replace" + " '" + pattern + "' " + "with? ")
                out_list = [replace if x == pattern else x for x in out_list]
                contin = input(
                    "Would you like to replace any more terms in the list? y/n"
                    )
            while contin == "n":
                return out_list
