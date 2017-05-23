# -*- coding: utf-8 -*-
"""tools for filtering DataFrames"""

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

import pandas as pd

class FilterTools(object):

    @classmethod
    def master_cleanse(cls, protein_df):
        """Filter raw protein DataFrame for master proteins.

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
        clean = None
        try:
            mask = protein_df.Master.str.endswith("IsMasterProtein")
            clean = protein_df.ix[mask]
        except Exception:
            pass
        return clean

# FilterTools.master_cleanse

    @classmethod
    def one_percent_expq(cls, protein_df):
        """Filters raw protein DataFrame for proteins that are less than 1% the
        expected q-value.

        Scans through the protein DataFrame selecting only the proteins with
        less than 1% of the experimental q-value.

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
        one_per = None
        try:
            mask = protein_df["Exp. q-value"] < .01
            one_per = protein_df.ix[mask]
            return one_per
        except KeyError:
            print("No 'Exp. q-value' column found.")
            return protein_df

    @classmethod
    def master_one(cls, protein_df):
        """Takes a raw protein DataFrame and filters it using first the
        'master_cleanse' function and 'one_percent_expq' function.

        Parameters
        ----------
        protein_df : DataFrame
            Raw proteins.

        Returns
        -------
        master_one : DataFrame
            Of master proteins with exp. q-value <1%
        """
        master = cls.master_cleanse(protein_df)
        master_one = cls.one_percent_expq(master)
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
    def vLook(cls, peptides=None, proteins=None):
        """Return a tuple of selected peptides and proteins.

        Takes raw peptides and protiens returns a tuple of selected peptides
        and proteins. The function can also select for a sigle modification or
        many modifications.

        Parameters
        ----------
        peptides : DataFrame
        proteins : DataFrame
        mods : list

        Returns
        -------
        peptide_select : DataFrame
        protein_select : DataFrame
        """
        fdr = cls.master_one(proteins)
        mpa = cls.masterPep(peptides)

        fdrdf = pd.DataFrame(fdr.Accession, index=fdr.index)

        peptide_select = mpa.merge(fdrdf, on="Accession",
                                   how="left", right_index=True)

        protein_select = mpa.merge(fdrdf, on="Accession",
                                   how="left", left_index=True)

        return peptide_select, protein_select
