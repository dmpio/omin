# -*- coding: utf-8 -*-
"""DEPRECATED: MOST FUNCTIONS HAVE MIGRATED TO THE pandomics PACKAGE."""

# DEPRECATE: EVERYTHING!

# import pandas as pd
#
#
# class FilterTools(object):
#     """Tools for filtering"""
#
#     @staticmethod
#     def is_master_protein(df):
#         "Filter raw protein DataFrame for master proteins."
#         try:
#             mask = df.Master == "IsMasterProtein"
#             result = df.loc[mask]
#             return result
#         except Exception as err:
#             print('Could not Filter for master proteins.', err)
#             return df
#
#     @staticmethod
#     def high_confidence(df):
#         """
#         Return a DataFrame of high confidence proteins.
#         """
#         result = None
#         try:
#             try: # PD 2.1
#                 mask = df["Exp. q-value"] < .01
#                 result =df.loc[mask]
#                 return result
#             except KeyError: # PD 2.2
#                 mask = df['Protein FDR Confidence: Combined'] == "High"
#                 result = df.loc[mask]
#                 return result
#         except Exception as err:
#             print("Could not filter for FDR")
#             return df
#
#     @staticmethod
#     def filterRow(dataframe_in=None, on=None, term=None, *args, **kwargs):
#         """Return DataFrame that contains a given term on specific column.
#         DEPRECATED
#         ----------
#         Parameters
#         ----------
#         dataframe : DataFrame
#         term : str
#         on : str
#
#         Returns
#         -------
#         filtered : DataFrame
#         """
#         filtered = None
#         if on in dataframe_in.columns:
#             try:
#                 mask = dataframe_in[on].str.contains(term, *args, **kwargs)
#                 filtered = dataframe_in[mask]
#                 return filtered
#             except Exception:
#                 print("omin.FilterTools.FilterRow FAILED")
#         else:
#             try:
#                 mask = dataframe_in.filter(regex=on).ix[:, 0]
#                 mask = mask.str.contains(term, *args, **kwargs).fillna(False)
#                 filtered = dataframe_in[mask]
#                 return filtered
#             except Exception:
#                 print("omin.FilterTools.FilterRow FAILED")
#
#     @classmethod
#     def master_cleanse(cls, protein_df):
#         """Filter raw protein DataFrame for master proteins.
#
#         DEPRECATED
#         ----------
#
#         The raw protein data from Proteome Discoverer there is a column with
#         the title 'Master' this funtion scans through that column and selects
#         only the proteins that end with the string "IsMasterProtein".
#
#         Parameters
#         ----------
#         protein_df : DataFrame
#             Raw protein DataFrame
#
#         Returns
#         -------
#         clean : DataFrame
#             DataFrame that contains only proteins with 'IsMasterProtein' in
#             'Master' column of protein_df.
#         """
#         clean = None
#         try:
#             mask = protein_df.Master.str.endswith("IsMasterProtein")
#             clean = protein_df.ix[mask]
#
#         except Exception:
#             pass
#         return clean
#
#
#     @classmethod
#     def one_percent_expq(cls, protein_df):
#         """Filters raw protein DataFrame for proteins that are less than 1% the
#         expected q-value.
#
#         Scans through the protein DataFrame selecting only the proteins with
#         less than 1% of the experimental q-value.
#
#         Parameters
#         ----------
#         protein_df : DataFrame
#             Raw protein DataFrame
#
#         Returns
#         -------
#         clean : DataFrame
#             Protein data that contains only proteins with proteins only less
#             than 1% of the expected q-value.
#         """
#         one_per = None
#         try:
#             mask = protein_df["Exp. q-value"] < .01
#             one_per = protein_df.ix[mask]
#             return one_per
#         except KeyError:
#             print("No 'Exp. q-value' column found.")
#             return protein_df

    # @classmethod
    # def master_one(cls, protein_df):
    #     """Takes a raw protein DataFrame and filters it using first the
    #     'master_cleanse' function and 'one_percent_expq' function.
    #
    #     Parameters
    #     ----------
    #     protein_df : DataFrame
    #         Raw proteins.
    #
    #     Returns
    #     -------
    #     master_one : DataFrame
    #         Of master proteins with exp. q-value <1%
    #     """
    #     master = cls.master_cleanse(protein_df)
    #     master_one = cls.one_percent_expq(master)
    #     return master_one
    #
    # @classmethod
    # def masterPep(cls, peptide_df):
    #     """Takes a peptide DataFrame and returns just the first master protein
    #     accession for each peptide.
    #
    #     Notes
    #     -----
    #     Assumes the first uniprot ID list is the correct one. Peptides with no
    #     master protein accession will be lost however the index of peptide_df will
    #     be preserved.
    #
    #     Parameters
    #     ----------
    #     peptide_df : DataFrame
    #
    #     Returns
    #     -------
    #     master_prot_acc : DataFrame
    #
    #     """
    #     master_prot_acc = [i.split(';')[0] for i in peptide_df['Master Protein Accessions'].dropna()]
    #
    #     master_prot_acc = pd.DataFrame(master_prot_acc,
    #                                    index=peptide_df['Master Protein Accessions'].dropna().index, columns=['Accession'])
    #     return master_prot_acc

    # @classmethod
    # def vLook(cls, peptides=None, proteins=None):
    #     """Return a tuple of selected peptides and proteins.
    #
    #     Takes raw peptides and protiens returns a tuple of selected peptides
    #     and proteins. The function can also select for a sigle modification or
    #     many modifications.
    #
    #     Parameters
    #     ----------
    #     peptides : DataFrame
    #     proteins : DataFrame
    #     mods : list
    #
    #     Returns
    #     -------
    #     peptide_select : DataFrame
    #     protein_select : DataFrame
    #     """
    #     fdr = cls.master_one(proteins)
    #     mpa = cls.masterPep(peptides)
    #
    #     fdrdf = pd.DataFrame(fdr.Accession, index=fdr.index)
    #
    #     peptide_select = mpa.merge(fdrdf, on="Accession",
    #                                how="left", right_index=True)
    #
    #     protein_select = mpa.merge(fdrdf, on="Accession",
    #                                how="left", left_index=True)
    #
    #     return peptide_select, protein_select

    # @classmethod
    # def first_mpa(cls, peptide_df):
    #     """Return the first master protein.
    #
    #     Notes
    #     -----
    #     Assumes the first uniprot ID list is the correct one.
    #
    #     Parameters
    #     ----------
    #     peptide_df : DataFrame
    #
    #     Returns
    #     -------
    #     master_prot_acc : DataFrame
    #
    #     """
    #     mpa = peptide_df['Master Protein Accessions'].dropna()
    #     master_prot_acc = [i.split(';')[0] for i in mpa]
    #
    #     master_prot_acc = pd.DataFrame(master_prot_acc,
    #                                    index=mpa.index,
    #                                    columns=['Accession'])
    #     master_prot_acc = master_prot_acc.reindex(peptide_df.index)
    #     return master_prot_acc
    #
    # @staticmethod
    # def first_filter(dataframe, on, column_name=None):
    #     """Filter the filter element of a string.
    #
    #     FIXME: FIX THIS FUNCTION
    #     """
    #     # TODO: Deprecate in favor of omin.utils.pandas_tools.first_member
    #     ser = dataframe[on].dropna()
    #     # first = [i.split(';')[0] for i in ser]
    #     first = ser.apply(lambda x: x.split(";")[0])
    #     if column_name is not None:
    #         first = pd.DataFrame(first, columns=column_name)
    #     else:
    #         first = pd.DataFrame(first, columns=[on])
    #     first = first.reindex(index=dataframe.index)
    #     return first

    # @classmethod
    # def bridge(cls, peptides=None, proteins=None):
    #     """Return a tuple of selected peptides and proteins.
    #
    #     Parameters
    #     ----------
    #     peptides : DataFrame
    #     proteins : DataFrame
    #     mods : list
    #
    #     Returns
    #     -------
    #     peptide_select : DataFrame
    #     protein_select : DataFrame
    #     """
    #     mpa = cls.first_mpa(peptides)
    #     prot_acc = pd.DataFrame(proteins.Accession, index=proteins.index)
    #
    #     peptide_select = mpa.merge(prot_acc, on="Accession",
    #                                how="left", right_index=True)
    #
    #     protein_select = mpa.merge(prot_acc, on="Accession",
    #                                how="left", left_index=True)
    #
    #     return peptide_select, protein_select
