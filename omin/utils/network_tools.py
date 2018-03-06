# -*- coding: utf-8 -*-
"""omin.intermine_tools

Provides
--------
Tools for querying the intermine database.
Tools for investigating with fasta and uniprot.
"""
# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

# ==========
# TO DO LIST
# ==========
# FIXME: Intermine rescue function fails from wierd indexing of the master_index.
# FIXME: Extent to other intermine APIs.

# ----------------
# EXTERNAL IMPORTS
# ----------------
import re
import pandas as pd
import numpy as np
from urllib.error import HTTPError
from urllib.request import urlopen
# Try to import the intermine package.
try:
    from intermine.webservice import Service
except ImportError as err:
    print("Cannot import Intermine:", err)


# ==============
# INTERMINETOOLS
# ==============

class IntermineTools(object):

    @staticmethod
    def mousemine_accession_lookup(accession, verbose=False):
        """Return an intermine query object for a given protien accession number.

        The query is design to find Entrez gene IDs for proetiens of a given accession
        number. It only considers the cannonical isoform of whole mouse proteins(not fragments)
        presnt in the Swiss-Prot or TrEMBL databases.

        Parameters
        ----------
        accession: str
            Accession numbers in the form: XXXXXX-N, will be shortented to: XXXXXX.

        verbose: bool
            Defaults to False set to True for printed error messages.

        Returns
        -------
        result: (:obj)
            Intermine query object.

        See Also
        --------
        mousemine_query_format
        mousemine_accession_lookup_reduce
        mousemine_accession_to_entrez

        """
        query = None
        # Remove the isoform number if present.
        if "-" in accession:
            accession = accession.split('-')[0]
        try:
            # Begin intermine call.
            service = Service("http://www.mousemine.org/mousemine/service")
            query = service.new_query("Gene")

            query.add_view("primaryIdentifier",
                           "ncbiGeneNumber",
                           "proteins.primaryAccession",
                           "symbol",
                           "proteins.synonyms.value",
                           "proteins.length",
                           "proteins.isFragment",
                           "proteins.dataSets.name")

            # Declare sort order of results.
            query.add_sort_order("Gene.symbol", "ASC")
            query.add_sort_order("Gene.proteins.dataSets.name", "ASC")
            query.add_sort_order("Gene.proteins.length", "DESC")

            # Declare constraints.
            query.add_constraint("organism.taxonId", "=", "10090", code="B")
            query.add_constraint("proteins.isFragment", "=", "False", code="C")
            query.add_constraint("proteins.dataSets.name",
                                 "ONE OF",
                                 ["Swiss-Prot data set", "TrEMBL data set"],
                                 code="D")

            # Main constraint: find all with accession number X.
            query.add_constraint("Gene", "LOOKUP", accession, code="A")
            return query

        except Exception as err:
            if verbose:
                print(err)
            return query

    @staticmethod
    def mousemine_query_format(query, format_list=None, verbose=False):
        """Return a DataFrame from an intermine query object.

        Parameters
        ----------
        query: (:obj)
            Intermine query object.

        format_list: list
            List of information to include in the result.

        Returns
        -------
        result: pandas.DataFrame
            With columns labels based on format list.
        """
        result = None
        if format_list is None:
            format_list = ["primaryIdentifier", "ncbiGeneNumber", "proteins.primaryAccession",
                           "symbol", "proteins.synonyms.value", "proteins.isFragment"]
        # Remove the dots from the terms in format for DataFrame safety.
        column_labels = list(map(lambda x: x.replace('.', ''), format_list))

        result = list()
        try:
            for row in query.rows():
                row_out = [row[j] for j in format_list]
                result.append(row_out)

            result = pd.DataFrame(result)
            result.columns = column_labels
            return result

        except Exception as err:
            if verbose:
                print(err)
            return result

    @staticmethod
    def mousemine_accession_lookup_reduce(df, verbose=False):
        """Returns Entrez Gene ID as string from intermine query as pandas.DataFrame.

        Will return a string if the given query df has a single primary identifier
        and only one Entrez Gene ID. In cases where there are multiple primary
        identifiers and/or Entrez Gene IDs a np.NaN is returned instead.

        Parameters
        ----------
        df: pandas.DataFrame

        Returns
        -------
        result: str or np.nan
            Entrez Gene ID

        """
        result = np.nan
        if df is not None:
            try:
                if len(df.primaryIdentifier.unique()) == 1:
                    if len(df.ncbiGeneNumber.unique()) == 1:
                        result = df.ncbiGeneNumber.unique()[0]
                        return result
                else:
                    return result

            except Exception as err:
                if verbose:
                    print(err)
                return result

    @classmethod
    def mousemine_accession_to_entrez(cls, accession, verbose=False):
        """Return an Entrez Gene ID as a string for a given proteins accession number.
        Parameters
        ----------
        accession: str

        Returns
        -------
        entrez_gene_id: str or np.nan
            Entrez Gene ID

        """
        entrez_gene_id = np.nan
        try:
            query = cls.mousemine_accession_lookup(accession, verbose=verbose)
            df = cls.mousemine_query_format(query, verbose=verbose)
            entrez_gene_id = cls.mousemine_accession_lookup_reduce(df, verbose=verbose)
            return entrez_gene_id

        except Exception as err:
            if verbose:
                print(err)
            return entrez_gene_id

    @classmethod
    def rescue_entrez_ids(cls, protein_master_index):
        """Attempts to fill in missing Entrez Gene IDs from Intermine.
        """
        print("Fetching missing Entrez Gene IDs from intermine. This may take a few minutes...")
        c=0
        for i in protein_master_index.index:
            if protein_master_index.loc[i].EntrezGeneID is np.nan:
                entrez_id = cls.mousemine_accession_to_entrez(protein_master_index.loc[i].Accession)
                if entrez_id is not np.nan:
                    c+=1
                    protein_master_index.loc[i].EntrezGeneID = entrez_id

        print(c,"Retrieved Entrez Gene IDs from availible protein accession numbers.")
        print(protein_master_index.EntrezGeneID.isnull().sum(), "Entrez Gene IDs could not be determined.")


# ==========
# FASTATOOLS
# ==========

# FIXME: find the relevant files in method development and rewrite this class from scratch.

class FastaTools(object):
    """A collection of tools for handling FASTA formatted strings."""

    @classmethod
    def fasta2Seq(cls, fasta):
        """Join strings in list into one string.

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
        """Take a sequence and returns the number for each amino acid.

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

# ============
# UNIPROTTOOLS
# ============

class UniProtTools(object):
    """Tools for querying the UniProt database."""

    id_rx = re.compile('ID;\s(\d+);')

    @classmethod
    def get_uniprot(cls, accession, file_format='.fasta'):
        """Retrieve UniProt data for a given accession."""
        if file_format[0] != '.':
            file_format = '.' + file_format

        url = "http://www.uniprot.org/uniprot/"
        req = ''.join([url, accession, file_format])
        try:
            data = urlopen(req)
            return data.read().decode()
        except HTTPError:
            print('Entry :', accession, 'not found.')
            return None

    @classmethod
    def get_gene_id(cls, accession):
        """Retrieve the gene id for a given accession."""
        result = cls.get_uniprot(accession, file_format='.txt')
        result = cls.id_rx.findall(result)
        result = result[0]
        return result

    @classmethod
    def go_anno(cls, accession, desc=False):
        """Retrieve GO annotations for a given accession."""
        txt = cls.get_uniprot(accession, file_format='.txt')
        if desc:
            rx = re.compile('GO:\d+;.+\n')
        else:
            rx = re.compile('(GO:\d+);')
        results = rx.findall(txt)
        return results
