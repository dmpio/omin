import pandas as pd
import numpy as np
# FIXME: Add to install reqs.
# FIXME: Add import try and except.
# FIXME: Test running w/o installed.
# FIXME: Extent to other intermine APIs.
from intermine.webservice import Service


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