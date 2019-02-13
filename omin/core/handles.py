# -*- coding: utf-8 -*-
"""
omin.core.handles
-----------------

Provides handles for omin. A handle in this context is a class composed of
several pandas DataFrames, and other varibles that are either derived from the
DataFrames or provided by the user.
"""
# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

# ----------------
# EXTERNAL IMPORTS
# ----------------
import re
import os
import numpy as np
from pandomics import pandas as pd
# Get the version numbers of the hard dependencies.
from guipyter import __version__ as guipyter_version
from pandomics import __version__ as pandomics_version

# ----------------
# INTERNAL IMPORTS
# ----------------
from .base import repr_dec
from .containers import PeptideGroups, Proteins, Occupancy, Normalized

from ..databases import mitocarta


# Ugly hack to find this module's version number.
# FIXME: use bump for versioning instead.
__module_path__ = os.path.dirname(os.path.realpath(__file__))
__module_path__ = os.path.split(__module_path__)[0]

with open(os.path.join(__module_path__, '__version__')) as f:
    omin_version = f.read().strip()


# =============
# PROJECT CLASS
# =============
@repr_dec
class Project(object):
    """
    """
    def __init__(self, file_list=None, peptides_file=None, proteins_file=None, rescue_entrez_ids=None, verbose=None, *args, **kwargs):
        """Load data for peptides_file and proteins_file

        Parameters
        ----------
        file_list: str
            A list of files

        peptides_file: str or _io.TextWrapper

        proteins_file: str or _io.TextWrapper

       rescue_entrez_ids: bool
            Attempts to query the Intermine database for proteins that have
            Master Protein Accessions but no Entrez Gene ID. This process can
            take several minutes. Defaults to False.

        verbose: bool
            Defaults to True.
        """
        # If verbose is not declared then set it to True.
        if verbose is None:
            verbose = True

        # Set all of the version info for object verification.
        self._version_info = dict()
        self._version_info["omin"] = omin_version
        self._version_info["guipyter"] = guipyter_version
        self._version_info["pandomics"] = pandomics_version

        # FIXME: For some reason even if a list is provided guipyter is stil triggered.
        if file_list is not None:
            rx = re.compile("[Pp]eptide")
            peptides_file = list(filter(rx.findall, file_list))[0]

            rx = re.compile("[Pp]roteins")
            proteins_file = list(filter(rx.findall, file_list))[0]

        # FIXME: Put in kwargs switch like in guipyter.
        # Load the Peptide Groups file.
        if peptides_file is not None:
            try:
                self.peptide_groups = PeptideGroups(filepath_or_buffer=peptides_file)
            except Exception as err:
                if verbose:
                    print(err)
        else:
            self.peptide_groups = PeptideGroups(*args, **kwargs)

        # Load the Proteins file.
        if proteins_file is not None:
            try:
                self.proteins = Proteins(filepath_or_buffer=proteins_file, rescue_entrez_ids=rescue_entrez_ids)
            except Exception as err:
                if verbose:
                    print(err)
        else:
            self.proteins = Proteins(rescue_entrez_ids=rescue_entrez_ids, *args, **kwargs)

# =============
# PROCESS CLASS
# =============

class Process(Project):
    """Attempts several normalization steps.
    """

    # PROTIP: Wait util the last possible moment to link databases.
    # microprotip: Exceptions are the rule here.
    # PROTIP: TRY NOT TO SET VARS AT THIS LEVEL.

    def __init__(self, *args, **kwargs):
        """Initalize Process class.
        """
        # Initialize the Project class.
        if "verbose" in kwargs:
            verbose = kwargs['verbose']
        else:
            verbose = False

        Project.__init__(self, *args, **kwargs)
        # Connect master index from peptide groups to proteins.
        self._peptide_groups_master_index_update()

        # Fill the NaNs in mitocarta.
        self._peptide_groups_mitocarta_fillna()

        # Reset the master_index.
        self._reset_proteins_master_index(verbose=verbose)

        # Link proteins to peptides
        self._link_proteins_to_peptides()

        # Reset the peptides groups master index
        self._reset_peptide_groups_master_index()

        # Remove the extra gene id column that is introduced from the function above.
        self._peptide_groups_entrez_gene_clean_up()

        # # Attempt to calculate the relative occupancy.
        self._calculate_relative_occupancy(verbose=verbose)

        # Add mitocarta info to metadata.
        self._add_mitocarta_metadata(verbose=verbose)

        # Extract the gene names.
        self.peptide_groups._gene_name_extractor()

        self.proteins._gene_name_extractor()

    def _peptide_groups_master_index_update(self):
        """Merge the proteins.master_index with the peptide_groups.master_index.
        """
        try:
            updated_master_index = self.peptide_groups.master_index.merge(self.proteins.master_index, on="Accession", how="left")
            self.peptide_groups.master_index = updated_master_index

        except Exception as err:
            if verbose:
                print(err)


    def _peptide_groups_mitocarta_fillna(self):
        """Attempts to fill missing values (NaNs) created by merging the proteins.master_index with the peptide_groups.master_index.
        """
        if "MitoCarta2_List" in self.peptide_groups.master_index:
            try:
                # ASSUMPTION: If MitoCarta2_List is present then IMS and Matrix will be aswell.
                self.peptide_groups.master_index["MitoCarta2_List"].fillna(False, inplace=True)
                self.peptide_groups.master_index["IMS"].fillna(False, inplace=True)
                self.peptide_groups.master_index["Matrix"].fillna(False, inplace=True)
            except Exception as err:
                if verbose:
                    print(err)


    def _reset_proteins_master_index(self, verbose=True):
        """Work-around to return rows with missing values to protein.master_index.
        """
        # Create a copy of the old_master_index.
        old_master_index = self.proteins._old_master_index.copy()

        # Retrieve the accession from the old_master_index.
        old_master_index = old_master_index.filter(regex="Accession")
        # Merge the old_master_index with the "new" master_index
        result = old_master_index.merge(self.proteins.master_index, how="left", on="Accession")
        # Set the index.
        result.index = self.proteins._old_master_index.index
        # Fill the missing values in the Mitocarta columns with False.
        # FIXME: Generalize this so that it isn't so mitocarta specific
        result["MitoCarta2_List"].fillna(False, inplace=True)
        result["Matrix"].fillna(False, inplace=True)
        result["IMS"].fillna(False, inplace=True)
        self.proteins.master_index = result


    def _reset_peptide_groups_master_index(self):
        """Work-around that fill missing descriptions in the peptide groups master index.
        """
        ind = self.proteins.raw.copy()[["Accession", "EntrezGeneID", "Description"]]

        ind.EntrezGeneID = ind.EntrezGeneID.astype(str)

        mito = mitocarta.MitoCartaTwo.essential.copy()

        mito.EntrezGeneID = mito.EntrezGeneID.astype(str)

        ind = ind.merge(mito, on="EntrezGeneID", how="left")

        trun = self.peptide_groups.master_index.copy().iloc[:, :5]

        result = trun.merge(ind, on="Accession", how="left")

        result[["MitoCarta2_List", "Matrix", "IMS"]] = result[["MitoCarta2_List", "Matrix", "IMS"]].fillna(False)

        self.peptide_groups.master_index = result

    def _peptide_groups_entrez_gene_clean_up(self):
        """Workaround that cleans up the extra EntrezGeneID column introduced by the function above.
        """
        if "EntrezGeneID_y" in self.peptide_groups.master_index:
            self.peptide_groups.master_index.drop("EntrezGeneID_y", axis=1, inplace=True)
        if "EntrezGeneID_x" in self.peptide_groups.master_index:
            self.peptide_groups.master_index.rename({"EntrezGeneID_x":"EntrezGeneID"}, axis=1, inplace=True)

    def _link_proteins_to_peptides(self):
        """Links protiens to peptides.

        Creates a series that serves as a map between proteins and peptides
        groups. The series has the same number of rows as the peptides groups
        but the index corresponds to the proteins.
        """

        link_to_peptides = self.peptide_groups.master_index.merge(self.proteins.master_index,
                                                                  on="Accession",
                                                                  how="left",
                                                                  left_index=True,
                                                                  validate="m:1")
        link_to_peptides = link_to_peptides.Accession
        self.proteins.link_to_peptides = link_to_peptides


    def _normalize_input_protiens_to_input_peptides(self, verbose=False):
        """Normalize input proteins to input peptides.
        """
        self.proteins.load_normalized = Normalized()

        if self.proteins.input_number > 0:
            if verbose:
                print("Input fractions found calculating relative occupancy...")
            input_mask = self.proteins.study_factor_table[self.proteins.study_factor_with_input].str.contains("[Ii]nput")
            number_input_fractions = len(self.proteins.study_factor_table.loc[input_mask]._Fn.unique())
            # isolate the input fractions study factors.
            inps = self.proteins.study_factor_table.loc[input_mask]
            inps = [inps.loc[inps._Fn.str.contains(i)] for i in inps._Fn.unique()]
            normalized = dict()
            for inp in inps:
                if len(inp._Fn.unique()) == 1:
                    inp_fn = inp._Fn.unique()[0]
                    inp_tag = self.proteins.fraction_tag(inp_fn)

                    if inp_tag in self.peptide_groups.load_normalized.__dict__:
                        inp_pr_linked = self.proteins.Abundance[self.proteins.Abundance.columns[inp.index]]
                        inp_pg_linked = self.peptide_groups.Abundance[self.peptide_groups.Abundance.columns[inp.index]]
                        load_norm = inp_pr_linked.normalize_to(inp_pg_linked)
                        normalized[inp_tag] = load_norm
            self.proteins.load_normalized = Normalized(**normalized)


    def _link_related_load_normalized_protein_abundances_to_peptides(self, verbose=False):
        """
        """
        # Check for the presence of input fractions.
        if self.proteins.input_number > 0:
            # Normalize input proteins to input peptides.
            self._normalize_input_protiens_to_input_peptides()
            related_proteins_dict = dict()
            for k,v in self.proteins.load_normalized.__dict__.items():
                # related_proteins = v.iloc[self.proteins.link_to_peptides.index]
                related_proteins = v.loc[self.proteins.link_to_peptides.index]
                # related_proteins.index = self.peptide_groups.raw.index
                related_proteins.index = self.peptide_groups.master_index.index
                related_proteins_dict[k] = related_proteins
            self.peptide_groups.load_normalized_related_proteins = Normalized(**related_proteins_dict)


    def _filter_out_false_hit_related_proteins(self):
        """Removes the false hits in the related proteins.

        For some unknown reason false hits are generated from merging the
        peptides and proteins data.
        """
        accessions_in_proteins = self.peptide_groups.master_index.Accession.isin(self.proteins.master_index.Accession)

        for k,v in self.peptide_groups.load_normalized_related_proteins.__dict__.items():
            related_proteins_filtered = v.loc[accessions_in_proteins]
            related_proteins_filtered = related_proteins_filtered.reindex(v.index)
            self.peptide_groups.load_normalized_related_proteins.__dict__[k] = related_proteins_filtered


    def _calculate_relative_occupancy(self, verbose=False):
        """Calculate the relative occupancy if possible.
        """
        # Initalize the empty relative occupancy container.
        self.peptide_groups.relative_occupancy = Occupancy()
        # Check for the presence of input fractions.
        if self.proteins.input_number > 0:
            self._link_related_load_normalized_protein_abundances_to_peptides()
            self._filter_out_false_hit_related_proteins()

            # Calculate relative occupancy.
            occupancy = dict()
            for k,v in self.peptide_groups._linked_fractions.items():
                peptide_norm = self.peptide_groups.load_normalized.__dict__[k].log2_normalize()
                proteins_related_norm = self.peptide_groups.load_normalized_related_proteins.__dict__[v].log2_normalize()
                result = peptide_norm.subtract_by_matrix(proteins_related_norm, prepend_cols="Relative Occupancy: ")
                occupancy[k] = result
            self.peptide_groups.relative_occupancy = Occupancy(**occupancy)

        # No input fractions found.
        else:
            if verbose:
                print("Could not generate relative occupancy from data.")
            pass


    def _add_mitocarta_metadata(self, verbose=False):
        for mod in self.peptide_groups._in_vivo_modifications:
            # Filter for given in vivo modification.
            try:
                result = self.peptide_groups.master_index.filter_rows(on="Modifications", term=mod)

                self.peptide_groups.metadata[mod.lower() + "_mitocarta_total_hits"] = result.MitoCarta2_List.sum()
                self.peptide_groups.metadata[mod.lower() + "_mitocarta_matrix"] = result.Matrix.sum()
                self.peptide_groups.metadata[mod.lower() + "_mitocarta_ims"] = result.IMS.sum()

            except Exception as err:
                if verbose:
                    print("Could not add Mitocarta info to metadata.", err)


    def comparision(self, on=None, where=None, fraction_key=None, mask=None, right=None,
                    numerator=None, denominator=None, filter_out_numerator=None, filter_out_denominator=None):
        """Creates a comparision dataframe from a Process object.

        Parameters
        ----------
        on: str
            either: peptide_groups or proteins.

        where: str
            either: load_normalized or relative_occupancy.

        fraction_key: str
            Can be whatever fraction keys that are availible in the process object.

        mask: pd.Index
            Uses this index to select for.

        right: DataFrame
            When comparing to another separate DataFrame taken as the denominator.
            Defaults to None.

        numerator: str
            The term to filter for the numerator.

        denominator: str
            The term to filter for the denominator.

        filter_out_numerator: bool
            Defaults to False.

        filter_out_denominator: bool
            Defaults to False.

        Returns
        ------
        result: pd.DataFrame
        """
        # FIXME: Add ANDing of masks
        # FIXME: ADD THIS FUNCTION TO THE NORMALIZED and/or OCCUPANCY classes.
        # FIXME: Consider adding more kwargs like fraction_key for fine tuning.
        result=None

        if mask is not None:
            # Isolate just the masked peptides.
            try:
                shown = self.__getattribute__(on).__getattribute__(where).__dict__[fraction_key].iloc[mask]

            except Exception as err:
                print(err)
                print("Cannot use mask defaulting to whole data set.")
                shown = self.__getattribute__(on).__getattribute__(where).__dict__[fraction_key]

        else:
            shown = self.__getattribute__(on).__getattribute__(where).__dict__[fraction_key]

        ## Drop the rows with only missing values in the Abundance columns.
        abun = shown.dropna(how="all")

        ## Log2 Normalize the load normalized.
        if where == "relative_occupancy":
            log2load = abun
        else:
            log2load = abun.log2_normalize(prepend_cols="Log2 Normalized: ")

        ## Create the comparison.
        ## FIXME: Add kwargs switch for right ect.
        comp = log2load.fold_change_with_ttest(right=right,
                                               numerator=numerator,
                                               denominator=denominator,
                                               filter_out_numerator=filter_out_numerator,
                                               filter_out_denominator=filter_out_denominator)

        ## Take the -log10(p-value)
        comp["negative_log10_pvalue"] = -np.log10(comp.pvalue)

        #truncated_index = self.__getattribute__(on).master_index

        ## Grab the related annotations.
        truncated_index = self.__getattribute__(on).master_index.loc[comp.index]

        ## Package for export.
        if where == "load_normalized":
            result = pd.concat([truncated_index, abun, log2load, comp], axis=1)

        if where == "relative_occupancy":
            result = pd.concat([truncated_index, log2load, comp], axis=1)

        return result
