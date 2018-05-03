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

# Ugly hack to find this module's version number.
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
    """A metaclass that uses PreProcess attempting several normalization steps.

    WARNING: This class is under construction switch to stable branch if you need to work.
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

        # Link proteins to peptides
        self._link_proteins_to_peptides()

        # # Attempt to calculate the relative occupancy.
        self._calculate_relative_occupancy(verbose=verbose)

        # Reset the master_index
        self._reset_master_index(verbose=verbose)


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


    def _link_proteins_to_peptides(self):
        """Links protiens to peptides.
        """
        link_to_peptides = self.peptide_groups.master_index.merge(self.proteins.master_index, on="Accession", how="left", left_index=True)
        link_to_peptides = link_to_peptides.Accession
        self.proteins.link_to_peptides = link_to_peptides


    def _calculate_relative_occupancy(self, verbose=False):
        """Calculate the relative occupancy if possible.
        """
        self.peptide_groups.relative_occupancy = Occupancy()
        self.proteins.load_normalized = Normalized()
        self.proteins.relative_occupancy = Occupancy()

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
            related_proteins_dict = dict()
            for k,v in self.proteins.load_normalized.__dict__.items():
                related_proteins = v.iloc[self.proteins.link_to_peptides.index]
                related_proteins.index = self.peptide_groups.raw.index
                related_proteins_dict[k] = related_proteins
            self.peptide_groups.load_normalized_related_proteins = Normalized(**related_proteins_dict)
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


    def _reset_master_index(self, verbose=True):
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
            The term to filterfor the denominator.

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
