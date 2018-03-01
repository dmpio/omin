# -*- coding: utf-8 -*-
"""
Omin core handles.

Handle in this context is a class composed of several pandas DataFrames, and
other varibles that are either derived from the DataFrames or provided by the
user.
"""
# -------
# LiCENSE
# -------
# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

# ----------
# TO DO LIST
# ----------
# FIXME: DOCUMENT OR DIE #DOD
# FIXME: Add relative occupancy method to Process

# ----------------
# EXTERNAL IMPORTS
# ----------------
import re

# ----------------
# INTERNAL IMPORTS
# ----------------
from .containers import PeptideGroups, Proteins, Occupancy, Normalized

# =============
# PROJECT CLASS
# =============

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
        #rescue_entrez_ids =rescue_entrez_ids or False
        verbose = verbose or True
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
        self.peptide_groups_master_index_update()
        self.peptide_groups_mitocart_fillna()
        # # Attempt to calculate the relative occupancy
        self.calculate_relative_occupancy(verbose=verbose)

    def peptide_groups_master_index_update(self):
        """Merge the proteins.master_index with the peptide_groups.master_index.
        """
        try:
            updated_master_index = self.peptide_groups.master_index.merge(self.proteins.master_index, on="Accession", how="left")
            self.peptide_groups.master_index = updated_master_index

        except Exception as err:
            if verbose:
                print(err)


    def peptide_groups_mitocart_fillna(self):
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


    def calculate_relative_occupancy(self, verbose=False):
        """Calculate the relative occupancy is possible.
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
                    # print(inp_tag)
                    if inp_tag in self.peptide_groups.load_normalized.__dict__:

                        inp_pr_linked = self.proteins.Abundance[self.proteins.Abundance.columns[inp.index]]

                        inp_pg_linked = self.peptide_groups.Abundance[self.peptide_groups.Abundance.columns[inp.index]]

                        load_norm = inp_pr_linked.normalize_to(inp_pg_linked)
                        normalized[inp_tag] = load_norm

            self.proteins.load_normalized = Normalized(**normalized)
        # No input fractions found.
        else:
            pass
