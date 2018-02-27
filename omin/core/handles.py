# -*- coding: utf-8 -*-
# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

"""Omin core handles.

Handle in this context is a class composed of several pandas DataFrames, and
other varibles that are either derived from the DataFrames or provided by the
user.

"""

# ----------------
# EXTERNAL IMPORTS
# ----------------
import re
import os
import guipyter as gptr

# ----------------
# INTERNAL IMPORTS
# ----------------
# import the Handle super class.
from .base import Handle
from .containers import PeptideGroups, Proteins

# -------------
# UTILS IMPORTS
# -------------
from ..utils import IOTools
from ..utils import StringTools
from ..utils import SelectionTools
from ..utils import FilterTools

# -------------
# NORMALIZATION
# -------------
from ..normalize.toPool import NormalizedToPool
from ..normalize.toInput import NormalizedToInput

# --------
# DATBASES
# --------
from ..databases import mitoCartaCall
from ..utils.pandas_tools import pd



class RawData(Handle):
    """Handler for raw data.
    """

    def __init__(self, file_list=None, peptides_file=None, proteins_file=None):
        """Load data for peptides_file and proteins_file as pandas DataFrames.
        """
        super(RawData, self).__init__()
        if file_list is not None:
            rx = re.compile("[Pp]eptide")

            peptides_file = list(filter(rx.findall, file_list))[0]

            rx = re.compile("[Pp]roteins")
            proteins_file = list(filter(rx.findall, file_list))[0]

        self.peptide_groups = PeptideGroups(peptides_file)
        self.proteins = Proteins(proteins_file)



class PreProcess(RawData):
    """A metaclass that uses RawData attempting several filtering steps.

    DEPRECATED
    ----------

    Most of these steps will be carried out by:

        omin.core.containers.PeptideGroups

        omin.core.containers.Proteins

    The rest will be done in the Process class.

    """
    def __init__(self, file_list=None, peptides_file=None, proteins_file=None, modifications=None, *args, **kwargs):
        """Initialize instance of PreProcess class.
        """
        # Initalize the RawData base class.
        # super(PreProcess, self).__init__(file_list, peptides_file, proteins_file, *args, **kwargs)
        RawData.__init__(self, file_list, peptides_file, proteins_file, *args, **kwargs)


        # Find invivo modifications
        invivos = SelectionTools.findInVivoModifications(self.raw_peptides)
        self._invivo_modifications = invivos

        # Set modifications with given or derived.
        modifications = modifications or self._invivo_modifications

        # Create 2 Dataframes that map specific peptide or protein Uniprot ID to it's relevent mitocarta index.
        pep_sel, prot_sel = FilterTools.bridge(self.peptide_groups.raw, self.proteins.master_high_confidence)

        # Used to index filter to statistically relevant PeptideGroups
        self.pep_sel = pep_sel

        # Used to index filter to statistically relevant Proteins
        self.prot_sel = prot_sel

        # MitoCarta calls made
        mito, nonmito = mitoCartaCall.mitoCartaPepOut(self,
                                                      mods=modifications,
                                                      dex=True)
        # Mitocarta hits DataFrame
        self.mitodex = mito
        # Mitocarta non hits DataFrame
        self.nonmitodex = nonmito
        # Create a unified index of mitocarta calls.
        self.unidex = None
        try:
            unidex = pd.concat([self.mitodex, self.nonmitodex]).sort_index()
            self.unidex = unidex.reindex(index=self.peptide_groups.raw.index)
        except Exception:
            pass
        self.master_index = None
        try:
            # Get the Gene Symbol
            if "Gene ID" in self.proteins.master_high_confidence: # PD2.1
                entrez_gene_id = self.proteins.master_high_confidence["Gene ID"].iloc[self.prot_sel.index]
                entrez_gene_id.index = self.peptide_groups.raw.index
            if "Entrez Gene ID" in self.proteins.master_high_confidence: # PD2.2
                entrez_gene_id = self.proteins.master_high_confidence["Entrez Gene ID"].iloc[self.prot_sel.index]
                entrez_gene_id.index = self.peptide_groups.raw.index

            # Get Gene Description.
            ga = self.proteins.master_high_confidence["Description"].ix[self.prot_sel.index]
            ga.index = self.peptide_groups.raw.index
            mdex = pd.concat([self.unidex.ix[:, -4],
                              self.unidex.ix[:, -2:]],
                             axis=1)
            mdex = pd.DataFrame(mdex.fillna(0.0), dtype="bool")

            self.master_index = pd.concat([self.pep_sel,
                                           entrez_gene_id,
                                           ga,
                                           self.peptide_groups.raw.Modifications,
                                           self.peptide_groups.raw["Modifications in Proteins"],
                                           self.peptide_groups.raw.Sequence,
                                           mdex], axis=1)
            self.numbers['mitocarta_hits'] = self.master_index.MitoCarta2_List.sum()

        except Exception:
            print("Could not create master_index")
            pass
        # Find the number of inputs.
        self._input_number = SelectionTools.find_number_input(self.raw_peptides)
        # Find the number of PTMs present in the PeptideGroups
        self._ptm_fraction_numbers = SelectionTools.find_fractions(self.raw_peptides)


class Process(PreProcess):
    """A metaclass that uses PreProcess attempting several normalization steps.

    WARNING: This class is under construction switch to stable branch it you need to work.

    """

    # PROTIP: Wait util the last possible moment to link databases.
    # microprotip: Exceptions are the rule here.
    # PROTIP: TRY NOT TO SET VARS AT THIS LEVEL.

    def __init__(self, file_list=None, peptides_file=None, proteins_file=None, modifications=None, *args, **kwargs):
        """Initalize Process class.
        """
        # self.normalized = None
        super(Process, self).__init__(file_list, peptides_file, proteins_file, *args, **kwargs)
        # FIXME: Make the selection more specifically target abundance columns
        # FIXME: UPDATE this may be fixed now. Kill the vars set at ths level
        # FIXME: Remove all trace of the vile PreProcess class.
        # FIXME: Please dark lord of omics guide my hand as I destroy this class and rebuild it from scratch.

        # if self._input_number > 0:
        #     inp_notify = "{} input fraction(s) found. Normalizing now..."
        #     inp_notify = inp_notify.format(self._input_number)
        #     print(inp_notify)
        #     try:
        #         self.normalized = NormalizedToInput(self)
        #     except Exception:
        #         print("omin.normalize.toInput.NormalizedToInput FAILED.")
        #
        # elif self.raw_peptides.columns.str.contains("pool|control",
        #                                             case=False).any():
        #
        #     print("Pool columns. Omin will attempt to compare the data to it.")
        #     try:
        #         self.normalized = NormalizedToPool(self.raw_peptides)
        #     except Exception:
        #         print("omin.normalize.toPool.NormalizedToPool FAILED.")
        # else:
        #     # FIXME: Make this a place where the user could specify.
        #     print("Cannot find anything to normalize to.")
