# -*- coding: utf-8 -*-
# Copyright 2017 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach,
# Blair Chesnut, and Elizabeth Hauser.

"""Omin core handles.

Handle in this context is a class composed of several pandas DataFrames, and
other varibles that are either derived from the DataFrames or provided by the
user.

"""

# FIXME: Investigate a SQLite/Json file stradegy.
# FIXME: Store each handle class as SQLite database in same parent dir.
# FIXME: Define varibles used at the lowest possible class level.
# FIXME: Add type checking.
# FIXME: Add more try and excepts but try to put them on function level
# FIXME: Include paragraph description of the types of filtering.

import re
import gc
# import pandas as pd
import numpy as np
from omin.utils import StringTools
from omin.utils import SelectionTools
from omin.utils import FilterTools
from omin.normalize.toPool import NormalizedToPool
from omin.normalize.toInput import NormalizedToInput
from omin.databases import mitoCartaCall
from omin.utils.pandas_tools import pd

# Hotpatched functions for pd.DataFrame
# -------------------------------------
# Hotpatching aka monkey-patching functions to an existing class is suggested
# over subclassing functions. I prefer to use this method instead of pandas
# pd.DataFrame.pipe command because I think that it yields a simpler syntax for
# for the lab rats.

#
# def log2(self):
#     """Take the Log2 of all values in a DataFrame."""
#     return self.apply(np.log2)
#
#
# setattr(pd.DataFrame, "log2", log2)
#
#
# def log2_normalize(self):
#     """Return a Log2 normalized DataFrame Log2(values)-mean(Log2(values))."""
#     return self.log2().sub(self.log2().mean(axis=1), axis=0)
#
#
# setattr(pd.DataFrame, "log2_normalize", log2_normalize)
#
#
# def compare_to(self, *args, **kwargs):
#     """Subtract all columns of a DataFrame by specified column."""
#     result = self.subtract(self.filter(*args, **kwargs).iloc[:, 0], axis=0)
#     return result
#
#
# setattr(pd.DataFrame, "compare_to", compare_to)
#
#
# def normalization_factors(self):
#     """Compute normalization factors for a DataFrame."""
#     norm_factors = self.sum() / self.sum().mean()
#     return norm_factors
#
#
# setattr(pd.DataFrame, "normalization_factors", normalization_factors)
#
#
# def normalize_to(self, normal):
#     """Normalize a DataFrame to another DataFrame."""
#     # Sort by columns.
#     self_sort = self.columns.tolist()
#     self_sort.sort()
#     self = self[self_sort]
#     # Sort the normal dataframe by it's columns.
#     normal_sort = normal.columns.tolist()
#     normal_sort.sort()
#     normal = normal[normal_sort]
#     # Divide by the normalization factors.
#     normalized = self / normal.normalization_factors().as_matrix()
#     normalized.columns = self.columns + ": Normalized to: " + normal.columns
#     return normalized
#
#
# setattr(pd.DataFrame, "normalize_to", normalize_to)


# Omin's core handle objects
# ---------------------------
# These are essentially containers for DataFrames.

class Handle(object):
    """The core omin handle base class."""

    def __init__(self):
        """Initalize the core handle."""
        self.numbers = dict()

    def __repr__(self):
        """Show all attributes."""
        return "Attributes: "+", ".join(list(self.__dict__.keys()))


class ProteomeDiscovererRaw(Handle):
    """Base class for Proteome Discoverer raw files."""

    def __init__(self, raw_data):
        """Initalize base class for Proteome Discoverer raw files."""
        super(ProteomeDiscovererRaw, self).__init__()
        self.raw = raw_data.copy()
        self.abundance = self.raw.filter(regex="Abundance:")
        self.numbers['total_ids'] = self.raw.shape[0]


class PeptideGroups(ProteomeDiscovererRaw):
    """Base class for Peptide Groups."""

    def __init__(self, *args, **kwargs):
        """Initialize the base class."""
        super(PeptideGroups, self).__init__(*args, **kwargs)
        # Find invivo modifications.
        in_vivo_mods = SelectionTools.findInVivoModifications(self.raw)
        # Declare varible
        self._in_vivo_modifications = []
        # Check if None.
        if in_vivo_mods is not None:
            # If the list is greater than zero then set varible.
            if len(in_vivo_mods) > 0:
                self._in_vivo_modifications = in_vivo_mods
                # For each modification create a filtered dataframe attribute.
                for i in in_vivo_mods:
                    mod = StringTools.remove_punctuation(i).lower()
                    df = SelectionTools.filterRow(self.raw,
                                                  on="Modifications",
                                                  term=i)
                    self.numbers[mod] = df.shape[0]
                    self.__dict__[mod] = df
            else:
                pass
        else:
            pass


class Proteins(ProteomeDiscovererRaw):
    """Base class for Proteins."""

    def __init__(self, *args, **kwargs):
        """Initalize the base class."""
        super(Proteins, self).__init__(*args, **kwargs)


class RawData(Handle):
    """Converts Proteome Discoverer .txt files into pandas DataFrames.

    Attributes
    ----------
    peptides :  DataFrame
        Raw data from Proteome Discoverer peptides data.
    proteins : DataFrame
        Raw data Proteome Discoverer corresponding proteins data.
    _numbers : tuple
    """

    def __init__(self, file_list=None, peptides_file=None, proteins_file=None):
        """Load data for peptides_file and proteins_file as pandas DataFrames.

        Parameters
        ----------
        file_list : list
            A list of files.
        peptides_file : str
            Name of Proteome Discoverer Peptides file as string.
        proteins_file : str
            Name of Proteome Discoverer Proteins file as string.

        Examples
        --------
        Load file strings:
        >>>peptides_file = "mydatafolder/peptides.txt"
        >>>proteins_file = "mydatafolder/proteins.txt"

        Create RawData object:
        >>>raw_data = RawData(peptides_file,proteins_file)

        """
        super(RawData, self).__init__()
        if file_list is not None:
            rx = re.compile("[Pp]eptide")

            peptides_file = list(filter(rx.findall, file_list))[0]

            rx = re.compile("[Pp]roteins")
            proteins_file = list(filter(rx.findall, file_list))[0]

        # Load your peptide groups file as a pandas DataFrame.
        self.raw_peptides = pd.read_csv(peptides_file,
                                        delimiter="\t",
                                        low_memory=False)
        # Load your protein file as a pandas DataFrame.
        self.raw_proteins = pd.read_csv(proteins_file,
                                        delimiter="\t",
                                        low_memory=False)
        # Load the RawData into their respective classes.
        self.peptide_groups = PeptideGroups(raw_data=self.raw_peptides)
        self.proteins = Proteins(raw_data=self.raw_proteins)
        # Store the shape of the respective DataFrames.
        self._numbers = (self.raw_peptides.shape, self.raw_proteins.shape)

    def __repr__(self):
        """Show all attributes."""
        return "Attributes: "+", ".join(list(self.__dict__.keys()))


class PreProcess(RawData):
    """A metaclass that uses RawData attempting several filtering steps.

    Attributes
    ----------
    _invivo_modifications : list
        A list of invivo modifications.
    pep_sel : DataFrame
        DataFrame to filter raw peptide groups DataFrame by indexing.
    prot_sel : DataFrame
        DataFrame to filter raw proteins DataFrame by indexing.
    mitodex : DataFrame
        For filtering mitochondrial peptide groups DataFrame by indexing.
    nonmitodex : DataFrame
        For filtering non-mitochondrial peptide groups DataFrame by indexing.

    Notes
    -----

    See Also
    --------
    omin.utils.SelectionTools.vLook
    omin.utils.SelectionTools.masterCleanse
    """

    """Handles Proteome Discoverer search results.
    """
    def __init__(self, file_list=None, peptides_file=None, proteins_file=None,
                 modifications=None):
        """Initalize instance of PreProcess class.

        Parameters
        ----------
        peptides_file : str
            The location of peptide groups file.
        proteins_file : str
            The location of peptide groups file.
        modifications : list
            List of derived or given modifications.
        """
        # Initalize the RawData base class.
        super(PreProcess, self).__init__(file_list, peptides_file,
                                         proteins_file)
        # Find invivo modifications
        invivos = SelectionTools.findInVivoModifications(self.raw_peptides)
        self._invivo_modifications = invivos
        # Set modifications with given or derived.
        modifications = modifications or self._invivo_modifications
        # Create 2 Dataframes that map specific peptide or protien uniprot ID
        # to it's relevent mitocarta index. Simillar to vlookup in excel.
        # pep_sel, prot_sel = FilterTools.vLook(self.raw_peptides,
        #                                       self.raw_proteins)
        pep_sel, prot_sel = FilterTools.bridge(self.peptide_groups.raw,
                                               self.proteins.raw)

        # Used to index filter to statistically relevent PeptideGroups
        self.pep_sel = pep_sel
        # Used to index filter to statistically relevent Proteins
        self.prot_sel = prot_sel
        # MitoCarta calls made
        mito, nonmito = mitoCartaCall.mitoCartaPepOut(self,
                                                      mods=modifications,
                                                      dex=True)
        # Mitocarta hits DataFrame
        self.mitodex = mito
        # Mitocarta non hits DataFrame
        self.nonmitodex = nonmito
        # Delete the temporary varibles.
        del mito, nonmito, pep_sel, prot_sel
        # Garbage collect.
        gc.collect()
        # Create a unified index of mitocarta calls.
        self.unidex = None
        try:
            unidex = pd.concat([self.mitodex, self.nonmitodex]).sort_index()
            self.unidex = unidex.reindex(index=self.peptide_groups.raw.index)
            del unidex
            gc.collect()
        except Exception:
            pass
        self.master_index = None
        try:
            # Get the Gene Symbol
            gi = self.proteins.raw["Gene ID"].ix[self.prot_sel.index]
            gi.index = self.peptide_groups.raw.index
            # Get Gene Description.
            ga = self.proteins.raw["Description"].ix[self.prot_sel.index]
            ga.index = self.peptide_groups.raw.index
            mdex = pd.concat([self.unidex.ix[:, -4],
                              self.unidex.ix[:, -2:]],
                             axis=1)
            mdex = pd.DataFrame(mdex.fillna(0.0), dtype="bool")

            self.master_index = pd.concat([self.pep_sel,
                                           gi,
                                           ga,
                                           self.peptide_groups.raw.Modifications,
                                           self.peptide_groups.raw["Modifications in Proteins"],
                                           self.peptide_groups.raw.Sequence,
                                           mdex], axis=1)
            self.numbers['mitocarta_hits'] = self.master_index.MitoCarta2_List.sum()
            del gi, mdex
            gc.collect()

        except Exception:
            print("Could not create master_index")
            pass
        # Find the number of inputs.
        self._input_number = SelectionTools.find_number_input(self.raw_peptides)
        # Find the number of PTMs present in the PeptideGroups
        self._ptm_fraction_numbers = SelectionTools.find_fractions(self.raw_peptides)


class Process(PreProcess):
    """A metaclass that uses PreProcess attempting several normalization steps.

    Attributes
    ----------
    _invivo_modifications : list
        A list of invivo modifications.
    pep_sel : DataFrame
        DataFrame to filter raw peptide groups DataFrame by indexing.
    prot_sel : DataFrame
        DataFrame to filter raw proteins DataFrame by indexing.
    mitodex : DataFrame
        For filtering mitochondrial peptide groups DataFrame by indexing.
    nonmitodex : DataFrame
        For filtering non-mitochondrial peptide groups DataFrame by indexing.

    Notes
    -----
    FIXME: Include paragraph description of the types of filtering.

    See Also
    --------
    omin.utils.SelectionTools.vLook
    omin.utils.SelectionTools.masterCleanse
    """

    def __init__(self, file_list=None, peptides_file=None, proteins_file=None,
                 modifications=None):
        """Initalize Process class.

        Parameters
        ----------
        peptides_file : str
            The location of peptide groups file.
        proteins_file : str
            The location of peptide groups file.
        modifications : list
            List of derived or given modifications.
        """
        self.normalized = None
        super(Process, self).__init__(file_list, peptides_file, proteins_file)

        # FIXME: Make the selection more specifically target abundance columns
        if self._input_number > 0:
            inp_notify = "{} input fraction(s) found. Normalizing now..."
            inp_notify = inp_notify.format(self._input_number)
            print(inp_notify)
            try:
                self.normalized = NormalizedToInput(self)
            except Exception:
                print("omin.normalize.toInput.NormalizedToInput FAILED.")

        elif self.raw_peptides.columns.str.contains("pool|control",
                                                    case=False).any():

            print("Pool columns. Omin will attempt to compare the data to it.")
            try:
                self.normalized = NormalizedToPool(self.raw_peptides)
            except Exception:
                print("omin.normalize.toPool.NormalizedToPool FAILED.")
        else:
            # FIXME: Make this a place where the user could specify.
            print("Cannot find anything to normalize to.")
