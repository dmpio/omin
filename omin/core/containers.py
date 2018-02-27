# -*- coding: utf-8 -*-
"""

"""
# FIXME: DOCUMENT OR DIE #DOD

# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

# ----------------
# EXTERNAL IMPORTS
# ----------------
import re
import os
import numpy as np
from guipyter import DataLoader

# ----------------
# INTERNAL IMPORTS
# ----------------
from .base import Handle
# -------------
# UTILS IMPORTS
# -------------
from ..utils import IOTools
from ..utils import StringTools
from ..utils import SelectionTools
from ..utils import FilterTools
from ..utils import IntermineTools
# --------
# DATBASES
# --------
from ..databases import MitoCartaTwo
# import pandas with the pandomics plug-in
from pandomics import pandas as pd

# FIXME: Add the same try and except as seen in the guipyter module.
# # The line below could be used as a backup.
# from ..utils.pandas_tools import pd


# FIXME: Find out what words are on limits.

class Normalized(object):
    """Empty class that catches results from normalization methods.
    """
    def __init__(self, *args, **kwargs):
        for k,v in kwargs.items():
            self.__dict__[k] = v

    def __repr__(self):
        """Show all attributes."""
        return "Attributes: "+", ".join(list(self.__dict__.keys()))


class Container(DataLoader, Handle):
    """Base class for Proteome Discoverer raw files.
    """

    def __init__(self, *args, **kwargs):
        """Initialize base class for Proteome Discoverer raw files.

        This class inherits attributes from the DataLoader and Handle classes respectfully.

        Attributes
        ----------
        raw : pandas.DataFrame
            With pandomics plugin.
        file_path : str
        file_name : str
        file_ext : str
        metadata : dict

        """

        # Initialize the DataLoader class with args and kwargs
        DataLoader.__init__(self, *args, **kwargs)
        # Initialize the Handle class.
        Handle.__init__(self)

        self.metadata["file_name"] = self.file_name
        self.metadata["file_path"] = self.file_path
        self.metadata["file_ext"] = self.file_ext


class MaxQuantRaw(Container):
    """Base class for MaxQuant raw files.

    This class inherits attributes from the Container class.

    HARD HAT AREA: Under construction.
    """
    def __init__(self, *args, **kwargs):
        # FIXME: Add functional algorithms.
        # FIXME: Add as mixin class later on.
        Container.__init__(self, *args, **kwargs)


class ProteomeDiscovererRaw(Container):
    """Base class for Proteome Discoverer raw files.

    This class inherits attributes from the Container class.
    """

    def __init__(self, *args, **kwargs):
        """Initialize base class for Proteome Discoverer raw files.
        """
        # super(ProteomeDiscovererRaw, self).__init__(low_memory=False, delimiter='\t', **kwargs)
        Container.__init__(self, low_memory=False, delimiter='\t', **kwargs)
        # Calling the function below.
        self.expose_thermo_categories()
        # METADATA: Number of IDs no filter.
        self.metadata['total_ids'] = self.raw.shape[0]
        # METADATA: Number of IDs with quantification.
        self.metadata['total_ids_with_quant'] = self.Abundance.dropna(axis=0, how='all').shape[0]


    def expose_thermo_categories(self):
        """Creates attribute DataFrames based on Thermo's categories.
        """
        thermo_category_values = set([i.split(":")[0] for i in self.raw.columns])
        thermo_category_keys = list(map(StringTools.remove_punctuation, thermo_category_values))
        thermo_category_keys = list(map(lambda x: x.strip().replace(" ", "_"), thermo_category_keys))
        thermo_category = dict(zip(thermo_category_keys, thermo_category_values))
        # This loop creates DataFrames for each category that Thermo has decided is important.
        #
        for k,v in thermo_category.items():
            filter_with_colon = self.raw.filter(regex=v+":")
            if filter_with_colon.shape[-1] > 0:
                self.__dict__[k] = filter_with_colon
            else:
                filter_without_colon = self.raw.filter(regex=v)
                self.__dict__[k] = filter_without_colon

    # ------------------
    # STUDY FACTOR TOOLS
    # ------------------
    @property
    def study_factor_table(self):
        """Attempt to retrieve the study factors from Abundance column headers.

        Study Factors should be in the format:
            KO (Genotype), WT (Genotype), Input (Fraction), Phospho (Fraction) ect.

        If study factor cannot be identified then StudyFactor_<ABC...> will be assigned.
        """
        column_names = [i.split(":")[-1].strip() for i in self.Abundance.columns]

        column_names = pd.DataFrame(list(map(lambda x: x.split(','), column_names)))

        # Remove the sample column if present.
        if len(column_names.iloc[:, 1].unique()) < 2:
            column_names = pd.concat([column_names.iloc[:, 0], column_names.iloc[:, 2:]], axis=1)

        # Compile the regular expression to catch things in (parenthesis).
        # FIXME: BIG ASSUMPTION HERE. DOCUMENT OR DIE. #DOD
        rx = re.compile("\((.+)\)")
        study_factors = []
        for i in range(1, len(column_names.columns)):
            ucols = column_names.iloc[:, i].apply(lambda x: rx.findall(x)[0])
            ucols = ucols.unique()
            if len(ucols) == 1:
                study_factor = ucols[0]
                study_factors.append(study_factor)
            else:
                # If study factor cannot be identified then StudyFactor_<ABC...> will be assigned.
                study_factor = "StudyFactor_" + string.ascii_uppercase[i-1]
                study_factors.append(study_factor)

        column_names.columns = ['_TMT_tag'] + study_factors
        # Prepend the fraction number.
        fn = [i.split(":")[1].strip() for i in self.Abundance.columns]
        fn = pd.DataFrame(fn, columns=["_Fn"])

        column_names = pd.concat([fn, column_names],axis=1)
        return column_names

    @property
    def study_factor_dict(self):
        """Return the study_factor_table as a flattened dict.
        """
        study_factor_dict = dict()
        for i in range(2, self.study_factor_table.columns.shape[0]):
            usf = self.study_factor_table.iloc[:, i].unique()
            usf = list(map(lambda x: x.strip(), usf))
            study_factor_dict[self.study_factor_table.iloc[:, i].name] = usf
        return study_factor_dict

    @property
    def tmt_plex_number(self):
        """Returns the TMT plex number.
        """
        plex_number = 0

        try:
            plex_number = self.study_factor_table._TMT_tag.unique().shape[0]

        except Exception as err:
            print(err)

        return plex_number

    @property
    def input_number(self):
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
            plex_number = self.tmt_plex_number
            # FIXME: Proof this against Inputase ect.
            # To do that that load the regex with ~ [Ii]nput[\w[punctuation]]
            # FIXME: Redundant filtering of raw DataFrames = memory leak
            input_abundance = self.Abundance.filter(regex="[Ii]nput").shape[1]
            number_input = input_abundance/plex_number
        except Exception:
            print("utils.SelectionTools.find_number_input failed")
        return number_input

    def fraction_tag(self, fraction_number=None):
        """Return a characteristic string for a fraction.
        """
        by_fn = self.study_factor_table.loc[self.study_factor_table._Fn == fraction_number].iloc[:, 2:]

        rx = re.compile('(.+)\s\(.+\)')

        tag_for_fraction = ""
        for i in by_fn:
            terms = by_fn[i].unique()
            terms = list(map(lambda x:x.strip(), terms))
            if len(terms) == 1:
                term = terms[0]
                term = rx.findall(term)
                if len(term) == 1:
                    term = term[0].lower()
                    tag_for_fraction = "_".join([tag_for_fraction, term])

        return tag_for_fraction


class PeptideGroups(ProteomeDiscovererRaw):
    """Base class for Peptide Groups.

    Derived from the ProteomeDiscovererRaw class
    """

    def __init__(self, filepath_or_buffer, *args, **kwargs):
        """Initialize the base class."""
        filepath_or_buffer = filepath_or_buffer or None
        ProteomeDiscovererRaw.__init__(self, filepath_or_buffer=filepath_or_buffer, title="Select peptide groups file", *args, **kwargs)

        # Fill all of the modifications with NaNs with a blank string
        try:
            self.raw.Modifications.fillna('', inplace=True)
        except Exception as err:
            print("No Modifications column found in Peptide Groups data.", err)
        # Create the master_index
        self.master_index = None
        self.set_master_index()

        # Declare variable
        self._in_vivo_modifications = []
        self.set_in_vivo_modifications()


    def set_master_index(self):
        """Attempt to set the master index for the Proteins.
        """
        # FIXME: Add try and except for each of these columns.
        master_index_components = ["Sequence", "Modifications", "Modifications in Proteins"]

        # Start the master index with the the first master protein accession.
        self.master_index = pd.DataFrame(self.raw["Master Protein Accessions"].first_member())
        # Change the name of "Master Protein Accessions" to "Accession".
        try:
            self.master_index.rename(columns={"Master Protein Accessions":"Accession"}, inplace=True)
        except Exception as err:
            pass

        for i in master_index_components:
            try:
                self.master_index = pd.concat([self.master_index, self.raw[i]],axis=1)
            except Exception as err:
                # print(err)
                pass
        return


    def set_in_vivo_modifications(self):
        """FIXME: Add docs.
        """
        # FIXME: Rethink this part. Should this be done in the ProteomeDiscovererRaw class?
        # Find invivo modifications.
        in_vivo_mods = SelectionTools.findInVivoModifications(self.raw)
        # Check if in_vivo modifications is None.
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
                    # METADATA: Total Peptide Group IDs
                    self.metadata[mod+"_total_peptide_ids"] = df.shape[0]
                    self.__dict__[mod] = df
            else:
                pass
        else:
            pass

    @property
    def study_factor_with_input(self):
        """Return the study factor that contains "Input" or "input".

        PROTIP: In your when creating study factors stick with the conventions: Input (Fraction), Acetyl (Fraction), ect.
        NOTE: For this to work correctly there must be only one study factor that contains the terms: "Input" or "input".

        """
        #FIXME: Inputase-proof this function or provide informative error handling with messages.
        rx = re.compile("[Ii]nput")
        study_factor_with_input = None
        for k,v in self.study_factor_dict.items():
            #FIXME: BIG ASSUMPTION HERE -> There will be only one study factor that contains a term including [Ii]nput.
            li = bool(sum(list(map(lambda x:bool(len(rx.findall(x))), v))))
            if li:
                study_factor_with_input = k
            else:
                pass
        return study_factor_with_input

    @property
    def load_normalized(self):
        """Return object with load normalzed pd.DataFrames as attributes.
        """
        input_mask = self.study_factor_table[self.study_factor_with_input].str.contains("[Ii]nput")
        # FIXME: Replace the bobo method to find input fraction numbers
        number_input_fractions = len(self.study_factor_table.loc[input_mask]._Fn.unique())

        inps = self.study_factor_table.loc[input_mask]
        inps = [inps.loc[inps._Fn.str.contains(i)] for i in inps._Fn.unique()]
        frcs = self.study_factor_table.loc[~input_mask]
        frcs = [frcs.loc[frcs._Fn.str.contains(i)] for i in frcs._Fn.unique()]

        linked = []
        for i in inps:
            for j in frcs:
                score = sum(sum(j.as_matrix() == i.as_matrix()))
                link = pd.DataFrame(i.index, columns=["Link"], index=j.index)
                score_df = pd.DataFrame([i.shape[0]*[score]]).T
                score_df.columns = ["Score"]
                score_df.index = j.index
                j_prime = pd.concat([j, link, score_df], axis=1)
                linked.append(j_prime)

        scores = np.array([i.Score.unique()[0] for i in linked])
        linked = list(filter(lambda x:x.Score.unique()[0] == scores.max(), linked))

        normalized = dict()
        for link in linked:
            inp = self.Abundance[self.Abundance.columns[link.Link]]
            other = self.Abundance[self.Abundance.columns[link.index]]
            other_label = self.fraction_tag(link._Fn.unique()[0])
            #other_label = [i.split(":")[1] for i in other.columns]
            load_normalized = other.normalize_to(inp)
            normalized[other_label] = load_normalized
        # Throwing the normalized dict into the Normalized class
        return Normalized(**normalized)


class Proteins(ProteomeDiscovererRaw):
    """Base class for Proteins.
    Derived from the ProteomeDiscovererRaw class
    """

    def __init__(self, filepath_or_buffer, attempt_rescue_entrez_ids=True, *args, **kwargs):
        """

        Parameters
        ----------
        filepath_or_buffer: str or _io.TextWrapper

        attempt_rescue_entrez_ids: Bool
            Defaults to True.

        """
        filepath_or_buffer = filepath_or_buffer or None
        # FIXME: Add verbose argument.
        ProteomeDiscovererRaw.__init__(self, filepath_or_buffer=filepath_or_buffer, title="Select proteins file", *args, **kwargs)

        # Filter for master proteins and high confidence.
        self._high_confidence = FilterTools.high_confidence(self.raw)
        self.master_high_confidence = FilterTools.is_master_protein(self._high_confidence).copy()
        # METADATA: Number of high confidence protein.
        self.metadata['high_confidence_ids'] = self.master_high_confidence.shape[0]

        # Find and relabel the Entrez Gene ID column.
        self.set_entrez()

        # Create the master_index
        self.set_master_index()
        # Attempt to rescue Entrez Gene IDs
        if attempt_rescue_entrez_ids:
            IntermineTools.rescue_entrez_ids(self.master_index)
        # Attach MitoCarta2 data to the master_index.
        self.add_database(MitoCartaTwo.essential)

    def set_entrez(self):
        """Find and relabel the Entrez Gene ID column.
        """
        if "Gene ID" in self.master_high_confidence: # PD2.1
            self.master_high_confidence.rename(columns={'Gene ID':'EntrezGeneID'}, inplace=True)

        if "Entrez Gene ID" in self.master_high_confidence: # PD2.2
            self.master_high_confidence.rename(columns={'Entrez Gene ID':'EntrezGeneID'}, inplace=True)

    def set_master_index(self):
        """Attempt to set the master index for the Proteins.
        """
        # FIXME: Add try and except for each of these columns.
        try:
            self.master_index = pd.concat([self.master_high_confidence.Accession, self.master_high_confidence.EntrezGeneID, self.master_high_confidence.Description], axis=1)
        except Exception as err:
            print(err)

    def add_database(self, DataFrame):
        """Add databases to master_index.
        """
        try:
            self.master_index.dropna(inplace=True)
            self.master_index.EntrezGeneID = self.master_index.EntrezGeneID.first_member().apply(np.int64)
            self.master_index = self.master_index.merge(DataFrame.copy(), on="EntrezGeneID", how="left")

            self.master_index.index = self.master_index.index
            self.master_index.fillna(False, inplace=True)
        except Exception as err:
            print(err)
