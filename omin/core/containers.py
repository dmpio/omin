# -*- coding: utf-8 -*-
# Copyright 2017 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

import re
import os
import numpy as np
from guipyter import DataLoader

# ----------------
# INTERNAL IMPORTS
# ----------------
from .base import Handle
# UTILS IMPORTS
from ..utils import IOTools
from ..utils import StringTools
from ..utils import SelectionTools
from ..utils import FilterTools
from ..utils import IntermineTools

# DATBASES
from ..databases import MitoCartaTwo
# import pandas with the pandomics plug-in
from pandomics import pandas as pd

# FIXME: Add the same try and except as seen in the guipyter module.
# # The line below could be used as a backup.
# from ..utils.pandas_tools import pd


# FIXME: Find out what words are on limits.


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
        super(ProteomeDiscovererRaw, self).__init__(low_memory=False, delimiter='\t', **kwargs)

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

        self.metadata['total_ids_no_filter'] = self.raw.shape[0]
        self.metadata['total_id_with_quant'] = self.Abundance.dropna(axis=0, how='all').shape[0]

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
        return column_names

    @property
    def study_factor_dict(self):
        """Return the study_factor_table as a flattened dict.
        """
        study_factor_dict = dict()
        for i in range(1, self.study_factor_table.columns.shape[0]):
            usf = self.study_factor_table.iloc[:, i].unique()
            usf = list(map(lambda x: x.strip(), usf))
            study_factor_dict[self.study_factor_table().iloc[:, i].name] = usf
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


class PeptideGroups(ProteomeDiscovererRaw):
    """Base class for Peptide Groups.

    Derived from the ProteomeDiscovererRaw class
    """

    def __init__(self, filepath_or_buffer, *args, **kwargs):
        """Initialize the base class."""
        filepath_or_buffer = filepath_or_buffer or None
        ProteomeDiscovererRaw.__init__(self, filepath_or_buffer=filepath_or_buffer, title="Select peptide groups file", *args, **kwargs)
        # Create the master_index
        self.master_index = None
        # Fill all of the modifications with NaNs with a blank string
        try:
            self.raw.Modifications.fillna('', inplace=True)
        except Exception as err:
            print("No Modifications column found in Peptide Groups data.", err)

        # FIXME: Rethink this part. Should this be done in the ProteomeDiscovererRaw class?
        # Find invivo modifications.
        in_vivo_mods = SelectionTools.findInVivoModifications(self.raw)
        # Declare variable
        self._in_vivo_modifications = []
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
        # Create the master_index
        self.master_index = None
        # Filter for master proteins and high confidence.
        self._high_confidence = FilterTools.high_confidence(self.raw)
        self.master_high_confidence = FilterTools.is_master_protein(self._high_confidence).copy()
        # METADATA: Number of high confidence protein.
        self.metadata['high_confidence_ids'] = self.master_high_confidence.shape[0]
        # Find and relabel the Entrez Gene ID column.
        if "Gene ID" in self.master_high_confidence: # PD2.1
            self.master_high_confidence.rename(columns={'Gene ID':'EntrezGeneID'}, inplace=True)

        if "Entrez Gene ID" in self.master_high_confidence: # PD2.2
            self.master_high_confidence.rename(columns={'Entrez Gene ID':'EntrezGeneID'}, inplace=True)

        # FIXME: Add try and except for each of these columns.
        try:
            self.master_index = pd.concat([self.master_high_confidence.Accession, self.master_high_confidence.EntrezGeneID, self.master_high_confidence.Description], axis=1)
        except Exception as err:
            print(err)

        if attempt_rescue_entrez_ids:
            IntermineTools.rescue_entrez_ids(self.master_index)

        # Attach MitoCarta2 data to the master_index.
        try:
            self.master_index.dropna(inplace=True)
            self.master_index.EntrezGeneID = self.master_index.EntrezGeneID.first_member().apply(np.int64)
            self.master_index = self.master_index.merge(MitoCartaTwo.essential.copy(), on="EntrezGeneID", how="left")

            self.master_index.index = self.master_index.index
            self.master_index.fillna(False, inplace=True)
        except Exception as err:
            print(err)
