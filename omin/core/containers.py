# -*- coding: utf-8 -*-
# Copyright 2017 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach,
# Blair Chesnut, and Elizabeth Hauser.

import re
import os
# import guipyter as gptr
from guipyter import DataLoader

# import the Handle super class.
from .base import Handle

from ..utils import IOTools
from ..utils import StringTools
from ..utils import SelectionTools
from ..utils import FilterTools
# from ..normalize.toPool import NormalizedToPool
# from ..normalize.toInput import NormalizedToInput
# from ..databases import mitoCartaCall
from ..databases import MitoCartaTwo
from ..utils.pandas_tools import pd
from ..utils.intermine_tools import IntermineTools

# FIXME: Find out what words are on limits.


class Container(DataLoader, Handle):
    """Base class for Proteome Discoverer raw files."""

    def __init__(self, *args, **kwargs):
        """Initalize base class for Proteome Discoverer raw files.
        """
        # super(Container, self).__init__(**kwargs)
        DataLoader.__init__(self, *args, **kwargs)
        Handle.__init__(self)

        self.metadata["file_name"] = self.file_name
        self.metadata["file_path"] = self.file_path
        self.metadata["file_ext"] = self.file_ext


class ProteomeDiscovererRaw(Container):
    """Base class for Proteome Discoverer raw files."""

    def __init__(self, *args, **kwargs):
        """Initalize base class for Proteome Discoverer raw files.
        """
        super(ProteomeDiscovererRaw, self).__init__(low_memory=False, delimiter='\t', **kwargs)

        thermo_catagory_values = set([i.split(":")[0] for i in self.raw.columns])
        thermo_catagory_keys = list(map(StringTools.remove_punctuation, thermo_catagory_values))
        thermo_catagory_keys = list(map(lambda x: x.strip().replace(" ", "_"), thermo_catagory_keys))
        thermo_catagory = dict(zip(thermo_catagory_keys, thermo_catagory_values))
        for k,v in thermo_catagory.items():
            filter_with_colon = self.raw.filter(regex=v+":")
            if filter_with_colon.shape[-1] > 0:
                self.__dict__[k] = filter_with_colon
            else:
                filter_without_colon = self.raw.filter(regex=v)
                self.__dict__[k] = filter_without_colon

        self.metadata['total_ids_no_filter'] = self.raw.shape[0]
        self.metadata['total_id_with_quant'] = self.Abundance.dropna(axis=0, how='all').shape[0]


class PeptideGroups(ProteomeDiscovererRaw):
    """Base class for Peptide Groups."""

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

        # FIXME: Rethink this part.
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
                    # METADATA: Total Peptide Group IDs
                    self.metadata[mod+"_total_peptide_ids"] = df.shape[0]
                    self.__dict__[mod] = df
            else:
                pass
        else:
            pass


class Proteins(ProteomeDiscovererRaw):
    """Base class for Proteins."""

    def __init__(self, filepath_or_buffer, attempt_rescue_entrez_ids=True, *args, **kwargs):
        """Initalize the base class."""
        filepath_or_buffer = filepath_or_buffer or None
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

        self.master_index = pd.concat([self.master_high_confidence.Accession, self.master_high_confidence.EntrezGeneID], axis=1)

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
