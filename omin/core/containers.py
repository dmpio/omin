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
from ..normalize.toPool import NormalizedToPool
from ..normalize.toInput import NormalizedToInput
from ..databases import mitoCartaCall
from ..utils.pandas_tools import pd


# FIXME: Find out what words are on limits.


class Container(DataLoader, Handle):
    """Base class for Proteome Discoverer raw files."""

    def __init__(self, *args, **kwargs):
        """Initalize base class for Proteome Discoverer raw files.
        """
        super(Container, self).__init__(**kwargs)

        self.metadata["file_name"] = self.file_name
        self.metadata["file_path"] = self.file_path
        self.metadata["file_ext"] = self.file_ext


class ProteomeDiscovererRaw(Container):
    """Base class for Proteome Discoverer raw files."""

    def __init__(self, *args, **kwargs):
        """Initalize base class for Proteome Discoverer raw files.
        """
        super(ProteomeDiscovererRaw, self).__init__(low_memory=False, delimiter='\t', **kwargs)

        self.thermo_catagory_values = set([i.split(":")[0] for i in self.raw.columns])
        self.thermo_catagory_keys = list(map(StringTools.remove_punctuation, self.thermo_catagory_values))
        self.thermo_catagory_keys = list(map(lambda x: x.strip().replace(" ", "_"), self.thermo_catagory_keys))
        self.thermo_catagory = dict(zip(self.thermo_catagory_keys, self.thermo_catagory_values))
        for k,v in self.thermo_catagory.items():
            self.__dict__[k] = self.raw.filter(regex=v+":")

        self.metadata['total_ids_no_filter'] = self.raw.shape[0]
        self.metadata['total_id_with_quant'] = self.Abundance.dropna(axis=0, how='all').shape[0]


class PeptideGroups(ProteomeDiscovererRaw):
    """Base class for Peptide Groups."""

    def __init__(self, *args, **kwargs):
        """Initialize the base class."""
        ProteomeDiscovererRaw.__init__(self, title="Select peptide groups file", *args, **kwargs)
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

    def __init__(self, *args, **kwargs):
        """Initalize the base class."""
        ProteomeDiscovererRaw.__init__(self, title="Select peptide groups file", *args, **kwargs)
        # Create the master_index
        self.master_index = None
        # Filter for master proteins and high confidence.
        self._high_confidence = FilterTools.high_confidence(self.raw)
        self.master_high_confidence = FilterTools.is_master_protein(self._high_confidence)
        self.metadata['high_confidence_ids'] = self.master_high_confidence.shape[0]
