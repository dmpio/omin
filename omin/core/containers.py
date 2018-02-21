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
class ProteomeDiscovererRaw(DataLoader, Handle):
    """Base class for Proteome Discoverer raw files."""

    def __init__(self, *args, **kwargs):
        """Initalize base class for Proteome Discoverer raw files.
        """
        super(ProteomeDiscovererRaw, self).__init__(**kwargs)
        # DataLoader.__init__(self, *args, **kwargs)
        # self.raw = raw_data.copy()

        #self.abundance = self.raw.filter(regex="Abundance:")
        #self.numbers['total_ids'] = self.raw.shape[0]


class PeptideGroups(ProteomeDiscovererRaw):
    """Base class for Peptide Groups."""

    def __init__(self, *args, **kwargs):
        """Initialize the base class."""
        # super(PeptideGroups, self).__init__(*args, **kwargs)
        ProteomeDiscovererRaw.__init__(self, *args, **kwargs)
        if 'Modifications' in self.raw.columns:
            # Replace any NaNs that might be present in Modifications.
            self.raw.Modifications.fillna('', inplace=True)
        else:
            print("No Modifications column found in Peptide Groups data.")
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
        # Filter for master proteins and high confidence.
        self._high_confidence = FilterTools.high_confidence(self.raw)
        self.master_high_confidence = FilterTools.is_master_protein(self._high_confidence)
