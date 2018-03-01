# -*- coding: utf-8 -*-
"""Call the Mitocarta 2.0 database."""

# LICENSE
# -------

# Copyright 2017 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach,
# Blair Chesnut, and Elizabeth Hauser.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files, (the software)), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions: The above copyright
# notice and this permission notice shall be included in all copies or
# substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS",
# WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
# TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM. OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import os
# import pickle
# import pandas as pd
from ...utils.pandas_tools import pd


class MitoCartaTwo(object):
    """MitoCarta2 database calls."""

    # Load the string for the mitocarta database.
    db_name = 'complete_mitocarta_2.p'
    # Get this files dir as a string.
    this_dir, _ = os.path.split(__file__)
    # Create the path string.
    carta_file_path = os.path.join(this_dir, db_name)
    # load the MitoCarta 2.0 database
    data = pd.read_pickle(carta_file_path)
    # Change the incorrectly labeled column from MouseGeneID to EntrezGeneID
    data.rename(columns={"MouseGeneID": "EntrezGeneID"}, inplace=True)
    # Collect the essential columns.
    essential = data[["EntrezGeneID", "MitoCarta2_List", "Matrix", "IMS"]]
    # Load the modified mitocarta2.0 database.
    ukb = 'uniprotkb_mitocarta2.p.gz'
    ukb_file_path = os.path.join(this_dir, ukb)
    ukb2mito = pd.read_pickle(ukb_file_path, compression='gzip')

    @classmethod
    def look_up(cls, syn_list):
        """Look-up a given accession number in MitoCarta 2.0."""
        return cls.data[cls.data.Synonyms.str.contains('|'.join(syn_list))]
