# -*- coding: utf-8 -*-
"""Call the Mitocarta 2.0 database."""

# LICENSE
# -------

# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

import os
from pandomics import pandas as pd


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
