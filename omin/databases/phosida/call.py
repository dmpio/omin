# LICENSE
# -------

# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

import os
import json


class Phosida(object):
    """Load phosphorylation motifs from phosida.
    """
    # Load the string for the phosida database.
    db_name = 'phosida/phosida_annotated.json'
    # Get this files dir as a string.
    this_dir, _ = os.path.split(__file__)
    # Create the path string.
    file_path = os.path.join(this_dir, db_name)
    # Load the JSON database.
    data = json.load(open(file_path, "r"))
