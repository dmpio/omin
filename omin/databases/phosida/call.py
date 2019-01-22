# LICENSE
# -------

# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

import os
import re

import pandas as pd
import numpy as np

from glob import glob

from urllib.request import urlopen, Request
from urllib.parse import urlencode

from bs4 import BeautifulSoup as bs


class Phosida(object):
    """Load phosphorylation motifs from phosida.
    """
    # Load the string for the phosida database.
    db_name = 'phosida/PHOSIDA posttranslational modification database.html'
    # Get this files dir as a string.
    this_dir, _ = os.path.split(__file__)
    # Create the path string.
    file_path = os.path.join(this_dir, db_name)

    with open(file_path, "r") as f:
        data = f.read()

    soup = bs(data, "lxml")

    motifs = list(filter(lambda x:len(x)>0,[i.text.strip() for i in soup.find_all("tr")]))
    motifs = motifs[1:-16]
    motifs = [i.split("\n") for i in motifs]


    phosida = dict()

    for i in motifs:
        if len(i)>1:
            ref = i
            phosida[ref[0]] = [ref[1]]
        else:
            phosida[ref[0]] = phosida[ref[0]] + i

    # Remove all digits
    for k in phosida.keys():
        phosida[k] = [re.sub("\d+", "", i) for i in phosida[k]]


    def process(list_to_process):
        result=[]
        list_to_process = [i.split("-") for i in list_to_process]

        for k in list_to_process:
            for j in range(len(k)):
                if "/" in k[j]:
                    j_prime = k[j].replace("/", "")
                    k[j] = "[{}]".format(j_prime)

            k_prime = "".join(k)

            k_prime = k_prime.replace("X", ".")

            result.append(k_prime)

        return result

    # process(phosida["PKA"])

    for k in phosida.keys():
        phosida[k] = process(phosida[k])
