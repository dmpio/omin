# -*- coding: utf-8 -*-
"""Call the KEGG Database"""

# LICENSE
# -------

# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

from urllib.request import urlopen, Request
from bs4 import BeautifulSoup as bs


def make_soup(url):
    """Returns soup for url."""
    req = urlopen(Request(url))

    data = req.read()

    soup = bs(data, "lxml")

    return soup


def retrieve_kegg_database():
    url = "https://www.kegg.jp/kegg/pathway.html"

    soup = make_soup(url)

    tbls = [i for i in soup.find_all("div",  attrs={"class":"list"})]

    kegg_pathways = dict()

    for i in range(1, len(tbls)):
        keys = [i.text for i in tbls[i].find_all("dt")]

        values = [i.text for i in tbls[i].find_all("dd")]

        d = dict(zip(keys, values))

        kegg_pathways = {**kegg_pathways, **d}

    return kegg_pathways
