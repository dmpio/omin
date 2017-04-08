# -*- coding: utf-8 -*-
"""
Copyright 2017 James Draper, Paul Grimsrud, Deborah Muoio

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files, Omics Modeling Integrating
Normalization (OMIN), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom
the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM.
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

# FIXME: Combine the contents of this script as a class of it's own in the routines.py file.

import os
import pickle
import gzip
import urllib.request


uniprot = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000589_10090.fasta.gz"

for fold in uniprot.split("/"):
    if fold.startswith("UP000000589"):
        if fold.endswith("gz"):
            zipped_database_name = fold
            database_name = ".".join(zipped_database_name.split(".")[:-1])
        else:
            print("Could not find:", fold)


def unGzipper(compressed_datafile, decompressed_data_file_name, parent=None):
    """Decompresses '.gz' files.

    Parameters
    ----------
    compressed_datafile : str
    decompressed_data_file_name : str
    parent : str

    Returns
    -------
    ungzipped file
    """
    with gzip.open(parent + "\\" + compressed_datafile, 'rb') as infile:
        with open(parent + "\\" + decompressed_data_file_name, 'wb') as outfile:
            for line in infile:
                outfile.write(line)
    return

if database_name not in os.listdir("databases"):
    # print("Downloading", database_name)
    urllib.request.urlretrieve(uniprot, "databases\\"+zipped_database_name)
    unGzipper(zipped_database_name, database_name, "databases")

if zipped_database_name in os.listdir("databases"):
    os.remove("databases\\"+zipped_database_name)

# === BUILD PROTEOME FASTA SEQUENCE DICTIONARY ===

# Create a dict with Uniprot IDs as keys and fasta sequences as values
genedict = {}

parent = "databases"

if database_name in os.listdir(parent):
    fastafile = "\\".join([parent, database_name])
    with open(fastafile, "r") as infile:
        for line in infile:
            # Tests each line to see if it is a new sequence
            if line[0] == ">":
                genename = line.split("|")[1]  # Grabs the genename
                # Makes the first entry for key
                genedict[genename] = list([line[:-1]])
            else:
                # Appeneds fasta lines for key
                genedict[genename].append(line[:-1])

if len(genedict) is not 0:
    pickle.dump(genedict, open("databases\\mouse_proteome.p", "wb"))
