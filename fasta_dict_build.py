# -*- coding: utf-8 -*-

import os
import urllib.request
import gzip

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

if database_name not in os.listdir("Databases"):
    # print("Downloading", database_name)
    urllib.request.urlretrieve(uniprot, "Databases\\"+zipped_database_name)
    unGzipper(zipped_database_name, database_name, "Databases")

# if zipped_database_name in os.listdir("Databases"):
os.remove("Databases\\"+zipped_database_name)

# import pickle
#
# #BUILD PROTEOME FASTA SEQUENCE DICTIONARY
# fastafile = "Databases\Mouse_UP000000589_20160623.fasta"#Whole mouse proteome @ 31.8 MB
#
# genedict = {}#Dictionary with Uniprot IDs as keys and fasta sequences as values
#
# with open(fastafile,"r") as infile:
#     for line in infile:
#         if line[0] == ">":#Tests each line to see if it is a new sequence
#             genename = line.split("|")[1]#Grabs the genename
#             genedict[genename] = list([line[:-1]])#Makes the first entry for key
#         else:
#             genedict[genename].append(line[:-1])#Appeneds fasta lines for key
#
# pickle.dump(genedict,open("genedict.p","wb"))
