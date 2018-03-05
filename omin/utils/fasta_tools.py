# -*- coding: utf-8 -*-
"""Tools for investigating with fasta and uniprot."""

# LICENSE
# -------

# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

import re
from urllib.error import HTTPError
from urllib.request import urlopen

# FIXME: DEPRECATE ALL

class FastaTools(object):
    """A collection of tools for handling FASTA formatted strings."""

    @classmethod
    def fasta2Seq(cls, fasta):
        """Join strings in list into one string.

        Parameters
        ----------
        fasta : str

        Returns
        -------
        one_string : str
        """
        one_string = "".join(fasta[1:])
        return one_string

    @staticmethod
    def seqNum(cls, seq):
        """Take a sequence and returns the number for each amino acid.

        Notes
        -----
            Review this functions useage and rewrite annotation.

        Parameters
        ----------
        seq : str

        Returns
        -------
        seq_num : str

        """
        seq_num = [seq[i]+str(i+1) for i in range(len(seq))]
        return seq_num

# === UNIPROT TOOLS ===


class UniProtTools(object):
    """Tools for querying the UniProt database."""

    id_rx = re.compile('ID;\s(\d+);')

    @classmethod
    def get_uniprot(cls, accession, file_format='.fasta'):
        """Retrieve UniProt data for a given accession."""
        if file_format[0] != '.':
            file_format = '.' + file_format

        url = "http://www.uniprot.org/uniprot/"
        req = ''.join([url, accession, file_format])
        try:
            data = urlopen(req)
            return data.read().decode()
        except HTTPError:
            print('Entry :', accession, 'not found.')
            return None

    @classmethod
    def get_gene_id(cls, accession):
        """Retrieve the gene id for a given accession."""
        result = cls.get_uniprot(accession, file_format='.txt')
        result = cls.id_rx.findall(result)
        result = result[0]
        return result

    @classmethod
    def go_anno(cls, accession, desc=False):
        """Retrieve GO annotations for a given accession."""
        txt = cls.get_uniprot(accession, file_format='.txt')
        if desc:
            rx = re.compile('GO:\d+;.+\n')
        else:
            rx = re.compile('(GO:\d+);')
        results = rx.findall(txt)
        return results
