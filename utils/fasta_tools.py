# -*- coding: utf-8 -*-
"""Tools for investigating with fasta and uniprot."""

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

from urllib.error import HTTPError
from urllib.request import urlopen


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


class UniprotTools(object):

    @classmethod
    def getFasta(cls, upID):
        """Takes a uniprot ID and returns a string.

        Parameters
        ----------
        upID : str

        Returns
        -------
        fasta : str

        """
        try:
            fasta = str
            link = "http://www.uniprot.org/uniprot/"+upID+".fasta"
            try:
                urlout = urlopen(link)
                # print("hey")

            except HTTPError:
                print(
                      "Try again, the formatting may be off.",
                      "Try removing dashes."
                      )

            # fasta = urlout.read().decode("utf-8")
            fasta = urlout.read().decode("utf-8").split("\n")
            return fasta
        except TypeError as err:
            print("Please enter UniProt ID as string.")
