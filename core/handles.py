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

import pandas as pd
from ..utils import SelectionTools
from ..normalize.toPool import NormalizedToPool
from ..normalize.toInput import NormalizedToInput
from ..databases import mitoCartaCall


class RawData(object):
    """Converts Proteome Discoverer .txt files into pandas DataFrames

    Attributes
    ----------
    peptides :  DataFrame
        Raw data from Proteome Discoverer peptides data.
    proteins : DataFrame
        Raw data Proteome Discoverer corresponding proteins data.
    """
    def __init__(self, peptides_file, proteins_file):
        """Loads the raw data for peptides_file and proteins_file as DataFrames
        that are contained as attributes.

        Note
        ----
            Please make sure that your files are in your current working
            directory. If you are working in jupyter notebook please put the a
            copy of the peptides and proteins in the same directory as the
            notebook file.

        Parameters
        ----------
        peptides_file : str
            Name of Proteome Discoverer Peptides file as string.
        proteins_file : str
            Name of Proteome Discoverer Proteins file as string.

        Examples
        --------
        Load file strings:
        >>>peptides_file = "mydatafolder/peptides.txt"
        >>>proteins_file = "mydatafolder/proteins.txt"

        Create RawData object:
        >>>raw_data = RawData(peptides_file,proteins_file)

        """
        self.raw_peptides = pd.read_csv(peptides_file,
                                        delimiter="\t",
                                        low_memory=False)

        self.raw_proteins = pd.read_csv(proteins_file,
                                        delimiter="\t",
                                        low_memory=False)

        self._numbers = (self.raw_peptides.shape, self.raw_proteins.shape)

    def __repr__(self):
        """Show all attributes.
        """
        return "Attributes: "+", ".join(list(self.__dict__.keys()))

    def showQuant(self):
        """Returns tuple (peptides DataFrame shape, proteins DataFrame shape).
        """
        return self._numbers


class PreProcess(RawData):
    """
    """
    def __init__(self, peptides_file, proteins_file, modifications=None,
                 genotype=None, treatments=None):

        # Initalize the RawData base class.
        super(PreProcess, self).__init__(peptides_file, proteins_file)

        # modifications = modifications or ["Acetyl", "Phospho"]
        modifications = modifications or SelectionTools.findInVivoModifications(self.raw_peptides)
        # Create 2 Dataframes that map specific peptide or protien uniprot ID
        # to it's relevent mitocarta index. Simillar to vlookup in excel
        pep_sel, prot_sel = SelectionTools.vLook(self.raw_peptides,
                                                 self.raw_proteins,
                                                 modifications)
        self.pep_sel = pep_sel
        self.prot_sel = prot_sel

        mito, nonmito = mitoCartaCall.mitoCartaPepOut(self,
                                                      mods=modifications,
                                                      dex=True)
        self.mitodex = mito
        self.nonmitodex = nonmito


class Process(PreProcess):
    """Formerly omin.Experiment
    """
    def __init__(self, peptides_file, proteins_file, modifications=None,
                 genotype=None, treatments=None):
        """
        """
        #
        # modifications = modifications or ["Acetyl", "Phospho"]
        # # Initalize the RawData base class.
        super(Process, self).__init__(peptides_file, proteins_file)
        # # Create 2 Dataframes that map specific peptide or protien uniprot ID
        # # to it's relevent mitocarta index. Simillar to vlookup in excel
        # pep_sel, prot_sel = SelectionTools.vLook(self.raw_peptides,
        #                                          self.raw_proteins,
        #                                          modifications)
        # self.pep_sel = pep_sel
        # self.prot_sel = prot_sel

        # NORMALIZING TO INPUT
        if self.raw_peptides.columns.str.contains("Input", case=False).any():

            print("Input fraction found. omin will attempt to normalize the data to it.")

            try:
                self.input_fraction_numbers = SelectionTools.find_number_input(raw_peptides)
                self.normalized = NormalizedToInput(self.raw_peptides, self.raw_proteins)

            except Exception:
                print("Something went wrong. Please check to make sure you data is formatted correctly.")

        elif self.raw_peptides.columns.str.contains("pool|control", case=False).any():

            print("Pool or control columns found. omin will attempt to normalize the data to it.")

            try:
                self.normalized = NormalizedToPool(self.raw_peptides)

            except Exception:
                print("Something went wrong. Please check to make sure you data is formatted correctly.")
        else:
            # FIXME: Make this a place where the use could specify.
            print("Cannot find anything to normalize to. Check to see that your data fits omins conventions.")

            self.normalized = None
