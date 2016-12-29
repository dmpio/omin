# -*- coding: utf-8 -*-
import pandas as pd
from omin.utils import SelectionTools

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
        self.peptides = pd.read_csv(peptides_file,
                                    delimiter="\t",
                                    low_memory=False)

        self.proteins = pd.read_csv(proteins_file,
                                    delimiter="\t",
                                    low_memory=False)

        self._numbers = (self.peptides.shape, self.proteins.shape)

    def __repr__(self):
        """Show all attributes.
        """
        return "Attributes: "+", ".join(list(self.__dict__.keys()))

    def showQuant(self):
        """Returns tuple (peptides DataFrame shape, proteins DataFrame shape).
        """
        return self._numbers


class Process(object):
    """Formerly omin.Experiment
    """
    def __init__(self, peptides_file, proteins_file, modifications=None):
        """
        """
        modifications = modifications or ["Acetyl", "Phospho"]
        self.raw_data = RawData(peptides_file, proteins_file)
        pep_sel, prot_sel = SelectionTools.vLook(self.raw_data.peptides,
                                                 self.raw_data.proteins,
                                                 modifications)
