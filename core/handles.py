# -*- coding: utf-8 -*-
import pandas as pd
from ..utils import SelectionTools
from ..normalize.toPool import NormalizedToPool
from ..normalize.toInput import NormalizedToInput


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


class Process(RawData):
    """Formerly omin.Experiment
    """
    def __init__(self, peptides_file, proteins_file, modifications=None,
                 genotype=None, treatments=None):
        """
        """

        modifications = modifications or ["Acetyl", "Phospho"]
        # Initalize the RawData base class.
        super(Process, self).__init__(peptides_file, proteins_file)
        # Create 2 Dataframes that map specific peptide or protien uniprot ID
        # to it's relevent mitocarta index. Simillar to vlookup in excel
        pep_sel, prot_sel = SelectionTools.vLook(self.raw_peptides,
                                                 self.raw_proteins,
                                                 modifications)
        self.pep_sel = pep_sel
        self.prot_sel = prot_sel

        # FIXME: Make the selection more specifically target abundance columns
        if self.raw_peptides.columns.str.contains("Input", case=False).any():

            print("Input fraction found. omin will attempt to normalize the data to it.")

            try:
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
            
