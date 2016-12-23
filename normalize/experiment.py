# -*- coding: utf-8 -*-
"""
experiment.py contains a number of functions that are designed to normalize
proteomics data to an input fraction.

FIXME : Beef up this documentation. Look for other sample python and R scripts
published we can use as template.

Author: James Draper
Email: james.drape@duke.edu
Date: October 12, 2016
"""

# FIXME : Remove deprecated classes below.

# Load boilerplate modules.

# import omin
import re
import pandas as pd
import numpy as np
import pickle

from utils import SelectionTools


class RawData:
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
        >>>raw_data = omin.RawData(peptides_file,proteins_file)

        """
        self.peptides = pd.read_csv(peptides_file,
                                    delimiter="\t",
                                    low_memory=False)

        print("number of peptides:", self.peptides.shape[0])

        self.proteins = pd.read_csv(proteins_file,
                                    delimiter="\t",
                                    low_memory=False)

        print("number of proteins:", self.proteins.shape[0])

    def __repr__(self):
        return "Attributes: "+", ".join(list(self.__dict__.keys()))

# === MAIN CLASS ===


class Experiment:
    """

    Attributes
    ----------
    peptides : (:obj)
        WithInput object.
    proteins : (:obj)
        WithInput object.
    [modifications] : (:obj)
        These modification could be anything defined by the user. The defaults
        are Acetyl and Phospho

    """
    def __init__(self,
                 raw_file=None,
                 modifications=["Acetyl", "Phospho"],
                 genotypes=["ko", "wt"],
                 compare_in_order=True,
                 treatments=None):

        """
        Parameters
        ----------
        raw_file : (:obj)
            Takes RawData object
        modifications : list
            These can be defined by the user. Defaults to ["Acetyl","Phospho"]
        genotypes : list
            These can be defined by the user. Defaults to ["ko","wt"]
        compare_in_order : bool
            Should the genotypes be compared in their current order with the
            first element as the numerator and the second as the denominator.
            If False the list is reversed. Defaults to True.

        Examples
        --------
        Load file strings:
        >>>peptides_file = "mydatafolder/peptides.txt"
        >>>proteins_file = "mydatafolder/proteins.txt"

        Create RawData object:
        >>>raw_data = omin.RawData(peptides_file,proteins_file)

        Create Experiment object
        >>>new_exp = omin.Experiment(raw_data)

        Create Experiment object with non-standard modifications, and
        treatments.
        >>>weird_exp = omin.Experiment(weird_data, modifications=["HMG"], treatments=["Rest", "10 Post"])
        """
        # Verify raw_file
        if not isinstance(raw_file, omin.experiment.RawData):
            raise NotImplementedError("omin.experiment.RawData object required")

        # NORMALIZE TO INPUT
        if raw_file.peptides.columns.str.contains("Input", case=False).any():
            print("Normalizing to: Input...")
            pep_sel,prot_sel = SelectionTools.vLook(raw_file.peptides,raw_file.proteins,modifications)
            self.peptides = omin.PeptidesWithInput(raw_file.peptides,modifications,pep_sel)
            self.proteins = omin.ProteinsWithInput(raw_file.proteins,modifications)
        else:
            # NORMALIZE TO POOL
            print("Normalizing to: Pool...")
            if treatments is None:
                print("Treatments were not specified. Please select them now.")

                treatments = omin.treatmentSelect(omin.sep(raw_file.peptides, "Abundance:"))

            pep_sel, prot_sel = omin.vLook(raw_file.peptides,
                                           raw_file.proteins,
                                           modifications)

            self.peptides = omin.WithPool(raw_file.peptides,
                                          modifications,
                                          genotypes,
                                          treatments)

            self.proteins = omin.WithPool(raw_file.proteins,
                                          modifications,
                                          genotypes,
                                          treatments)
    def __repr__(self):
        return "Attributes: "+", ".join(list(self.__dict__.keys()))
        # return "".join(list(self.__dict__.keys()))
