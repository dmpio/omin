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

# Load boilerplate modules.
import re
import pandas as pd
import numpy as np
import pickle

from omin.utils import SelectionTools

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
