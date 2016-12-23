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

import omin
import re
import pandas as pd
import numpy as np
import pickle


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
        """Loads the raw data for peptides_file and proteins_file as DataFrames that are contained as attributes.

        Note
        ----
            Please make sure that your files are in your current working directory. If you are
            working in jupyter notebook please put the a copy of the peptides and proteins in
            the same directory as the notebook file.

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
        self.peptides = pd.read_csv(peptides_file, delimiter="\t", low_memory=False)
        print("number of peptides:", self.peptides.shape[0])
        self.proteins = pd.read_csv(proteins_file, delimiter="\t", low_memory=False)
        print("number of proteins:", self.proteins.shape[0])
    def __repr__(self):
        return "Attributes: "+", ".join(list(self.__dict__.keys()))

class Compare:
    """DEPRECATED
    Attributes
    ----------
    abundance : DataFrame
    lfc : DataFrame
        Log2 Fold Change in the given comparison.
    pval : DataFrame
        P-values in the given comparison.
    num_std : DataFrame
        The standard deviation of the numerator in the comparison.
    dem_std : DataFrame
        The standard deviation of the denominator in the comparison.
    num_ave : DataFrame
        The mean of the columns in the numerator of the comparison.
    dem_ave : DataFrame
        The mean of the columns in the denominator of the comparison.

    """
    def __init__(self, abundance, genotype_list):
        '''

        Parameters
        ----------
        abundance : DataFrame
            Abundance columns.
        genotype_list : list
            A list of the genotypes to be compared. The first element of the list is the numerator the second is the
            denominator in the comparison.

        '''
        self.abundance = abundance
        # Set the numerator
        num = omin.sep(abundance, genotype_list[0])
        # Set the denominator
        dem = omin.sep(abundance, genotype_list[1])
        # Get Log fold change
        self.lfc = omin.logFolder(num, dem)
        # Get the pvalue
        self.pval = omin.pvaln(num, dem)
        # Get the numerator standard deviation
        self.num_std = pd.DataFrame(num.std(axis=1), columns=[genotype_list[0] + " STDEV"], index=num.index)
        # Get the denominator standard deviation
        self.dem_std = pd.DataFrame(dem.std(axis=1), columns=[genotype_list[1] + " STDEV"], index=dem.index)
        # Get the numerator average
        self.num_ave = pd.DataFrame(num.mean(axis=1), columns=[genotype_list[0] + " AVE"], index=num.index)
        # Get the denominator average
        self.dem_ave = pd.DataFrame(dem.mean(axis=1), columns=[genotype_list[1] + " AVE"], index=dem.index)

class ModDetect:
    """DEPRECATED
    This class detects the modified peptides and proteins for a given modification.

    Attributes
    ----------
    raw : DataFrame
    mpa : DataFrame
    fdr : DataFrame
    fdr_mpa_prot : DataFrame
        FDR filterd master protein accessions with the same index as the proteins DataFrame.
    fdr_mpa_pep : DataFrame
        FDR filterd master protein accessions with the same index as the peptides DataFrame.
    correlated_protein_abundance : DataFrame
    carta_prot : DataFrame
    carta_pep : DataFrame
    bdex : DataFrame
        Just mitochondrial peptides.
    wdex : DataFrame
        Just the Non-mitochondrial proteins.
    norm_abundance : DataFrame
    relative_abundance : DataFrame
    relative_occupancy : DataFrame
    """
    def __init__(self, whole_set, modification, genotypes, processed_peptides, processed_proteins, compare_in_order):
        """

        Parameters
        ----------
        whole_set : DataFrame
            Set of all the peptides.
        modification : str
            The modification to detect.
        genotypes : list
            List of genotypes
        processed_peptides : DataFrame
            Normalized peptides.
        processed_proteins : DataFrame
            Normalized proteins.
        compare_in_order : bool
            True means the first element in the genotypes list the numerator in the comparison.

        """
        # Grab the raw peptide data for modification in question
        self.raw = manyModSel(whole_set, modification)[0]
        # Grab the master protein accession for the peptides
        self.mpa = omin.masterPep(self.raw)
        # Grab the protein accessions for the correlated protein abundances for the modifidied peptides
        fdr = pd.DataFrame(processed_proteins.fdr.Accession, index=processed_proteins.fdr.Accession.index)
        # Essentially the same as PAG's vlookup step
        fdr_mpa_prot = self.mpa.merge(fdr, on="Accession", how="left", left_index=True)
        fdr_mpa_pep = self.mpa.merge(fdr, on="Accession", how="left", right_index=True)
        # Grab the correlated protein abundances for the modifidied peptides
        self.correlated_protein_abundance = processed_proteins.norm_log.log_div_ave.ix[fdr_mpa_prot.index]
        #Make mitocarta 2.0 calls
        self.carta_prot = omin.mitoCartaCall.mitoProt(processed_proteins.raw)
        self.carta_pep = self.mpa.copy().merge(self.carta_prot.copy(), on="Accession", how="left")
        self.carta_pep.index = self.mpa.index

        self.bdex = self.carta_pep[self.carta_pep.MitoCarta2_List == 1.0]
        self.wdex = self.carta_pep[self.carta_pep.MitoCarta2_List != 1.0]

        # RELATIVE ABUNDANCE
        self.norm_abundance = processed_peptides.norm_log.log_div_ave.ix[fdr_mpa_pep.index]
        # Switch the order of the genotypes list if True
        if compare_in_order == True:
            pass
        else:
            genotypes = genotypes[::-1]
        #Set the relative abundance level
        self.relative_abundance = Compare(self.norm_abundance, genotypes)
        #RELATIVE OCCUPANCY
        occupancy_abundance = pd.DataFrame(self.norm_abundance.values - self.correlated_protein_abundance.values,
                                           index=fdr_mpa_pep.index)
        occupancy_abundance.columns = "Relative occupancy " + self.correlated_protein_abundance.columns
        self.relative_occupancy = Compare(occupancy_abundance, genotypes)

###MAIN CLASS###
#----------------------------------------------------------------------------------------------------------------------
class Experiment:
    """

    Attributes
    ----------
    peptides : (:obj)
        WithInput object.
    proteins : (:obj)
        WithInput object.
    [modifications] : (:obj)
        These modification could be anything defined by the user. The defaults are Acetyl and Phospho

    """
    def __init__(self, raw_file = None, modifications = ["Acetyl", "Phospho"], genotypes = ["ko", "wt"],
        compare_in_order = True, treatments = None):

        """

        Parameters
        ----------
        raw_file : (:obj)
            Takes RawData object
        modifications : list
            These can be defined by the user. If no list is entered then it defaults to ["Acetyl","Phospho"]
        genotypes : list
            These can be defined by the user. If no list is entered then it defaults to ["ko","wt"]
        compare_in_order : bool
            Should the genotypes be compared in their current order with the first element as the numerator and the
            second as the denominator. If False the list is reversed. Defaults to True.

        Examples
        --------
        Load file strings:
        >>>peptides_file = "mydatafolder/peptides.txt"
        >>>proteins_file = "mydatafolder/proteins.txt"

        Create RawData object:
        >>>raw_data = omin.RawData(peptides_file,proteins_file)

        Create Experiment object
        >>>new_exp = omin.Experiment(raw_data)

        Create Experiment object with non-standard modifications, and treatments.
        >>>weird_exp = omin.Experiment(weird_data, modifications=["HMG"], treatments=["Rest", "10 Post"])
        """
        #verify raw_file
        if not isinstance(raw_file,omin.experiment.RawData):
            raise NotImplementedError("raw_file must be a omin.experiment.RawData object.")
        #NORMALIZE TO INPUT
        if raw_file.peptides.columns.str.contains("Input", case=False).any():
            print("Normalizing to: Input...")
            pep_sel,prot_sel = omin.vLook(raw_file.peptides,raw_file.proteins,modifications)
            self.peptides = omin.PeptidesWithInput(raw_file.peptides,modifications,pep_sel)
            self.proteins = omin.ProteinsWithInput(raw_file.proteins,modifications)
        else:
        #NORMALIZE TO POOL
            print("Normalizing to: Pool...")
            if treatments == None:
                print("You didn't specify any treatments. Please select them now...")
                treatments = omin.treatmentSelect(omin.sep(raw_file.peptides,"Abundance:"))
            pep_sel,prot_sel = omin.vLook(raw_file.peptides,raw_file.proteins,modifications)
            self.peptides = omin.WithPool(raw_file.peptides,modifications,genotypes,treatments)
            self.proteins = omin.WithPool(raw_file.proteins,modifications,genotypes,treatments)
    def __repr__(self):
        return "Attributes: "+", ".join(list(self.__dict__.keys()))
        # return "".join(list(self.__dict__.keys()))
