# -*- coding: utf-8 -*-
"""
experiment.py contains a number of functions that are designed to normalize proteomics data to an input fraction.
FIXME : Beef up this documentation. Look for other sample python and R scripts published we can use as template.
Author: James Draper
Email: james.drape@duke.edu
Date: October 12, 2016
"""

# Dict of modifications found in proteome discoverer modifications column
# FIXME: Make normalization to input handle several modification types.
# FIXME: Make this a file in the module that can be updated by the user and automatically
# FIXME: make method that includes some amount of case handling.this may be better handled at manyModSel level.
#Load boilerplate modules.
import omin
import re
import pandas as pd
import numpy as np
import pickle
###LOAD MODIFICATION DICTIONARY###
#----------------------------------------------------------------------------------------------------------------------
# mod_dict = pickle.load(open("mod_dict.pickle","rb"))

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
        """
        self.peptides = pd.read_csv(peptides_file, delimiter="\t", low_memory=False)
        print("number of peptides:", self.peptides.shape[0])
        self.proteins = pd.read_csv(proteins_file, delimiter="\t", low_memory=False)
        print("number of proteins:", self.proteins.shape[0])


def masterCleanse(protein_df):
    """Filters raw protein DataFrame for master proteins.

    The raw protein data from Proteome Discoverer there is a column with the title 'Master' this funtion scans through
    that column and selects only the proteins that end with the string "IsMasterProtein"

    Parameters
    ----------
    protein_df : DataFrame
        Raw protein DataFrame

    Returns
    -------
    clean : DataFrame
        Protein DataFrame that contains only proteins with 'IsMasterProtein' in 'Master' column of protein_df
    """
    clean = protein_df.ix[protein_df.Master.str.endswith("IsMasterProtein")]
    return clean


def onePerQ(protein_df):
    """Filters raw protein DataFrame for proteins that are less than 1% the expected q-value.

    Scans through the protein DataFrame selecting only the proteins with less than 1% of the expected q-value.

    Parameters
    ----------
    protein_df : DataFrame
        Raw protein DataFrame

    Returns
    -------
    clean : DataFrame
        Protein data that contains only proteins with proteins only less than 1% of the expected q-value.
    """
    one_per = protein_df["Exp. q-value"] < .01
    one_per = protein_df.ix[one_per]
    return one_per


def masterOne(protein_df):
    """Takes a raw protein DataFrame and filters it using first the 'masterCleanse' function and 'onePerQ' function.

    Parameters
    ----------
    protein_df : DataFrame
        Raw proteins.

    Returns
    -------
    master_one : DataFrame
        Of master proteins with exp. q-value <1%
    """
    master = masterCleanse(protein_df)
    master_one = onePerQ(master)
    return master_one


def masterPep(peptide_df):
    """Takes a peptide DataFrame and returns just the first master protein accession for each peptide.

    Notes
    -----
    Assumes the first uniprot ID list is the correct one. Peptides with no master protein accession will be lost however
    the index of peptide_df will be preserved.

    Parameters
    ----------
    peptide_df : DataFrame

    Returns
    -------
    master_prot_acc : DataFrame

    """
    master_prot_acc = [i.split(';')[0] for i in peptide_df['Master Protein Accessions'].dropna()]

    master_prot_acc = pd.DataFrame(master_prot_acc,
                                   index=peptide_df['Master Protein Accessions'].dropna().index, columns=['Accession'])
    return master_prot_acc

class Compare:
    """
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
    """This class detects the modified peptides and proteins for a given modification.

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
        self.mpa = masterPep(self.raw)
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

class FracParse:
    """
    Attributes
    ----------
    abundance : DataFrame
    log_div_ave : DataFrame
    pool_normalized : DataFrame
    [selected conditions] : DataFrame
        These are created dynamically from the list the user selects.
    """
    def __init__(self,abundance,select_list):
        """
        Parameters
        ----------
        abundance : DataFrame
        select_list : list
        """
        self.abundance = abundance
        self.log_div_ave = omin.logNormToAve(self.abundance)
        self.pool_normalized = omin.normToPool(self.log_div_ave)

        for select in select_list:
            term = re.sub("_", " ", select)
            self.__dict__[select] = omin.sep(self.pool_normalized,term)

class PoolMod:
    """

    Attributes
    ----------
    abundance : DataFrame
    [genotype] : (:obj)
        FracParse object.
    """
    def __init__(self,abundance,mod,genotypes,select_list):
        self.abundance = omin.sep(abundance,mod)
        for geno in genotypes:
            self.__dict__[geno] = FracParse(omin.sep(self.abundance,geno),select_list)

class WithPool:
    """
    Attributes
    ----------
    raw : DataFrame
        The raw DataFrame with all information.
    abundance : DataFrame
        Abundance columns from raw DataFrame.
    """

    def __init__(self, raw,modifications,genotypes,select_list):
        """

        Parameters
        ----------
        raw: DataFrame

        """
        self.raw = raw
        self.abundance = omin.sep(raw, "Abundance:")

        for mod in modifications:
            self.__dict__[omin.mod_dict[mod]] = PoolMod(self.abundance,mod,genotypes,select_list)

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
    def __init__(self, raw_file, modifications=["Acetyl", "Phospho"], genotypes=["ko", "wt"], compare_in_order=True):
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
        """
        #NORMALIZE TO INPUT
        if raw_file.peptides.columns.str.contains("Input", case=False).any():
            print("Normalizing to: Input...")
            pep_sel,prot_sel = omin.vLook(raw_file.peptides,raw_file.proteins,modifications)
            self.peptides = omin.PeptidesWithInput(raw_file.peptides,modifications)
            #self.proteins = WithInput(raw_file.proteins,modifications)
            # # FIXME: the next couple of lines would probably be better handled in an outside function or class idk.
            # # Normalize the enriched peptides to the input peptides
            # setattr(self.peptides, 'norm', normalizeTo(self.peptides.enriched, self.peptides.inputs))
            # # Normlize the enriched proteins to the input peptides
            # setattr(self.proteins, 'norm', normalizeTo(self.proteins.inputs, self.peptides.inputs))
            # # Run normalized peptides through Logger class
            # setattr(self.peptides, 'norm_log', Logger(self.peptides.norm))
            # # Run normalized proteins through logger class
            # setattr(self.proteins, 'norm_log', Logger(self.proteins.norm))
            # # Grab the master proteins with <1% FDR
            # setattr(self.proteins, 'fdr', masterOne(raw_file.proteins))
            # fdr_mpa_list = pd.DataFrame(self.proteins.fdr.Accession, index=self.proteins.fdr.Accession.index)
            # # Dynamically set the modifications as class attributes
            # for mod in modifications:
            #     if len(manyModSel(raw_file.peptides, mod_dict[mod])) != 0:
            #         self.__dict__[mod] = ModDetect(raw_file.peptides, mod_dict[mod], genotypes, self.peptides,
            #                                        self.proteins, compare_in_order)
        else:
            print("Normalizing to: Pool...")
            select_list = conditionSelect(omin.sep(raw_file.peptides,"Abundance:"))
            self.peptides = WithPool(raw_file.peptides,modifications,genotypes,select_list)
            self.proteins = WithPool(raw_file.proteins,modifications,genotypes,select_list)
    def __repr__(self):
        return "Attributes: "+", ".join(list(self.__dict__.keys()))
        # return "".join(list(self.__dict__.keys()))

def conditionSelect(peptide_abundance):
    """Allows the user to select which columns contain conditions.

    Parameters
    ----------
    peptide_abundance : DataFrame

    Returns
    -------
    select_list : list
    """
    print(pd.DataFrame([i.split(",") for i in peptide_abundance.columns]))
    col_num = int(input("Which column has condition data?(Enter the number)"))
    select_set = set(pd.DataFrame([i.split(",") for i in peptide_abundance.columns]).ix[:, col_num])

    select_list = [re.sub(" ", "_", i.strip()) for i in select_set]
    #print(select_list)
    [print(n, i) for n, i in enumerate(select_list)]
    remove_num = int(input("Enter the number of any element that need to be removed."))
    select_list.remove(select_list[remove_num])
    return select_list
