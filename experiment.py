"""
experiment.py contains a number of functions that are designed to normalize proteomics data to an input fraction.
Author: James Draper
Email: james.drape@duke.edu
Date: October 12, 2016
"""

# Dict of modifications found in proteome discoverer modifications column
# FIXME: Make this a file in the module that can be updated by the user and automatically
# FIXME: make method that includes some amount of case handling.this may be better handled at manyModSel level.
mod_dict = {}
mod_dict["Acetyl"] = 'Acetyl'
mod_dict["Phospho"] = "Phospho"
mod_dict["HMG"] = "hydroxy...methyl.glutaryl"

#Load boilerplate modules.
import omin
import pandas as pd
import numpy as np

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

def manyModSel(pepdf, *terms):
    """Searches peptide DataFrame Modifications column for terms.

    Parameters
    ----------
        pepdf : DataFrame
            With peptides information.
        terms : str
            Can be any number of modifictions as a string. Case does not matter and regex special characters can be
            used e.g. 'acetyl', 'Phospho',hydroxy...methyl.glutaryl,'ect'
    Returns
    -------
        selected : DataFrame or np.NAN
            Contains the just the modified peptides.
    """
    for i in terms:
        moddex = pepdf.Modifications.str.contains(pat=i, case=False)
        if moddex.sum() > 0:
            selected = pepdf.ix[moddex]
            print(moddex.sum(), "peptides with", i, "modification found.")
            return selected
        else:
            print("No peptides with", i, "modification were found.")
            return np.nan

def normFactors(peptide_data):
    """Takes peptide abundance data and returns normalization factors.
    """
    return peptide_data.sum() / peptide_data.sum().mean()


def normalizeTo(different, normal):
    return different / normFactors(normal).as_matrix()


class Logger:
    """The class Logger preforms operations on normalized peptides.

    Attributes:
        log2 (:obj): DataFrame containing the log2 of the normalized data
        ave (:obj): DataFrame containing the average of the aformentioned DataFrame
        log_div_ave (:obj): DataFrame containing the log2-ave.
    """

    def __init__(self, normalized_data):
        """
        Args:
            normalized_data (:obj): Takes normalized peptide or protein in a DataFrame
        """
        logged = normalized_data.apply(np.log2)
        logged.columns = "Log2 " + logged.columns
        logged.index = normalized_data.index
        self.log2 = logged
        ave = self.log2.mean(axis=1)
        ave = pd.DataFrame(ave, columns=["Ave"])
        ave.index = normalized_data.index
        self.ave = ave
        log_div_ave = pd.DataFrame(self.log2.values - self.ave.values)
        log_div_ave.columns = "Log2-AVE " + normalized_data.columns
        log_div_ave.index = ave.index
        self.log_div_ave = log_div_ave


def masterCleanse(protein_df):
    """Filters raw protein DataFrame for master proteins.

    Args:
        protein_df (:obj): Raw protein DataFrame
    Return:
        clean (:obj): Protein DataFrame that contains only proteins with 'IsMasterProtein'
            in 'Master' column of protein_df
    """
    clean = protein_df.ix[protein_df.Master.str.endswith("IsMasterProtein")]
    return clean


def onePerQ(protein_df):
    """Filters raw protein DataFrame for proteins that are less than 1% the expected q-value.

    Args:
        protein_df (:obj): Raw protein DataFrame

    Return:
        clean (:obj): Protein DataFrame that contains only proteins with proteins only less
            than 1% of the expected q-value.
    """
    one_per = protein_df["Exp. q-value"] < .01
    one_per = protein_df.ix[one_per]
    return one_per


def masterOne(protein_df):
    """Takes a raw protein DataFrame and filters it using first the 'masterCleanse' function and then
    the 'onePerQ' function.

    Args:
        protein_df (:obj): Raw protein DataFrame

    Return:
        master_one (:obj): DataFrame of master proteins with exp. q-value <1%

    """
    master = masterCleanse(protein_df)
    master_one = onePerQ(master)
    return master_one


def masterPep(peptide_df):
    """Takes a peptide DataFrame and returns just the first master protein accession for each peptide.

    Note
    ----
    Peptides with no master protein accession will be lost however the index of peptide_df will be preseved.

    Parameters
    ----------
    peptide_df:

    Returns
    -------
    master_prot_acc:
    """
    master_prot_acc = [i.split(';')[0] for i in peptide_df['Master Protein Accessions'].dropna()]

    master_prot_acc = pd.DataFrame(master_prot_acc,
                                   index=peptide_df['Master Protein Accessions'].dropna().index, columns=['Accession'])
    return master_prot_acc

class WithInput:
    """

    """
    def __init__(self, raw):
        """
        Parameters
        ----------
        raw
        """
        self.raw = raw
        self.abundance = omin.sep(raw, 'Abundance:')
        self.inputs, self.enriched = omin.sepCon(self.abundance, 'Input')

class Compare:
    def __init__(self, abundance, genotype_list):
        '''

        Parameters
        ----------
        abundance
        genotype_list
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
    """
    This class detects the modified peptides and proteins for a given modification.
    """
    def __init__(self, whole_set, modification, genotypes, processed_peptides, processed_proteins, compare_in_order):
        # Grab the raw peptide data for modification in question
        self.raw = manyModSel(whole_set, modification)
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

class Experiment:
    """
    Takes instance of
    """
    def __init__(self, raw_file, modifications=["Acetyl", "Phospho"], genotypes=["ko", "wt"], compare_in_order=True):
        self.peptides = WithInput(raw_file.peptides)
        self.proteins = WithInput(raw_file.proteins)
        # FIXME: the next would probably be better handled as a outside fucntion or class idk
        # Normalize the enriched peptides to the input peptides
        setattr(self.peptides, 'norm', normalizeTo(self.peptides.enriched, self.peptides.inputs))
        # Normlize the enriched proteins to the input peptides
        setattr(self.proteins, 'norm', normalizeTo(self.proteins.inputs, self.peptides.inputs))
        # Run normalized peptides through Logger class
        setattr(self.peptides, 'norm_log', Logger(self.peptides.norm))
        # Run normalized proteins through logger class
        setattr(self.proteins, 'norm_log', Logger(self.proteins.norm))
        # Grab the master proteins with <1% FDR
        setattr(self.proteins, 'fdr', masterOne(raw_file.proteins))

        fdr_mpa_list = pd.DataFrame(self.proteins.fdr.Accession, index=self.proteins.fdr.Accession.index)

        # Dynamically set the modifications as class attributes
        for mod in modifications:
            if type(manyModSel(raw_file.peptides, mod_dict[mod])) != type(np.nan):
                self.__dict__[mod] = ModDetect(raw_file.peptides, mod_dict[mod], genotypes, self.peptides,
                                               self.proteins, compare_in_order)
