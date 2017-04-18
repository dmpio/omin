# -*- coding: utf-8 -*-
"""Omin core handles.

Handle in this context is a class composed of several pandas DataFrames, and
other varibles that are either derived from the DataFrames or provided by the
user.

Copyright 2017 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach,
Blair Chesnut, and Elizabeth Hauser.

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
from omin.utils import SelectionTools
from omin.normalize.toPool import NormalizedToPool
from omin.normalize.toInput import NormalizedToInput
from omin.databases import mitoCartaCall


# FIXME: Define varibles used at the highest possible class level.
# FIXME: Store each handle class as SQLite database in same parent dir.

class RawData(object):
    """
    Converts Proteome Discoverer .txt files into pandas DataFrames.

    Attributes
    ----------
    peptides :  DataFrame
        Raw data from Proteome Discoverer peptides data.
    proteins : DataFrame
        Raw data Proteome Discoverer corresponding proteins data.
    _numbers : tuple
    """

    def __init__(self, peptides_file, proteins_file):
        """Load data for peptides_file and proteins_file as pandas DataFrames.

        Note
        ----
        Please make sure that your files are in your current working directory.
        If you are working in jupyter notebook put copy your peptides groups
        and proteins files in the same directory as the notebook file.

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
        # FIXME: Add type checking.
        # FIXME: Add more try and excepts but try to put them on function level

        # Load your peptide groups file as a pandas DataFrame.
        self.raw_peptides = pd.read_csv(peptides_file,
                                        delimiter="\t",
                                        low_memory=False)
        # Load your protein file as a pandas DataFrame.
        self.raw_proteins = pd.read_csv(proteins_file,
                                        delimiter="\t",
                                        low_memory=False)
        # Store the shape of the respective DataFrames.
        self._numbers = (self.raw_peptides.shape, self.raw_proteins.shape)

    def __repr__(self):
        """Show all attributes."""
        return "Attributes: "+", ".join(list(self.__dict__.keys()))


class PreProcess(RawData):
    """A metaclass that uses RawData attempting several filtering steps.

    Attributes
    ----------
    _invivo_modifications : list
        A list of invivo modifications.
    pep_sel : DataFrame
        DataFrame to filter raw peptide groups DataFrame by indexing.
    prot_sel : DataFrame
        DataFrame to filter raw proteins DataFrame by indexing.
    mitodex : DataFrame
        For filtering mitochondrial peptide groups DataFrame by indexing.
    nonmitodex : DataFrame
        For filtering non-mitochondrial peptide groups DataFrame by indexing.

    Notes
    -----
    FIXME: Include paragraph description of the types of filtering.

    See Also
    --------
    omin.utils.SelectionTools.vLook
    omin.utils.SelectionTools.masterCleanse
    """

    """Handles Proteome Discoverer search results.
    """
    def __init__(self, peptides_file, proteins_file, modifications=None,
                 genotype=None, treatments=None):
        """Initalize instance of PreProcess class.

        Parameters
        ----------
        peptides_file : str
            The location of peptide groups file.
        proteins_file : str
            The location of peptide groups file.
        modifications : list
            List of derived or given modifications.
        genotype : list
            List of given genotypes. May be DEPRICATED in future.
        treatments: list
            List of given treatments. May be DEPRICATED in future.
        """
        # Initalize the RawData base class.
        super(PreProcess, self).__init__(peptides_file, proteins_file)
        # Find invivo modifications
        self._invivo_modifications = SelectionTools.findInVivoModifications(self.raw_peptides)
        # Set modifications with given or derived.
        modifications = modifications or self._invivo_modifications
        # Create 2 Dataframes that map specific peptide or protien uniprot ID
        # to it's relevent mitocarta index. Simillar to vlookup in excel.
        pep_sel, prot_sel = SelectionTools.vLook(self.raw_peptides,
                                                 self.raw_proteins,
                                                 modifications)

        # Used to index filter to statistically relevent PeptideGroups
        self.pep_sel = pep_sel
        # Used to index filter to statistically relevent Proteins
        self.prot_sel = prot_sel
        # MitoCarta calls made
        mito, nonmito = mitoCartaCall.mitoCartaPepOut(self,
                                                      mods=modifications,
                                                      dex=True)
        # Mitocarta hits DataFrame
        self.mitodex = mito
        # Mitocarta non hits DataFrame
        self.nonmitodex = nonmito
        # Find the number of inputs.
        self._input_number = SelectionTools.find_number_input(self.raw_peptides)
        # Find the number of PTMs present in the PeptideGroups
        self._ptm_fraction_numbers = SelectionTools.find_fractions(self.raw_peptides)


class Process(PreProcess):
    """A metaclass that uses PreProcess attempting several normalization steps.

    Attributes
    ----------
    _invivo_modifications : list
        A list of invivo modifications.
    pep_sel : DataFrame
        DataFrame to filter raw peptide groups DataFrame by indexing.
    prot_sel : DataFrame
        DataFrame to filter raw proteins DataFrame by indexing.
    mitodex : DataFrame
        For filtering mitochondrial peptide groups DataFrame by indexing.
    nonmitodex : DataFrame
        For filtering non-mitochondrial peptide groups DataFrame by indexing.

    Notes
    -----
    FIXME: Include paragraph description of the types of filtering.

    See Also
    --------
    omin.utils.SelectionTools.vLook
    omin.utils.SelectionTools.masterCleanse
    """

    def __init__(self, peptides_file, proteins_file, modifications=None,
                 genotype=None, treatments=None):
        """Initalize Process class.

        Parameters
        ----------
        peptides_file : str
            The location of peptide groups file.
        proteins_file : str
            The location of peptide groups file.
        modifications : list
            List of derived or given modifications.
        genotype : list
            List of given genotypes. May be DEPRICATED in future.
        treatments: list
            List of given treatments. May be DEPRICATED in future.
        """
        self.normalized = None
        super(Process, self).__init__(peptides_file, proteins_file)

        # FIXME: Make the selection more specifically target abundance columns
        if self._input_number > 0:
            inp_notify = "{} input fraction(s) found. Normalizing now..."
            inp_notify = inp_notify.format(self._input_number)
            print(inp_notify)
            try:
                self.normalized = NormalizedToInput(self.raw_peptides,
                                                    self.raw_proteins)
            except Exception:
                print("omin.normalize.toInput.NormalizedToInput FAILED.")

        elif self.raw_peptides.columns.str.contains("pool|control",
                                                    case=False).any():

            print("Pool columns. Omin will attempt to copare the data to it.")
            try:
                self.normalized = NormalizedToPool(self.raw_peptides)
            except Exception:
                print("omin.normalize.toPool.NormalizedToPool FAILED.")
        else:
            # FIXME: Make this a place where the user could specify.
            print("Cannot find anything to normalize to.")
