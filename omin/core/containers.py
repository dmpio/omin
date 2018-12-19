# -*- coding: utf-8 -*-
"""
omin.core.containers
--------------------

Provides Container, Normalized, MaxQuantRaw, ProteomeDiscovererRaw,
PeptideGroups, and Proteins classes.

"""
# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

# ----------------
# EXTERNAL IMPORTS
# ----------------
import re
import os
import string
import itertools
import numpy as np

# Import pandas with the pandomics plug-in.
from pandomics import pandas as pd
# Import the guipyter DataLoader class.
from guipyter import DataLoader

# from pandomics import pandas

# ----------------
# INTERNAL IMPORTS
# ----------------
from .base import repr_dec, Handle

# -------------
# UTILS IMPORTS
# -------------
from ..utils import IOTools
from ..utils import StringTools
from ..utils import SelectionTools
from ..utils import IntermineTools

# --------
# DATBASES
# --------
from ..databases import MitoCartaTwo


# ================
# NORMALIZED CLASS
# ================
@repr_dec
class Normalized(object):
    """Empty class that catches results from normalization methods.
    """
    def __init__(self, *args, **kwargs):
        for k,v in kwargs.items():
            self.__dict__[k] = v


# ===============
# OCCUPANCY CLASS
# ===============
@repr_dec
class Occupancy(object):
    """Empty class that catches results from relative occupancy methods.

    NOTE: This is imported and applied by handles.Process.
    """
    def __init__(self, *args, **kwargs):
        for k,v in kwargs.items():
            self.__dict__[k] = v


# ===============
# CONTAINER CLASS
# ===============

class Container(DataLoader, Handle):
    """Base class for handling raw files.

    Accepted Raw files currently limited to

    """

    def __init__(self, *args, **kwargs):
        """Initialize base class for Proteome Discoverer raw files.

        This class inherits attributes from the DataLoader and Handle classes respectfully.

        Attributes
        ----------
        raw : pandas.DataFrame
            With pandomics plugin.
        file_path : str
        file_name : str
        file_ext : str
        metadata : dict

        """

        # Initialize the DataLoader class with args and kwargs
        DataLoader.__init__(self, *args, **kwargs)
        # Initialize the Handle class.
        Handle.__init__(self)
        self._title = "empty container"
        # METADATA: Collect basic file info
        self.metadata["file_name"] = self.file_name
        self.metadata["file_path"] = self.file_path
        self.metadata["file_ext"] = self.file_ext


    def show_metadata(self):
        """Print the metadata found in the container.
        """
        exculsions = set(["file_path", "file_ext"])
        header = " ".join([self._title, "Metadata"])
        print(header)
        print(len(header)*"=")
        for k,v in self.metadata.items():
            if k not in exculsions:
                print(k,":",v)


class MaxQuantRaw(Container):
    """Base class for MaxQuant raw files.

    This class inherits attributes from the Container class.

    HARD HAT AREA: Under construction.
    """
    def __init__(self, *args, **kwargs):
        # FIXME: Add functional algorithms.
        # FIXME: Add as mixin class later on.
        Container.__init__(self, *args, **kwargs)


class ProteomeDiscovererRaw(Container):
    """Base class for Proteome Discoverer raw files.

    This class inherits attributes from the Container class.
    """

    def __init__(self, *args, **kwargs):
        """Initialize base class for Proteome Discoverer raw files.
        """
        # Initialize the Container class
        Container.__init__(self, low_memory=False, delimiter='\t', **kwargs)
        # Calling the function below.
        self._expose_thermo_categories()
        # METADATA: Number of ids no filter.
        self.metadata['total_ids'] = self.raw.shape[0]
        # METADATA: Number of ids with quantification.
        self.metadata['total_ids_with_quant'] = self.Abundance.dropna(axis=0, how='all').shape[0]
        # Create an index of the ids with quantification
        self.quantified_index = self.Abundance.dropna(axis=0, how='all').index


    def _expose_thermo_categories(self):
        """Creates attributes that are pandas DataFrames based on Thermo's categories.
        """
        # FIXME: This method could potentially create collisions, figure out a way to safe gard against this.
        # fixme: Perhaps we could use one multi-index dataframe, in any case this loop is probably a substantial memory leak.
        thermo_category_values = set([i.split(":")[0] for i in self.raw.columns])
        thermo_category_keys = list(map(StringTools.remove_punctuation, thermo_category_values))
        thermo_category_keys = list(map(lambda x: x.strip().replace(" ", "_"), thermo_category_keys))
        thermo_category = dict(zip(thermo_category_keys, thermo_category_values))
        # This loop creates DataFrames for each category that Thermo has decided is important.
        for k,v in thermo_category.items():
            filter_with_colon = self.raw.filter(regex=v+":")
            if filter_with_colon.shape[-1] > 0:
                self.__dict__[k] = filter_with_colon
            else:
                filter_without_colon = self.raw.filter(regex=v)
                self.__dict__[k] = filter_without_colon

    # ------------------
    # STUDY FACTOR TOOLS
    # ------------------
    @property
    def study_factor_table(self):
        """Attempt to retrieve the study factors from Abundance column headers.

        Study Factors should be in the format:
            KO (Genotype), WT (Genotype), Input (Fraction), Phospho (Fraction) ect.

        If study factor cannot be identified then StudyFactor_<ABC...> will be assigned.
        """
        column_names = [i.split(":")[-1].strip() for i in self.Abundance.columns]

        column_names = pd.DataFrame(list(map(lambda x: x.split(','), column_names)))

        # Remove the sample column if present.
        if len(column_names.iloc[:, 1].unique()) < 2:
            column_names = pd.concat([column_names.iloc[:, 0], column_names.iloc[:, 2:]], axis=1)

        # Compile the regular expression to catch things in (parenthesis).
        # FIXME: BIG ASSUMPTION HERE. DOCUMENT OR DIE. #DOD
        rx = re.compile("\((.+)\)")
        study_factors = []
        for i in range(1, len(column_names.columns)):
            try:
                ucols = column_names.iloc[:, i].apply(lambda x: rx.findall(x)[0])
                ucols = ucols.unique()
                study_factor = ucols[0]
                study_factors.append(study_factor)

            except IndexError as err:
                # If study factor cannot be identified then StudyFactor_<ABC...> will be assigned.
                study_factor = "StudyFactor_" + string.ascii_uppercase[i-1]
                study_factors.append(study_factor)

        column_names.columns = ['_TMT_tag'] + study_factors
        # Prepend the fraction number.
        fn = [i.split(":")[1].strip() for i in self.Abundance.columns]
        fn = pd.DataFrame(fn, columns=["_Fn"])

        column_names = pd.concat([fn, column_names],axis=1)
        return column_names

    @property
    def study_factor_dict(self):
        """Return the study_factor_table as a flattened dict.
        """
        study_factor_dict = dict()
        for i in range(2, self.study_factor_table.columns.shape[0]):
            usf = self.study_factor_table.iloc[:, i].unique()
            usf = list(map(lambda x: x.strip(), usf))
            study_factor_dict[self.study_factor_table.iloc[:, i].name] = usf
        return study_factor_dict

    @property
    def tmt_plex_number(self):
        """Returns the TMT plex number.
        """
        plex_number = 0

        try:
            plex_number = self.study_factor_table._TMT_tag.unique().shape[0]

        except Exception as err:
            print(err)

        return plex_number


    @property
    def input_number(self):
        """Return number of inputs as a float.

        Parameters
        ----------
        dataframe : DataFrame

        Returns
        -------
        number_input : float
        """
        number_input = 0.0
        try:
            plex_number = self.tmt_plex_number
            # FIXME: Proof this against Inputase ect.
            # To do that that load the regex with ~ [Ii]nput[\w[punctuation]]
            # FIXME: Redundant filtering of raw DataFrames = memory leak
            input_abundance = self.Abundance.filter(regex="[Ii]nput").shape[1]
            number_input = input_abundance/plex_number
        except Exception:
            print("utils.SelectionTools.find_number_input failed")
        return number_input


    @property
    def study_factor_with_input(self):
        """Return the study factor that contains "Input" or "input".

        PROTIP: In your when creating study factors stick with the conventions:

            Input (Fraction), Acetyl (Fraction), ect.

        NOTE: For this to work correctly there must be only one study factor
        that contains the terms: "Input" or "input".

        """
        # FIXME: Ensure that the methods that use this know what to do with None.
        # FIXME: Inputase-proof this function or provide informative error handling with messages.
        rx = re.compile("[Ii]nput")
        study_factor_with_input = None
        for k,v in self.study_factor_dict.items():
            #FIXME: BIG ASSUMPTION HERE -> There will be only one study factor that contains a term including [Ii]nput.
            li = bool(sum(list(map(lambda x:bool(len(rx.findall(x))), v))))
            if li:
                study_factor_with_input = k
            else:
                pass
        return study_factor_with_input


    def fraction_tag(self, fraction_number=None):
        """Return a characteristic string for a fraction.

        Parameters
        ----------
        fraction_number: str
            Must be the form Fn e.g. F1, F2, ect.
        """
        # FIXME: Add try and excepts with some kind of unique string passed
        by_fn = self.study_factor_table.loc[self.study_factor_table._Fn == fraction_number].iloc[:, 2:]

        # Regex for isolating the study factor from it's catagory.
        rx = re.compile('(.+)\s\(.+\)')

        tag_for_fraction = ""
        for i in by_fn:
            terms = by_fn[i].unique()

            terms = list(map(lambda x:x.strip(), terms))
            #  FIXME: add some try and excepts with better docs.
            if len(terms) == 1:

                if any(["(" in i or ")" in i for i in terms]):
                    term = terms[0]
                    # Isolate the study factor from it's catagory.
                    term = rx.findall(term)
                else:
                    term = terms

                if len(term) == 1:
                    # Make everything lowercase.
                    term = term[0].lower()
                    # Replace any whitespace with underscore.
                    term = re.sub("\s", "_", term)
                    # print(term)
                    tag_for_fraction = "_".join([tag_for_fraction, term])

        return tag_for_fraction


    @property
    def fraction_number2fraction_tag(self):
        """Returns a dict with fraction numbers as keys and fraction tags as values.
        """
        result = dict([(i,self.fraction_tag(i)) for i in self.study_factor_table._Fn.unique()])
        return result


    def _gene_name_extractor(self):
        """
        """
        dataframe = self.master_index.copy()

        result = None

        if "Description" in dataframe:
            result = dataframe.Description.str.extract("GN=([A-Za-z0-9]+)")
            result = pd.concat([result, dataframe], axis=1)
            result.rename(columns={0:"Gene Name"}, inplace=True)
            self.master_index = result

        else:
            pass


# ===================
# PEPTIDEGROUPS CLASS
# ===================

class PeptideGroups(ProteomeDiscovererRaw):
    """Base class for Peptide Groups.

    Derived from the ProteomeDiscovererRaw class
    """

    def __init__(self, filepath_or_buffer=None, *args, **kwargs):
        """Initialize the base class."""
        # filepath_or_buffer = filepath_or_buffer or None
        ProteomeDiscovererRaw.__init__(self, filepath_or_buffer=filepath_or_buffer, title="Select peptide groups file", *args, **kwargs)
        # Set the title.
        self._title = "Peptide Groups"
        # Fill all of the modifications with NaNs with a blank string
        try:
            self.raw.Modifications.fillna('', inplace=True)
        except Exception as err:
            print("No Modifications column found in Peptide Groups data.", err)
        # Create the master_index
        self.master_index = None
        self._set_master_index()

        # Declare variable
        self._in_vivo_modifications = []
        self._set_in_vivo_modifications()


    def _set_master_index(self):
        """Attempt to set the master index for the Proteins.
        """
        # FIXME: Add try and except for each of these columns.
        master_index_components = ["Sequence", "Modifications", "Modifications in Proteins"]

        # Start the master index with the the first master protein accession.
        self.master_index = pd.DataFrame(self.raw["Master Protein Accessions"].first_member())
        # Change the name of "Master Protein Accessions" to "Accession".
        try:
            self.master_index.rename(columns={"Master Protein Accessions":"Accession"}, inplace=True)
        except Exception as err:
            pass

        for i in master_index_components:
            try:
                self.master_index = pd.concat([self.master_index, self.raw[i]],axis=1)
            except Exception as err:
                # print(err)
                pass
        return

    def _simplify_modifications(self):
        """Return a series of simplifed modifications.

        Parameters
        ----------
        df : DataFrame

        Returns
        -------
        simplified : Series
        """
        simplified = None
        if any(self.raw.columns == "Modifications"):
            # Compile the regular expression to be used in findall
            # rx = re.compile("x(\w+)\s")
            # simplified = dataframe.Modifications.apply(rx.findall)
            try:
                simplified = self.raw.Modifications.str.findall("x(\w+)\s")
                simplified = simplified.apply(lambda x: x if isinstance(x, list) else [])
                return simplified
            except Exception as err:
                print(err)
                return simplified
        else:
            print("No Modifications column found in Peptide Groups data.")


    def get_all_modifications(self):
        """Return set of ALL of TYPES of modifications found.

        Parameters
        ----------
        df : Dataframe

        Returns
        -------

        """

        found = set(itertools.chain.from_iterable(self._simplify_modifications()))

        return found


    def get_in_vivo_modifications(self):
        """Return list of the in vivo modifications.

        Parameters
        ----------
        df : Dataframe

        Returns
        -------
        invivo_modifications : set
        """

        present_modifications = self.get_all_modifications()
        chemical_modifications = {'Oxidation', 'Carbamidomethyl', 'TMT6plex', 'TMT10plex'}
        invivo_modifications = [modification for modification in present_modifications if modification not in chemical_modifications]
        return invivo_modifications

    def _simplify_modifications(self):
        """Return a series of simplifed modifications.

        Parameters
        ----------
        df : DataFrame

        Returns
        -------
        simplified : Series
        """
        simplified = None
        if any(self.raw.columns == "Modifications"):
            # Compile the regular expression to be used in findall
            # rx = re.compile("x(\w+)\s")
            # simplified = dataframe.Modifications.apply(rx.findall)
            try:
                simplified = self.raw.Modifications.str.findall("x(\w+)\s")
                simplified = simplified.apply(lambda x: x if isinstance(x, list) else [])
                return simplified
            except Exception as err:
                print(err)
                return simplified
        else:
            print("No Modifications column found in Peptide Groups data.")


    def get_all_modifications(self):
        """Return set of ALL of TYPES of modifications found.

        Parameters
        ----------
        df : Dataframe

        Returns
        -------

        """

        found = set(itertools.chain.from_iterable(self._simplify_modifications()))

        return found


    def get_in_vivo_modifications(self):
        """Return list of the in vivo modifications.

        Parameters
        ----------
        df : Dataframe

        Returns
        -------
        invivo_modifications : set
        """

        present_modifications = self.get_all_modifications()
        chemical_modifications = {'Oxidation', 'Carbamidomethyl', 'TMT6plex', 'TMT10plex'}
        invivo_modifications = [modification for modification in present_modifications if modification not in chemical_modifications]
        return invivo_modifications

    def _set_in_vivo_modifications(self):
        """FIXME: Add docs.
        """
        # FIXME: Rewirte for deprecarion of the omin.utils.SelectionTools.

        # Find invivo modifications.

        # in_vivo_mods = SelectionTools.findInVivoModifications(self.raw)
        in_vivo_mods = self.get_in_vivo_modifications()

        # Check if in_vivo modifications is None.
        if in_vivo_mods is not None:
            # If the list is greater than zero then set varible.
            if len(in_vivo_mods) > 0:
                self._in_vivo_modifications = in_vivo_mods
                # For each modification create a filtered dataframe attribute.
                # FIXME: This loop complates thing for a nuber of reasons it needs revision.
                for i in in_vivo_mods:
                    mod = StringTools.remove_punctuation(i).lower()
                    #FIXME: SelectionTools.filter_row will be migrated to pandomics.
                    # df = SelectionTools.filter_row(self.raw, on="Modifications", term=i)
                    df = self.raw.filter_rows(on= "Modifications", term=i)
                    # METADATA: Total Peptide Group IDs
                    self.metadata[mod+"_total_peptide_ids"] = df.shape[0]
                    self.__dict__[mod] = df
            else:
                pass
        else:
            pass


    def _separate_enriched_and_input(self):
        """Return object with load normalzed pd.DataFrames as attributes.
        """
        # FIXME: NORMALIZE INPUT TO ITS SELF.

        input_mask = self.study_factor_table[self.study_factor_with_input].str.contains("[Ii]nput")
        # FIXME: Replace the bobo method to find input fraction numbers
        number_input_fractions = len(self.study_factor_table.loc[input_mask]._Fn.unique())

        # isolate the input fractions study factors.
        inps = self.study_factor_table.loc[input_mask]
        inps = [inps.filter_rows(on="_Fn", term=i+"\Z") for i in inps._Fn.unique()]

        # isolate the other fractions study factors.
        frcs = self.study_factor_table.loc[~input_mask]
        frcs = [frcs.filter_rows(on="_Fn", term=i+"\Z") for i in frcs._Fn.unique()]

        return inps, frcs


    def _link_enriched_to_input(self):
        # Create a list of linkage DataFrames.
        inps, frcs = self._separate_enriched_and_input()

        linked = []
        for i in inps:
            for j in frcs:
                # If the fractions have the same shape then proceed with linking.
                if j.shape == i.shape:
                    # Generate a integer score of the similarity.
                    score = sum(sum(j.values == i.values))
                    # Create a DataFrame of the links.
                    link = pd.DataFrame(i.index, columns=["Link"], index=j.index)
                    score_df = pd.DataFrame([i.shape[0]*[score]]).T
                    score_df.columns = ["Score"]
                    score_df.index = j.index
                    j_prime = pd.concat([j, link, score_df], axis=1)
                    linked.append(j_prime)


        scores = np.array([i.Score.unique()[0] for i in linked])
        # # Isolate the highest scoring linked fractions
        # linked = list(filter(lambda x:x.Score.unique()[0] == scores.max(), linked))

        # # Create a cutoff list that is approximately half of the scores.
        # cut_off = scores[(-scores).argsort()][:int(len(scores)/2)]
        # linked = list(filter(lambda x:x.Score.unique()[0] in cut_off, linked))

        return linked


    def _set_linked_fractions(self):
        # linked = _link_enriched_to_input(self)
        linked = self._link_enriched_to_input()

        # Normalize the linked fractions
        self._linked_fractions = dict()
        for link in linked:
            # Get the input label.
            inp_label = self.fraction_tag(self.study_factor_table.iloc[link.Link]._Fn.unique()[0])

            other_label = self.fraction_tag(link._Fn.unique()[0])

            self._linked_fractions[other_label] = inp_label


    def load_normalized(self):
        """Normalize the inputs to themselves and add them to the normalized dict.
        NOTE: Inputs are normalized to themselves in order to calculate relative
        occupancy and looking into the whole proteome. Inputs that have been
        normalized to themselves are not used to normalize enriched fractions
        non-normalized inputs are.
        """

        self._set_linked_fractions()

        linked = self._link_enriched_to_input()
        inps, frcs = self._separate_enriched_and_input()


        normalized = dict()
        for inp in inps:
            inp_df = self.Abundance[self.Abundance.columns[inp.index]]
            inp_df = inp_df.normalize_to(inp_df)
            inp_label = self.fraction_tag(inp._Fn.unique()[0])
            normalized[inp_label] = inp_df

        for link in linked:
            inp = self.Abundance[self.Abundance.columns[link.Link]]
            other = self.Abundance[self.Abundance.columns[link.index]]
            other_label = self.fraction_tag(link._Fn.unique()[0])

            load_normalized = other.normalize_to(inp)
            normalized[other_label] = load_normalized

        # Throwing the normalized dict into the Normalized class
        return Normalized(**normalized)



    # FIXME: Define the two following functions at the PeptideGroups and Proteins Level.
    def std_out(self, selected_mod, comparisons):
        """
        Parameters
        ----------
        selected_mod: str

        comparisons: list
            In the format: [[numerator, denominator],...]
        """
        mod_mask = self.master_index.Modifications.str.contains(selected_mod)

        peptides_metadata = self.master_index.loc[mod_mask]

        peptides_metadata = pd.set_super_columns([peptides_metadata], ["Metadata"])

        results_list = []
        for k in self._linked_fractions.keys():
            norm_peptides = self.load_normalized.__dict__[k].loc[mod_mask]
            norm_peptides_log2 = norm_peptides.log2_normalize()

            comparison_labels = []
            relative_abundance_comparisons = []

            for i in comparisons:
                numerator = i[0]
                denominator = i[1]

                comparision_label = "Relative_Abundance"+k+"_Stats_{}_vs_{}".format(numerator, denominator)

                comparison_labels.append(comparision_label)

                relative_abundance_comp = norm_peptides_log2.fold_change_with_ttest(numerator=numerator,
                                                                                    denominator=denominator,
                                                                                    missing_values=True)
                relative_abundance_comparisons.append(relative_abundance_comp)

            rel_abun = pd.set_super_columns([norm_peptides, norm_peptides_log2],
                                            ["Load_Normalized"+k, "Load_Normalized_Log2_Normalized"+k,])

            rel_comp = pd.set_super_columns(relative_abundance_comparisons, comparison_labels)

            # Collect the related proteins. These are just for the readers that want to check the work.
            related_proteins = self.load_normalized_related_proteins.__dict__[self._linked_fractions[k]].loc[mod_mask]

            related_proteins_log2_normalized = related_proteins.log2_normalize()

            norm_peptides_occ = self.relative_occupancy.__dict__[k].loc[mod_mask]

            comparison_labels = []
            relative_occupancy_comparisons = []

            for i in comparisons:
                numerator = i[0]
                denominator = i[1]
                comparision_label = "Relative_Occupancy"+k+"_Stats_{}_vs_{}".format(numerator, denominator)
                comparison_labels.append(comparision_label)

                relative_occupancy_comp = norm_peptides_occ.fold_change_with_ttest(numerator=numerator,
                                                                                   denominator=denominator,
                                                                                   missing_values=True)
                relative_occupancy_comparisons.append(relative_occupancy_comp)

            rel_occ = pd.set_super_columns([norm_peptides_occ], ["Relative_Occupancy"+k])
            rel_occ = pd.set_super_columns([related_proteins, related_proteins_log2_normalized, norm_peptides_occ],
                                           ["Load_Normalized_Related_Protein_Abundances"+self._linked_fractions[k], "Log2_Normalized_Load_Normalized_Related_Protein_Abundances"+self._linked_fractions[k], "Relative_Occupancy"+k])

            rel_occ_comp = pd.set_super_columns(relative_occupancy_comparisons, comparison_labels)
            result = pd.concat([rel_abun, rel_comp, rel_occ, rel_occ_comp], axis=1)
            results_list.append(result)

        results_list.insert(0, peptides_metadata)

        result = pd.concat(results_list, axis=1)

        return result

# ==============
# PROTEINS CLASS
# ==============

class Proteins(ProteomeDiscovererRaw):
    """Base class for Proteins.
    Derived from the ProteomeDiscovererRaw class
    """

    def __init__(self, filepath_or_buffer=None, rescue_entrez_ids=False, *args, **kwargs):
        """

        Parameters
        ----------
        filepath_or_buffer: str or _io.TextWrapper

        rescue_entrez_ids: Bool
            Defaults to False.

        Attributes
        ----------
        metadata: dict

        master_high_confidence: pandas.DataFrame
            Filtered for high confidence master proteins (USE THIS ONE FOR COMPUTATION).

        """
        # filepath_or_buffer = filepath_or_buffer or None
        # FIXME: Add verbose argument.
        ProteomeDiscovererRaw.__init__(self, filepath_or_buffer=filepath_or_buffer, title="Select proteins file", *args, **kwargs)
        # Set the title.
        self._title = "Proteins"
        # ----------------------------------------------
        # Filter for master proteins and high confidence
        # ----------------------------------------------
        # FIXME: This may be redone in the future.
        self._high_confidence = self.high_confidence
        self.master_high_confidence = self.is_master_protein
        # METADATA: Number of high confidence protein.
        self.metadata['high_confidence_ids'] = self.master_high_confidence.shape[0]

        # Find and relabel the Entrez Gene ID column.
        self._set_entrez()
        # Select the first master protein accession.
        self._set_master_protein_accession()
        # Create the master_index
        self._set_master_index()
        # Attempt to rescue Entrez Gene IDs
        if rescue_entrez_ids:
            # ------------------------
            # IntermineTools Used here
            # ------------------------
            IntermineTools.rescue_entrez_ids(self.master_index)

        # self.gene_name_extractor()
        # Attach MitoCarta2 data to the master_index.
        self.add_database(MitoCartaTwo.essential)
        # Filter the abundance by master_high_confidence
        self.filter_abundance()

    @property
    def high_confidence(self):
        """
        Return a DataFrame of high confidence proteins.
        """
        result = self.raw
        try:
            try: # PD 2.1
                mask = self.raw["Exp. q-value"] < .01
                result = self.raw.loc[mask]
                return result
            except KeyError: # PD 2.2
                mask = self.raw['Protein FDR Confidence: Combined'] == "High"
                result = self.raw.loc[mask]
                return result
        except Exception as err:
            print("Could not filter for FDR")
            return result


    @property
    def is_master_protein(self):
        "Filter raw protein DataFrame for master proteins."



        try:
            # BIG ASSUMPTION: self.high_confidence is able to work
            # mask = self.high_confidence.Master == "IsMasterProtein"

            # For PD version compatibility.
            master_protein_terms = {"IsMasterProtein", "Master Protein"}
            mask = self.high_confidence.Master.apply(lambda x: x in master_protein_terms)

            result = self.high_confidence.loc[mask]
            return result

        except Exception as err:
            print('Could not Filter for master proteins.', err)
            return self.high_confidence


    @property
    def gene_symbols_from_description(self):
        """Return the gene symbols derived from the Fasta description.
        """
        result = self.master_high_confidence.Description.str.extract("GN=([a-zA-Z0-9]*)", expand=True).iloc[:, 0]
        return result


    def _set_entrez(self):
        """Find and relabel the Entrez Gene ID column.
        """
        if "Gene ID" in self.master_high_confidence: # PD2.1 Workaround
            #FIXME: In PD2.1 Gene ID is actually the gene name. Make sure that this will not break anything with PD2.2.
            self.master_high_confidence.rename(columns={'Gene ID':'GeneName'}, inplace=True)

        if "Entrez Gene ID" in self.master_high_confidence: # PD2.1 Workaround
            self.master_high_confidence.rename(columns={'Entrez Gene ID':'EntrezGeneID'}, inplace=True)


    def _set_master_protein_accession(self):
        """Find and relabel the Master Protein Accession column as Accession.
        """
        if "Accession" in self.master_high_confidence: # PD2.1
            pass

        if "Master Protein Accessions" in self.master_high_confidence: # PD2.2
            self.master_high_confidence.rename(columns={'Master Protein Accessions':'Accession'}, inplace=True)


    def _set_master_index(self):
        """Attempt to set the master index for the Proteins.
        """
        # FIXME: Add try and except for each of these columns.
        master_index_components = ["EntrezGeneID", "Description"]

        # Start the master index with the the first master protein accession.
        self.master_index = pd.DataFrame(self.master_high_confidence["Accession"].first_member())

        for i in master_index_components:
            try:
                result = pd.concat([self.master_index, self.master_high_confidence[i]], axis=1)
                result.index = self.master_index.index
                self.master_index = result
            except Exception as err:
                print(err)
                pass
        return


    def add_database(self, DataFrame):
        """Add databases to master_index.
        """
        # FIXME: Expose this function to users.
        try:
            # -------------------------------------------
            # DROPNA: Lose rows that have missing values.
            # -------------------------------------------
            self._old_master_index = self.master_index.copy()
            self.master_index.dropna(inplace=True)
            self.master_index.EntrezGeneID = self.master_index.EntrezGeneID.first_member().apply(np.int64)
            result = self.master_index.merge(DataFrame.copy(), on="EntrezGeneID", how="left")
            self.master_index = result
            self.master_index.fillna(False, inplace=True)

        except Exception as err:
            print(err)


    def filter_abundance(self):
        """Filter the Abundance columns by the the high confidence master proteins.
        """
        if "Abundance" in self.__dict__:
            try:
                self.Abundance = self.Abundance.loc[self.master_high_confidence.index]
            except Exception as err:
                print("Could not filter protein abundance by high confidence master proteins.", err)
        else:
            print("Could not filter protein abundance by high confidence master proteins.")


    def std_out(self, selected_mod, comparisons):
        """
        Parameters
        ----------
        selected_mod: str

        comparisons: list
            In the format: [[numerator, denominator],...]
        """

        proteins_metadata = self.master_index

        proteins_metadata =  pd.set_super_columns([proteins_metadata], ["Metadata"])

        results_list = []
        for k in self.load_normalized.__dict__.keys():

            norm_proteins = self.load_normalized.__dict__[k]
            norm_proteins_log2 = norm_proteins.log2_normalize()

            comparison_labels = []
            relative_abundance_comparisons = []

            for i in comparisons:
                numerator = i[0]
                denominator = i[1]

                comparision_label = "Protein_Expression_Stats"+k+"_{}_vs_{}".format(numerator, denominator)

                comparison_labels.append(comparision_label)

                relative_abundance_comp = norm_proteins_log2.fold_change_with_ttest(numerator=numerator,
                                                                                    denominator=denominator,
                                                                                    missing_values=True)
                relative_abundance_comparisons.append(relative_abundance_comp)

            rel_abun = pd.set_super_columns([norm_proteins, norm_proteins_log2],
                                            ["Load_Normalized"+k, "Load_Normalized_Log2_Normalized"+k,])

            rel_comp = pd.set_super_columns(relative_abundance_comparisons, comparison_labels)

            result = pd.concat([rel_abun, rel_comp], axis=1)

            results_list.append(result)

        results_list.insert(0, proteins_metadata)

        final_result = pd.concat(results_list, axis=1)

        return final_result
