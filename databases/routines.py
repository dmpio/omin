# -*- coding: utf-8 -*-

# FIXME: Add a "forget" function.
# FIXME: Make learn and forget fuzzy find.
# FIXME: Intergrate phrasewasher
# FIXME: Add transgenic conventions to genotype_terms.
# FIXME: Refactor for dynamic adding of dictionaries.
# FIXME: Force this to correllate to study factors.
# FIXME: See if we can make this ship with MACHINE LEARNING ALGO!!!

import os
import re
import pickle
import numpy as np
from ..utils import SelectionTools
from ..utils import StringTools

this_dir, _ = os.path.split(__file__)

# print(this_dir)


class SkyNet(object):

    """Helping to reduce human error.
    """
    def __init__(self):
        """Define attributes.
        """
        self.modification_terms = None
        self.genotype_terms = None
        self.treatment_terms = None

    def begin(self):
        """Load attributes from pickled file sources.
        """
        self.modification_terms = pickle.load(open(this_dir+"\mod_dict.p",
                                                   "rb"))
        self.genotype_terms = pickle.load(open(this_dir+"\\genotype_terms.p",
                                               "rb"))
        self.treatment_terms = pickle.load(open(this_dir+"\\treatment_terms.p",
                                                "rb"))

    def stop(self):
        """Save attributes as pickle files.
        """

        if type(self.modification_terms) is dict:
            pickle.dump(self.modification_terms,
                        open(this_dir+"\mod_dict.p", "wb"))

        if type(self.genotype_terms) is dict:
            pickle.dump(self.genotype_terms,
                        open(this_dir+"\\genotype_terms.p", "wb"))

        if type(self.treatment_terms) is dict:
            pickle.dump(self.treatment_terms,
                        open(this_dir+"\\treatment_terms.p", "wb"))

    def learn(self, attribute, term, meaning):
        """Make any dict attribute learn a new term with given meaning.

        Parameters
        ----------
        attribute : str
        term : str
        meaning : str

        """
        try:
            self.__dict__[attribute][term] = meaning

        except Exception:
            print("I don't understand?")

    def train(self, attribute, term, new_def=None):
        """
        Parameters
        ----------
        attribute: str
        term: str
        new_def: str

        """
        if new_def is None:
            new_def = term

        rx = re.compile(term)
        newlist = list(filter(rx.findall, self.__dict__[attribute].keys()))
        [self.learn(attribute, i, new_def) for i in newlist]


def tagStudyFactors(abundance, ignore_phrases=None):
    """Returns a dict of tagged study factors.

    Notes
    -----
    FIXME : Make tagStudyFactors fail notification more robust.
    If column headers are not formated correctly however it says no abundance
    column is found.

    Parameters
    ----------
    abundance: DataFrame
        SPECIFICALLY peptide or protien abundance data.
        Meaning column headers that ALL begin with "Abundance:"

    ignore_phrases: str
        Can be a normal string or regex.

    Returns
    -------
    study_factors: dict
        With types of study factors as keys and different types of values.
    """

    ignore_phrases = ignore_phrases or "[Pp]ool"
    # If "control" needs to be ignored then it has to be added to the ignore phrases
    # ignore_phrases = ignore_phrases or "([Cc]ontrol)|([Pp]ool)"

    if not np.all(abundance.columns.str.contains("Abundance:")):
        print("One or more of these is not an abundance column.")
        return None

    geno_treat = [[j.strip() for j in i.split(",")] for i in abundance.columns if re.search(ignore_phrases,i) is None]

    # Create list of lists for study factors
    factors = []
    for i in range(len(geno_treat[0])):
        # Capture study factors as a list.
        factors.append(list(set([j[i] for j in geno_treat])))

    # Create a list of all the dicts in skynet
    skynet = SkyNet()
    skynet.begin()
    term_dict_list = [i for i in skynet.__dict__.keys() if i.endswith("_terms")]

    study_factors = dict()
    for fact in factors:
        # Filter out study factors that have 10 or more elements.
        if len(fact) < 10:
            # Filter out study factors that have one or fewer elements.
            if len(fact)>1:
                # In each of the dicts in skynet.
                for term_dict in term_dict_list:
                    # Get the term_dict name.
                    tl_name = term_dict.split("_")[0]
                    # Create the regex string of all of the terms in the term dict.
                    regex = StringTools.multiRegExOr(skynet.__dict__[term_dict].values())
                    # Test the study factor set for the presence of one of the terms.
                    if re.search(regex," ".join(fact)) is not None:
                        # If any of the terms is present add the study factor set to the dict.
                        study_factors[tl_name] = fact
                # Try to ID new study factors
                if fact not in study_factors.values() and re.findall("\([\w\s]+\)"," ".join(fact)) != None:
                    # Find the first occurence of the new study factor name.
                    new_factor_name = re.findall("\([\w\s]+\)", " ".join(fact))[0]
                    # Remove parentheses and spaces.
                    new_factor_name = re.sub("[\(\)]","",new_factor_name.replace(" ","_").lower())
                    # Add new study factor to the dict
                    study_factors[new_factor_name] = fact

    return study_factors

# if __name__ == "__main__":
#     import re
#     skynet = SkyNet()
#     skynet.begin()
    # print(skynet.treatment_terms)
    # skynet.learn("treatment_terms", "10 Post", "ten_post")
    # print(skynet.treatment_terms)
    # print(re.findall(skynet.genotype_terms["ko"], "WT, Wt, wt, KO, Ko, ko"))
    # skynet.stop()
