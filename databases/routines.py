# -*- coding: utf-8 -*-

# FIXME: Add a "forget" function.
# FIXME: Make learn and forget fuzzy find.
# FIXME: Intergrate phrasewasher
# FIXME: Add transgenic conventions to genotype_terms.
# FIXME: Refactor for dynamic adding of dictionaries.
# FIXME: Force this to correllate to study factors.
# FIXME: See if we can make this ship with MACHINE LEARNING ALGO!!!

import pickle
import os

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
        self.modification_terms = pickle.load(open(this_dir+"\mod_dict.p", "rb"))
        self.genotype_terms = pickle.load(open(this_dir+"\\genotype_terms.p", "rb"))
        self.treatment_terms = pickle.load(open(this_dir+"\\treatment_terms.p", "rb"))

    def stop(self):
        """Save attributes as pickle files.
        """

        if type(self.modification_terms) is dict:
            pickle.dump(self.modification_terms, open(this_dir+"\mod_dict.p", "wb"))

        if type(self.genotype_terms) is dict:
            pickle.dump(self.genotype_terms, open(this_dir+"\\genotype_terms.p", "wb"))

        if type(self.treatment_terms) is dict:
            pickle.dump(self.treatment_terms, open(this_dir+"\\treatment_terms.p", "wb"))

    def learn(self, attribute, term, meaning):
        """
        Make any dict attribute learn a new term with given meaning.
        """
        try:
            self.__dict__[attribute][term] = meaning

        except Exception:
            print("I don't understand?")



if __name__ == "__main__":
    import re
    skynet = SkyNet()
    skynet.begin()
    # print(skynet.treatment_terms)
    # skynet.learn("treatment_terms", "10 Post", "ten_post")
    # print(skynet.treatment_terms)
    # print(re.findall(skynet.genotype_terms["ko"], "WT, Wt, wT, wt, KO, Ko, kO, ko"))
    # skynet.stop()
