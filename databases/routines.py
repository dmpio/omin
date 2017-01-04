# -*- coding: utf-8 -*-

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
        self.__dict__[attribute][term] = meaning


if __name__ == "__main__":
    skynet = SkyNet()
    skynet.begin()
    # print(skynet.treatment_terms)
    # skynet.learn("treatment_terms", "10 Post", "ten_post")
    # print(skynet.treatment_terms)
    skynet.stop()
