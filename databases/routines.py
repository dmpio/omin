# -*- coding: utf-8 -*-

import pickle


class SkyNet(object):

    """Helping to reduce human error.
    """
    def __init__(self):
        """Define attributes.
        """
        self.modification_terms = None

    def begin(self):
        """Load attributes from pickled file sources.
        """
        self.modification_terms = pickle.load(open("databases\mod_dict.p", "rb"))

    def stop(self):
        """Save attributes as pickle files.
        """
        pickle.dump(self.modification_terms, open("databases\mod_dict.p", "wb"))

    def learn(self, term, meaning):
        self.modification_terms[term] = meaning


if __name__ == "__main__":
    skynet = SkyNet()
    skynet.begin()
    print(skynet.modification_terms)
