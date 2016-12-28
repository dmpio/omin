# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import re
from methods import *


class ModificationDelinate:

    def __init__(self, abundance, mod, notlist):
        self.abundance = omin.specSel(abundance, [mod], notlist)

    def addAttribute(self, attribute_name, attribute_data):

        self.__dict__[attribute_name] = attribute_data

    def __repr__(self):
        return "Attributes: "+", ".join(list(self.__dict__.keys()))


class Fractions:

    def __init__(self, modifications, abundance, pepsel):

        modifications.append("Input")
        for mod in modifications:
            notlist = list(np.array(modifications)[np.array([mod != i for i in modifications])])
            self.__dict__[mod] = omin.ModificationDelinate(abundance, mod, notlist)

        modifications.remove("Input")

        for mod in modifications:
            self.__dict__[mod].addAttribute("load_normalized", normalizeTo(self.__dict__[mod].abundance, self.Input.abundance))
            self.__dict__[mod].addAttribute("load_norm_log", Logger(self.__dict__[mod].load_normalized))
            self.__dict__[mod].addAttribute("filtered",self.__dict__[mod].load_norm_log.log_div_ave.ix[pepsel.index])

    def addAttribute(self, attribute_name, attribute_data):
        self.__dict__[attribute_name] = attribute_data

    def __repr__(self):
        return "Attributes: "+", ".join(list(self.__dict__.keys()))


class PeptidesWithInput:

    def __init__(self, raw, modifications, pepsel):
        self.raw_file = raw.copy()
        self.abundance = omin.sep(self.raw_file, 'Abundance:')
        self.fractions = Fractions(modifications, self.abundance, pepsel)

    def addAttribute(self, attribute_name, attribute_data):
        self.__dict__[attribute_name] = attribute_data

    def __repr__(self):
        return "Attributes: "+", ".join(list(self.__dict__.keys()))


class ProteinsWithInput:
    def __init__(self, raw, modifications):
        self.raw_file = raw
        self.abundance = omin.sep(raw, 'Abundance:')

    def addAttribute(self, attribute_name, attribute_data):
        self.__dict__[attribute_name] = attribute_data

    def __repr__(self):
        return "Attributes: "+", ".join(list(self.__dict__.keys()))
