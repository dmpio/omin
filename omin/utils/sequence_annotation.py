# -*- coding: utf-8 -*-
# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

import re
import pandas as pd
import numpy as np


class SequenceAnnotationTools(object):
    """
    """

    @staticmethod
    def modification_formatter(mod_list):
        """Return regex term for a list of modifications."""
        mod_regex = "|".join(mod_list)

        mod_regex = "(?:{})".format(mod_regex)
        return mod_regex

    @classmethod
    def modifications_of_interest(cls, dataframe, mod_list):
        """Return Modifications of interest."""
        mod_regex = cls.modification_formatter(mod_list)

        rx = re.compile("\d.{}.+?\]".format(mod_regex))

        mod_interest = dataframe.Modifications.apply(rx.findall)
        # Format the as DataFrame
        mod_interest = pd.DataFrame(mod_interest,
                                    index=dataframe.index)

        mod_interest = mod_interest.Modifications.apply("".join)

        return mod_interest

    @classmethod
    def find_best_position(cls, dataframe, mod_list=None):
        """Find the best postion of modification."""
        if mod_list is None:
            try:
                mod_list = cls.findInVivoModifications(dataframe)
            except Exception:
                print("Please enter a mod_list and tray again.")
                return

        # Find the modifications of interest
        mod_interest = cls.modifications_of_interest(dataframe, mod_list)

        # Compile the regex to find best postions for modifications.
        # rx = re.compile("\w(\d+)\(\d+\)")
        rx = re.compile("\w(\d+)\(\d+.\d+\)")

        # Apply the compiled regex using findall to Modifications column.
        best_pos = mod_interest.apply(rx.findall)

        # Apply all of the find results
        best_pos = best_pos.apply(lambda x: np.array(x, dtype="int32"))

        return best_pos

    @classmethod
    def mods_in_proteins_cleaner(cls, dataframe):
        """Return DataFrame of cleaned Modifications in proteins."""
        # Compile the regex to find best postions for modifications.
        rx = re.compile(".+\d+\(\d+\)\]")

        cln = dataframe.filter(regex="Mod.+Proteins$").iloc[:,0].fillna("").apply(rx.findall).apply("".join)
        clnStr = lambda x:[i.strip() for i in x]
        brStr = lambda x:"<br>".join(clnStr(x))

        rx = re.compile(";")
        cln = cln.apply(rx.split).apply(brStr)
        cln = pd.DataFrame(cln)
        cln.columns = ["MIP_clean"]
        return cln

    @classmethod
    def htmlTagGen(cls, tag, style=None, mod=None):
        """Return generated HTML tags."""
        mod = mod or "style"

        if not tag.startswith("<"):
            tag = "<{}>".format(tag)

        close_tag = tag.replace("<", "</")

        if style is not None:
            mod = " {} =".format(mod)
            style = mod+"\'{}\'".format(style)

        else:
            style = ""

        style_tag = tag.replace(">", style+">")

        return style_tag, close_tag

    @classmethod
    def seqTagInject(cls, seq, ind, tag=None, style=None):
        """Return sequence with HTML tags injected."""
        tag = tag or "span"
        style_tag, close_tag = cls.htmlTagGen(tag, style)

        tags = style_tag+"{}"+close_tag

        seq_list = list(seq)
        for i in ind:
            seq_list[i] = tags.format(seq_list[i])
        return "".join(seq_list)

    @classmethod
    def tagInjector(cls, dataframe, best_pos=None, tag=None, style=None):
        """Return dataframe of html injected sequences."""
        style = style or "color:red;text-transform:lowercase;"

        if best_pos is None:
            best_pos = cls.find_best_position(dataframe)

        # seqTagInject(dataframe.Sequence[1],best_pos[1]-1,style=style)
        dataframe = dataframe.Sequence
        highlight_PTM = [cls.seqTagInject(dataframe[i], best_pos[i]-1, tag=tag, style=style) for i in dataframe.index]
        highlight_PTM_DF = pd.DataFrame(highlight_PTM,
                                        index=dataframe.index,
                                        columns=["TaggedSequence"])
        return highlight_PTM_DF
