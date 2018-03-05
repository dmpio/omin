# -*- coding: utf-8 -*-
"""String tools for omin."""
# LICENSE
# -------

# Copyright 2017 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach,
# Blair Chesnut, and Elizabeth Hauser.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files, (the software)), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions: The above copyright
# notice and this permission notice shall be included in all copies or
# substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS",
# WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
# TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM. OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import string
from datetime import datetime


class StringTools(object):
    """Tools for processing strings or that have a string basis."""

    @staticmethod
    def time_stamp(posix=False):
        """Return datetime timestamp as string at the time called.

        Parameters
        ----------
        posix : bool
            If True returns the string in a POSIX format.

        Returns
        -------
        ts : str

        Examples
        --------
        >>>StingTools.time_stamp()
        '1483631005_38507'
        """
        if posix:
            ts = str(datetime.timestamp(datetime.now())).replace(".", "_")
        else:
            ts = "{:%I:%M:%S %p %A %B %d %Y}".format(datetime.now())
        return ts

    @staticmethod
    def remove_punctuation(start_str):
        """Remove punctuation from a given string."""
        table = str.maketrans({key: None for key in string.punctuation})
        non_punct = start_str.translate(table)
        return non_punct


    # @staticmethod
    # def regexNot(term):
    #     """Return a string that has been formatted to negate a given term.
    #
    #     Parameters
    #     ----------
    #     term : str or list
    #
    #     Returns
    #     -------
    #     results : str
    #     """
    #     results = None
    #     if type(term) == str:
    #         results = "^(?!.*"+term+").*$"
    #
    #     if type(term) == list:
    #         ored = "({})".format("|".join(term))
    #         results = "^(?!.*{}).*$".format(ored)
    #
    #     return results
    #
    # @staticmethod
    # def multiRegExOr(term_list):
    #     """Return a singles string of "OR"ed regex terms from a given list.
    #
    #     Parameters
    #     ----------
    #     term_list: list
    #         List of regular expressions.
    #
    #     Returns
    #     -------
    #     all_exp: str
    #         A string pattern ready for use with re methods.
    #
    #     Exaples
    #     -------
    #     >>>re.search(multiRegExOr(skynet.modification_terms.values()),
    #     >>>                       'Phospho (fraction)')
    #
    #     See Also
    #     --------
    #     re.search
    #     re.match
    #     """
    #     all_exp = "|".join(["(%s)" % (term) for term in term_list])
    #     return all_exp
    #
    # @staticmethod
    # def multiRegExAnd(ex_list):
    #     """Return a regex psuedo "AND" formatted string from list of strings.
    #
    #     Parameters
    #     ----------
    #     ex_list: list
    #
    #     Returns
    #     -------
    #     formatted_list: str
    #     """
    #     andList = list(map(lambda x: "(?=.*\\b{}\\b)".format(x), ex_list))
    #
    #     formatted_list = "^"+"".join(andList)+".*$"
    #
    #     return formatted_list

    # # @classmethod
    # @staticmethod
    # def phraseWasher(phrase, word_separator=" "):
    #     """Replace numerical portions of strings with words.
    #
    #     Parameters
    #     ----------
    #     phrase : str
    #     word_separator : str
    #         Defaults to " " but you can put in whatever you like.
    #
    #     See Also
    #     --------
    #     num2words.num2words
    #
    #     Examples
    #     --------
    #     >>>phraseWasher("10 min post", word_separator = "_")
    #     "ten_min_post"
    #     >>>phraseWasher("65 min post", word_separator = "_")
    #     "sixty-five_min_post"
    #
    #     """
    #     new_phrase = []
    #     for i in phrase.split():
    #         if i.isnumeric():
    #             phrase_part = num2words.num2words(int(i))
    #             new_phrase.append(phrase_part)
    #         else:
    #             new_phrase.append(i)
    #     washed = word_separator.join(new_phrase)
    #     return washed
