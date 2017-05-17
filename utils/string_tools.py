"""String tools for omin."""

import string
import num2words
from datetime import datetime


class StringTools(object):

    @staticmethod
    def time_stamp():
        """Takes no arguments and returns datetime timestamp as string.

        Returns
        -------
        ts : str

        Examples
        --------
        >>>StingTools.time_stamp()
        '1483631005_38507'
        """
        ts = str(datetime.timestamp(datetime.now())).replace(".", "_")
        return ts

    @staticmethod
    def regexNot(term):
        """Return a string that has been formatted to negate a given term.

        Parameters
        ----------
        term : str or list

        Returns
        -------
        results : str
        """
        results = None
        if type(term) == str:
            results = "^(?!.*"+term+").*$"

        if type(term) == list:
            ored = "({})".format("|".join(term))
            results = "^(?!.*{}).*$".format(ored)

        return results

    @staticmethod
    def multiRegExOr(term_list):
        """Return a singles string of "OR"ed regex terms from a given list.

        Parameters
        ----------
        term_list: list
            List of regular expressions.

        Returns
        -------
        all_exp: str
            A string pattern ready for use with re methods.

        Exaples
        -------
        >>>re.search(multiRegExOr(skynet.modification_terms.values()),
        >>>                       'Phospho (fraction)')

        See Also
        --------
        re.search
        re.match
        """
        all_exp = "|".join(["(%s)" % (term) for term in term_list])
        return all_exp

    @staticmethod
    def multiRegExAnd(ex_list):
        """Return a regex psuedo "AND" formatted string from list of strings.

        Parameters
        ----------
        ex_list: list

        Returns
        -------
        formatted_list: str
        """
        andList = list(map(lambda x: "(?=.*\\b{}\\b)".format(x), ex_list))

        formatted_list = "^"+"".join(andList)+".*$"

        return formatted_list

    # @classmethod
    @staticmethod
    def phraseWasher(phrase, word_separator=" "):
        """Replace numerical portions of strings with words.

        Parameters
        ----------
        phrase : str
        word_separator : str
            Defaults to " " but you can put in whatever you like.

        See Also
        --------
        num2words.num2words

        Examples
        --------
        >>>phraseWasher("10 min post", word_separator = "_")
        "ten_min_post"
        >>>phraseWasher("65 min post", word_separator = "_")
        "sixty-five_min_post"

        """
        new_phrase = []
        for i in phrase.split():
            if i.isnumeric():
                phrase_part = num2words.num2words(int(i))
                new_phrase.append(phrase_part)
            else:
                new_phrase.append(i)
        washed = word_separator.join(new_phrase)
        return washed

    @staticmethod
    def remove_punctuation(start_str):
        """Remove punctuation from a given string."""
        table = str.maketrans({key: None for key in string.punctuation})
        non_punct = start_str.translate(table)
        return non_punct
