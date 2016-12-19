# -*- coding: utf-8 -*-
import re
from urllib.request import urlopen
from urllib.error import HTTPError

# FIXME: Create tools to search against against databases the user specifies.

# === UNIPROT TOOLS ===


class FastaTools(object):
    """A collection of tools for handling FASTA formatted strings.
    """

    @classmethod
    def fasta2Seq(cls, fasta):

        """ Joins strings in list into one string. Neglecting the first which is
        generally the annotation.

        Parameters
        ----------
        fasta : str

        Returns
        -------
        one_string : str
        """
        one_string = "".join(fasta[1:])
        return one_string

    @classmethod
    def seqNum(cls, seq):
        """Takes a sequence and returns the number for each amino acid (I think).

        Notes
        -----
            Review this functions useage and rewrite annotation.

        Parameters
        ----------
        seq : str

        Returns
        -------
        seq_num : str

        """
        seq_num = [seq[i]+str(i+1) for i in range(len(seq))]
        return seq_num

# === UNIPROT TOOLS ===

class UniprotTools(object):

    @classmethod
    def getFasta(cls, upID):
        """Takes a uniprot ID and returns a string.

        Parameters
        ----------
        upID : str

        Returns
        -------
        fasta : str

        """
        try:
            fasta = str
            link = "http://www.uniprot.org/uniprot/"+upID+".fasta"
            try:
                urlout = urlopen(link)
                # print("hey")

            except HTTPError:
                print(
                      "Try again, the formatting may be off.",
                      "Try removing dashes."
                      )

            # fasta = urlout.read().decode("utf-8")
            fasta = urlout.read().decode("utf-8").split("\n")
            return fasta
        except TypeError as err:
            print("Please enter UniProt ID as string.")

# === STRING TOOLS ====


class StringTools(object):

    # @classmethod
    @staticmethod
    def int2word(num, separator="-"):
        """Transforms integers =< 999 into english words

        Parameters
        ----------
        num : int
        separator : str

        Returns
        -------
        words : str
        """
        ones_and_teens = {0: "Zero", 1: 'One', 2: 'Two', 3: 'Three',
                          4: 'Four', 5: 'Five', 6: 'Six', 7: 'Seven',
                          8: 'Eight', 9: 'Nine', 10: 'Ten', 11: 'Eleven',
                          12: 'Twelve', 13: 'Thirteen', 14: 'Fourteen',
                          15: 'Fifteen', 16: 'Sixteen', 17: 'Seventeen',
                          18: 'Eighteen', 19: 'Nineteen'}
        twenty2ninety = {2: 'Twenty', 3: 'Thirty', 4: 'Forty', 5: 'Fifty',
                         6: 'Sixty', 7: 'Seventy', 8: 'Eighty', 9: 'Ninety', 0: ""}

        if 0 <= num < 19:
            return ones_and_teens[num]
        elif 20 <= num <= 99:
            tens, below_ten = divmod(num, 10)
            if below_ten > 0:
                words = twenty2ninety[tens] + separator + \
                    ones_and_teens[below_ten].lower()
            else:
                words = twenty2ninety[tens]
            return words

        elif 100 <= num <= 999:
            hundreds, below_hundred = divmod(num, 100)
            tens, below_ten = divmod(below_hundred, 10)
            if below_hundred == 0:
                words = ones_and_teens[hundreds] + separator + "hundred"
            elif below_ten == 0:
                words = ones_and_teens[hundreds] + separator + \
                    "hundred" + separator + twenty2ninety[tens].lower()
            else:
                if tens > 0:
                    words = ones_and_teens[hundreds] + separator + "hundred" + separator + twenty2ninety[
                        tens].lower() + separator + ones_and_teens[below_ten].lower()
                else:
                    words = ones_and_teens[
                        hundreds] + separator + "hundred" + separator + ones_and_teens[below_ten].lower()
            return words

        else:
            print("num out of range")


    @classmethod
    def phraseWasher(cls, phrase, number_separator="_", word_separator=" "):
        """Replaces numerical portions of strings with words.

        Parameters
        ----------
        phrase : str
        number_separator : str
            separator to be used with int2word function.
        word_separator : str
            Defaults to " " but you can put in whatever you like.

        See Also
        --------
        int2word

        Examples
        --------
        >>>phraseWasher("10 min post",word_separator = "_")
        "Ten_min_post"
        >>>phraseWasher("65 min post",number_separator = "_", word_separator = "_")
        "Sixty_five_min_post"

        """
        new_phrase = []
        for i in phrase.split():
            if i.isnumeric():
                phrase_part = cls.int2word(int(i), number_separator)
                new_phrase.append(phrase_part)
            else:
                phrase_part = i
                new_phrase.append(phrase_part)
        washed = word_separator.join(new_phrase)
        return washed
