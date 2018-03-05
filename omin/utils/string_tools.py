# -*- coding: utf-8 -*-
"""omin.utils.string_tools

Provides
--------
String tools for omin.
"""
# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

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
