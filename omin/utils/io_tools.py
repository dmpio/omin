# -*- coding: utf-8 -*-
"""omin.io_tools

Provides
--------
Tools for file input output handling.
"""
# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

import os
import re


class IOTools(object):
    """Tools for file handling and dir building."""

    @staticmethod
    def sanitize_file_path(path):
        """Return an os safe file string.
        """
        result = re.sub('[^0-9a-zA-Z]+', '_', path)
        return result

    @classmethod
    def file_path_here(cls, path, ext=None):
        """Return a sanitized file path for the current working directory.
        """
        path = cls.sanitize_file_path(path)

        if ext is not None:
            path = '.'.join([path, ext])

        here = os.path.abspath('.')
        result = os.path.join(here, path)
        return result

    @staticmethod
    def mkdir(directory):
        """Create a directory.

        Parameters
        ----------
        directory : str
        """
        assert type(directory) == str
        # directory = StringTools.remove_punctuation(directory)
        directory = directory.replace(" ", "_")
        if not os.path.exists(directory):
            os.makedirs(directory)
            # print(directory, "is ready.")
            return directory
        else:
            # print(directory, "already exists.")
            return directory
