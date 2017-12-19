# -*- coding: utf-8 -*-
"""io tools"""
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

import os
import re
# from .object_tools import inspectObject


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
