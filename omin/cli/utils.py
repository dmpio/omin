# -*- coding: utf-8 -*-

# LICENSE
# -------

# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach,
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

import xlrd
# import pandas as pd
# from pandomics import pandas as pd
# from .core import pandomics as pd
from ..core.pandomics import pandas as pd

class CLITools(object):
    """Command Line Interface Tools.
    """
    @staticmethod
    def list_choice_cli(inp, message="Choose from the following:\n"):
        """Return the users selection from a list.

        Parameters
        ----------
        inp : list

        message : string

        Returns
        -------
        result : string
        """
        show = ["{}) {}\n".format(n, i) for n,i in enumerate(inp)]
        show = ''.join(show)
        show = message+show
        print(show)
        selection = int(input())
        result = inp[selection]
        return result

    @classmethod
    def read_excel_cli(cls, file_name):
        """Return a DataFrame from the selected sheet of an Excel file.

        Parameters
        ----------
        file_name : string

        Returns
        -------
        result : DataFrame

        """
        book = xlrd.open_workbook(file_name)
        if book.nsheets > 1:
            selected_sheet = cls.list_choice_cli(book.sheet_names())
            result = pd.read_excel(file_name, sheet_name=selected_sheet)
        else:
            result = pd.read_excel(file_name)
        return result
