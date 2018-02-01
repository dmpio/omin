# -*- coding: utf-8 -*-
"""
Copyright 2017 James Draper, Paul Grimsrud, Deborah Muoio

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files, Omics Modeling Integrating
Normalization (OMIN), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom
the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM.
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import re
import os
import pickle
import numpy as np
import pandas as pd
import guipyter
from pandas import ExcelWriter
import xlsxwriter
from xlsxwriter.utility import xl_col_to_name
from ..utils import objectWalker
from ..utils import StringTools


class StdOut(object):
    """omin's canonical output handler.

    Attributes
    ----------
    writer: :obj:
        A pandas.ExcelWriter object

    workbook: :obj:
        A xlsxwriter workbook object
    """
    def __init__(self,
                 the_object=None,
                 file_name=None):
        """
        Parameters
        ----------
        the_object: :obj:
            This can be almost any kind of object.

        file_name: str

        parent_file: str
        """

        file_name = file_name or guipyter.filedialog.asksaveasfilename()

        # Make sure the string has the file extension added on
        if not file_name.endswith(".xlsx"):
            file_name = ".".join([file_name, 'xlsx'])

        # Create writer.
        self.writer = pd.ExcelWriter(file_name, engine='xlsxwriter')
        # Create a workbook object.
        self.workbook = self.writer.book
        self.obj = the_object

    def makeXL(self, start_row=1):
        """Creates a .xlsx file from all of the DataFrames found in the object.

        Parameters
        ----------
        start_row: int
            Choose the starting row for all xlsx worksheets
        """

        workbook = self.workbook

        head_format = workbook.add_format({'font_name': 'Arial',
                                           'font_size': 10,
                                           'bold': True})
        head_format.set_align("justify")
        head_format.set_align("top")
        head_format.set_text_wrap()

        # pd.core.format.header_style = None
        pd.formats.format.header_style = None

        df_list = objectWalker(self.obj, "dataframe")

        for i in df_list:
            print("writing:",i[0])
            # Exclude dataframes with single column
            if i[-1].shape[1] > 1:
                i[-1].to_excel(self.writer,
                               sheet_name=i[0],
                               startrow=start_row,
                               index=False)

                worksheet = self.writer.sheets[i[0]]
                worksheet.set_row(start_row, 45, head_format)
                # The excel column number the corresponding to the df width
                end = xl_col_to_name(i[-1].shape[1])+"1"
                # Set excel columns format.
                worksheet.set_column('A1:'+end, 20)
        # Save the xlsx file
        self.writer.save()
        print("Your file is ready.")
        return
