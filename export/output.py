# -*- coding: utf-8 -*-
import re
import os
import pickle
import numpy as np
import pandas as pd
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
                 file_name=None,
                 parent_file=None):
        """
        Parameters
        ----------
        the_object: :obj:
            This can be almost any kind of object.

        file_name: str

        parent_file: str
        """

        # Get the current working directory
        cwd = os.getcwd()

        # If no parent file is given then current working directory is used.
        parent_file = parent_file or cwd

        # If No file name is give then timestamp is used.
        file_name = file_name or StringTools.time_stamp()

        # Make sure the string has the file extension added on
        if not file_name.endswith(".xlsx"):
            file_name = ".".join([file_name, 'xlsx'])

        # Create save location file path variable
        full_path = "\\".join([parent_file, file_name])

        # Create writer.
        self.writer = pd.ExcelWriter(full_path, engine='xlsxwriter')
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
