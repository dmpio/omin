# -*- coding: utf-8 -*-
import re
import os
import pickle
import numpy as np
import pandas as pd
from pandas import ExcelWriter
import xlsxwriter
from xlsxwriter.utility import xl_col_to_name
# from datetime import datetime

from ..utils import StringTools


class Standard(object):

    @staticmethod
    def excel(dataframe,
              file_name=None,
              parent_file=None,
              sheet_name=None,
              start_row=3):
        """

        Parameters
        ----------

        Returns
        -------
        excel file
        """
        # Get the current working directory
        cwd = os.getcwd()

        # If no parent file is given then current working directory is used.
        parent_file = parent_file or cwd

        # If No file name is give then timestamp is used.
        file_name = file_name or StringTools.time_stamp()

        sheet_name = sheet_name or "data"

        # Make sure the string has the file extension added on
        if not file_name.endswith(".xlsx"):
            file_name = ".".join([file_name, 'xlsx'])

        full_path = "\\".join([parent_file, file_name])

        # Turn off header styles.
        pd.core.format.header_style = None

        # Create writer.
        writer = ExcelWriter(full_path,
                             engine='xlsxwriter',
                             options={'nan_inf_to_errors': True})

        # Use to_excel method to use writer on dataframe.
        dataframe.to_excel(writer,
                           sheet_name=sheet_name,
                           startrow=start_row,
                           index=False)

        # Create xlsxwriter workbook object.
        workbook = writer.book
        # Create xlsxwriter worksheet object.
        worksheet = writer.sheets[sheet_name]
        # Create xlsxwriter format for column headers.
        head_format = workbook.add_format({'font_name': 'Arial',
                                           'font_size': 10,
                                           'bold': True})
        head_format.set_align("justify")
        head_format.set_align("top")
        head_format.set_text_wrap()
        # The excel column number the corresponding to the dataframe shape.
        end = xl_col_to_name(dataframe.shape[1])+"1"
        # Set excel columns format.
        worksheet.set_column('A1:'+end, 20)
        # Set headers format
        worksheet.set_row(start_row, 45, head_format)

        writer.save()

        print("File has been saved:", full_path)
        return

    @staticmethod
    def create_worksheet(dataframe,
                         file_name=None,
                         parent_file=None,
                         sheet_name=None,
                         start_row=3):
        """

        Parameters
        ----------

        Returns
        -------
        excel file
        """
        # Get the current working directory
        cwd = os.getcwd()

        # If no parent file is given then current working directory is used.
        parent_file = parent_file or cwd

        # If No file name is give then timestamp is used.
        file_name = file_name or StringTools.time_stamp()

        sheet_name = sheet_name or "data"

        # Make sure the string has the file extension added on
        if not file_name.endswith(".xlsx"):
            file_name = ".".join([file_name, 'xlsx'])

        full_path = "\\".join([parent_file, file_name])

        # Turn off header styles.
        pd.core.format.header_style = None

        # Create writer.
        writer = ExcelWriter(full_path,
                             engine='xlsxwriter',
                             options={'nan_inf_to_errors': True})

        # Use to_excel method to use writer on dataframe.
        dataframe.to_excel(writer,
                           sheet_name=sheet_name,
                           startrow=start_row,
                           index=False)

        # Create xlsxwriter workbook object.
        workbook = writer.book
        # Create xlsxwriter worksheet object.
        worksheet = writer.sheets[sheet_name]
        # Create xlsxwriter format for column headers.
        head_format = workbook.add_format({'font_name': 'Arial',
                                           'font_size': 10,
                                           'bold': True})
        head_format.set_align("justify")
        head_format.set_align("top")
        head_format.set_text_wrap()
        # The excel column number the corresponding to the dataframe shape.
        end = xl_col_to_name(dataframe.shape[1])+"1"
        # Set excel columns format.
        worksheet.set_column('A1:'+end, 20)
        # Set headers format
        worksheet.set_row(start_row, 45, head_format)

        return worksheet
