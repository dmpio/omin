# -*- coding: utf-8 -*-
"""
guipyter.utils.DataLoader
=========================

Provides tools to quickly funnel data in and out of Jupyter Notebooks.
"""

import _io
import os
import xlrd
import inspect

# Try to import pandas from pandomics
try:
    from panomics import pandas as pd
except ImportError as err:
    import pandas as pd

from .cli_tools import CLITools
from ..jtkinter import filedialog

class DataLoader(object):
    """Multi-purpose dataloading tools.
    """

    def __init__(self, filepath_or_buffer=None, *args, **kwargs):
        self.file_name = None
        self.raw = None
        self.filepath_or_buffer = filepath_or_buffer
        self.file_name = ''
        self.file_path = ''
        self.file_ext = ''

        # Sort out the parameters.
        read_table_params = inspect.signature(pd.read_table)
        read_excel_params = inspect.signature(pd.read_excel)
        read_excel_params = set(read_excel_params.parameters.keys())
        read_table_params = set(['filepath_or_buffer']) ^ set(read_table_params.parameters.keys())
        filedialog_params = set(["defaultextension", "filetypes", "initialdir", "initialfile", "multiple", "parent", "title", "typevariable"])

        # pd.read_table params collected here
        read_table_params_final = dict()
        for k,v in kwargs.items():
            if k in read_table_params:
                read_table_params_final[k]=v

        # pd.read_excel params collected here
        read_excel_params_final = dict()
        for k,v in kwargs.items():
            if k in read_excel_params:
                read_excel_params_final[k]=v

        # filedialog params collected here.
        filedialog_params_final = dict()
        for k,v in kwargs.items():
            if k in filedialog_params:
                filedialog_params_final[k]=v

        # Test if buffer has been passed.
        if isinstance(self.filepath_or_buffer, _io.TextIOWrapper):
            pass
            # Try to collect the file name.
            try:
                self.file_name = os.path.split(self.file_path)[-1]
            except Exception as err:
                print(err.args[0])


        # If self.filepath_or_buffer is filepath.
        elif isinstance(self.filepath_or_buffer, str):
            # FIXME: Add assert.
            self.file_path = self.filepath_or_buffer

        # If None Set the filepath_or_buffer from filedialog.
        elif self.filepath_or_buffer is None:
            self.filepath_or_buffer = filedialog.askopenfile(filedialog_params_final)

            # Try to collect the file path.
            try:
                self.file_path = self.filepath_or_buffer.name
            except AttributeError as err:
                print(err.args[0])

        # Try to collect the file name.
        try:
            self.file_name = os.path.split(self.file_path)[-1]
        except Exception as err:
            print(err.args[0])

        # Try to collect the file extension.
        try:
            self.file_ext = os.path.splitext(self.file_name)[-1]
        except Exception as err:
            print(err.args[0])

        if self.file_ext == ".xls" or self.file_ext == ".xlsx":
            if isinstance(self.filepath_or_buffer, _io.TextIOWrapper):
                self.filepath_or_buffer.close()
            else:
                pass

            try:
                # Excel files are so special.
                # The CLI will allow the user to select a page in an excel file.
                self.raw = CLITools.read_excel_cli(self.file_path,
                                                   **read_excel_params_final)

            except ValueError as e:
                print(e.args[0])

        else:
            try:
                # This should work with everything else.
                self.raw = pd.read_table(self.filepath_or_buffer,
                                         **read_table_params_final)

            except ValueError as e:
                print(e.args[0])

        if isinstance(self.filepath_or_buffer, _io.TextIOWrapper):
            self.filepath_or_buffer.close()

if __name__ == "__main__":
    dl = DataLoader()
