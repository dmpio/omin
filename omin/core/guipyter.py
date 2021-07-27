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

# from .cli_tools import CLITools
# from ..jtkinter import filedialog

# import os
import tkinter
from tkinter.filedialog import Open, Directory, SaveAs

# Find the cureent working directory.
here = os.getcwd()


def conditional_kwargs(**nkwargs):
    """Returns function with conditionally supplied kwargs.

    If a given keyword argument has not been supplied in the function
    definition then the keyword from the decorator will be substituted in.
    """
    def decorator(some_function):
        def wrapper(nkwargs=nkwargs, *args, **kwargs):
            for k, v in nkwargs.items():
                if k not in kwargs:
                    kwargs[k] = v
                else:
                    pass
            return some_function(*args, **kwargs)
        return wrapper
    return decorator


def root_topmost():
    """Return root as a withdrawn topmost window.
    """
    root = tkinter.Tk()
    # Hide the main window.
    root.withdraw()
    root.call('wm', 'attributes', ".", '-topmost', True)
    return root


class filedialog(object):
    """Jupyter Notebook friendly tkinter filedialogs."""
    # The following functions have been modified from the following;

    # https://github.com/python/cpython/blob/3.6/Lib/tkinter/filedialog.py

    @staticmethod
    @conditional_kwargs(parent=root_topmost(), initialdir=here)
    def askopenfilename(**options):
        "Ask for a filename to open."
        return Open(**options).show()

    @staticmethod
    @conditional_kwargs(parent=root_topmost(), initialdir=here)
    def asksaveasfilename(**options):
        "Ask for a filename to save as."
        return SaveAs(**options).show()

    @staticmethod
    @conditional_kwargs(parent=root_topmost(), initialdir=here)
    def askopenfilenames(**options):
        """Ask for multiple filenames to open.
        Returns a list of filenames or empty list if
        cancel button selected
        """
        options["multiple"] = 1
        return Open(**options).show()

    @staticmethod
    @conditional_kwargs(parent=root_topmost(), initialdir=here)
    def askopenfile(mode="r", **options):
        "Ask for a filename to open, and returned the opened file"

        filename = Open(**options).show()
        if filename:
            return open(filename, mode)
        return None

    @staticmethod
    @conditional_kwargs(parent=root_topmost(), initialdir=here)
    def askopenfiles(mode="r", **options):
        """Ask for multiple filenames and return the open file
        objects
        returns a list of open file objects or an empty list if
        cancel selected
        """

        files = askopenfilenames(**options)
        if files:
            ofiles = []
            for filename in files:
                ofiles.append(open(filename, mode))
            files = ofiles
        return files

    @staticmethod
    @conditional_kwargs(parent=root_topmost(), initialdir=here)
    def asksaveasfile(mode="w", **options):
        "Ask for a filename to save as, and returned the opened file"

        filename = SaveAs(**options).show()
        if filename:
            return open(filename, mode)
        return None

    @staticmethod
    @conditional_kwargs(parent=root_topmost(), initialdir=here)
    def askdirectory(**options):
        "Ask for a directory, and return the file name"
        return Directory(**options).show()



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
