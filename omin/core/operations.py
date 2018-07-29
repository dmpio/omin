# -*- coding: utf-8 -*-
"""
omin.core.operations
--------------------

Provides operations and possible a main class entry point for the CLI.

"""
# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

# ----------------
# EXTERNAL IMPORTS
# ----------------
import os
import pickle
from glob import glob

# ----------------
# INTERNAL IMPORTS
# ----------------
from .handles import Process

# -------------
# UTILS IMPORTS
# -------------
from ..utils import IOTools

class Operate(object):

    @staticmethod
    def load_process(process_path):
        """Load an omin process file.
        """
        proc = None
        if os.path.exists(process_path):
            proc = pickle.load(open(process_path, "rb"))
            print("Loaded", process_path)
        else:
            print(process_path, "does not exist.")
        return proc


    @staticmethod
    def _dump_into_pickle_file(process_obj, process_path):
        """Dump a process obj into a pickle file.

        WARNING: Will overwrite existing files.
        """
        pickle.dump(process_obj, open(process_path, "wb"))
        print("Process has been saved at:", process_path)


    @classmethod
    def create_process(cls, process_obj, process_path, usr_inp=None, **kwargs):
        """Dump a process obj into a pickle file.

        Gives options on overwrite.
        """
        if os.path.exists(process_path):
            if usr_inp is None:
                usr_inp = input("A process file with the same name exists. Would you like to overwrite it? Y/N")

            if usr_inp == "Y" or usr_inp == "y":
                print("Process file", process_path, "is being overwritten...")
                cls._dump_into_pickle_file(process_obj, process_path)

            else:
                print("Process file", process_path, "was not overwritten.")
        else:
            cls._dump_into_pickle_file(process_obj, process_path)


    @classmethod
    def create_or_load_process(cls, process_path, overwrite=None, **kwargs):

        # If the process file already exists...
        if os.path.exists(process_path):
            # Collect user input...
            if overwrite is None:
                message = """A process with the same name exits in this directory.
                Would you like to:
                1)Load the file
                2)Overwrite the file
                """
                overwrite = input(message)
                overwrite = int(overwrite)
                overwrite = bool(overwrite-1)

            # If overwrite is True then overwrite the file.
            if overwrite:
                proc = Process(**kwargs)
                cls.create_process(proc, process_path, usr_inp="Y")

            # Exit without over writting.
            else:
                proc = cls.load_process(process_path)

            return proc
        # The path DNE the process file is witten.
        else:
            proc = Process(**kwargs)
            cls.create_process(proc, process_path)
            return proc

    @classmethod
    def start(cls, file_list=None, rescue_entrez_ids=True, *args, **kwargs):
        if file_list is None:
            # FIXME: Open this very narrow edge case up.
            # FIXME: Add error messages.
            file_list = glob("Proteome_Discoverer_Results/*.bz2")

        cls.create_or_load_process(file_list=file_list,
                                   rescue_entrez_ids=rescue_entrez_ids,
                                   *args,
                                   **kwargs)
