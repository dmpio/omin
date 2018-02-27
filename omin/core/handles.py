# -*- coding: utf-8 -*-
"""
Omin core handles.

Handle in this context is a class composed of several pandas DataFrames, and
other varibles that are either derived from the DataFrames or provided by the
user.
"""
# -------
# LiCENSE
# -------
# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

# ----------
# TO DO LIST
# ----------
# FIXME: DOCUMENT OR DIE #DOD


# ----------------
# EXTERNAL IMPORTS
# ----------------
import re

# ----------------
# INTERNAL IMPORTS
# ----------------
from .containers import PeptideGroups, Proteins

# =============
# PROJECT CLASS
# =============

class Project(object):
    """
    """
    def __init__(self, file_list=None, peptides_file=None, proteins_file=None, verbose=None, *args, **kwargs):
        """Load data for peptides_file and proteins_file
        """

        verbose = verbose or True

        if file_list is not None:
            rx = re.compile("[Pp]eptide")
            peptides_file = list(filter(rx.findall, file_list))[0]

            rx = re.compile("[Pp]roteins")
            proteins_file = list(filter(rx.findall, file_list))[0]

        try:
            self.peptide_groups = PeptideGroups(peptides_file)
        except Exception as err:
            if verbose:
                print(err)
            self.peptide_groups = PeptideGroups(*args, **kwargs)

        try:
            self.proteins = Proteins(peptides_file)
        except Exception as err:
            if verbose:
                print(err)
            self.proteins = Proteins(*args, **kwargs)

# =============
# PROCESS CLASS
# =============

class Process(Project):
    """A metaclass that uses PreProcess attempting several normalization steps.

    WARNING: This class is under construction switch to stable branch it you need to work.

    """

    # PROTIP: Wait util the last possible moment to link databases.
    # microprotip: Exceptions are the rule here.
    # PROTIP: TRY NOT TO SET VARS AT THIS LEVEL.

    def __init__(self, file_list=None, peptides_file=None, proteins_file=None, modifications=None, *args, **kwargs):
        """Initalize Process class.
        """

        Project.__init__(self, *args, **kwargs)
