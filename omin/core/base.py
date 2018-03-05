# -*- coding: utf-8 -*-
"""omin.core.base

Provides
--------
Base class used in the Container classes.
"""
# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

# ----------
# TO DO LIST
# ----------
# FIXME: DOCUMENT OR DIE #DOD

# ----------------
# EXTERNAL IMPORTS
# ----------------
import re
import os
import guipyter as gptr
# FIXME: Should this be a try and except for pandas
from pandomics import pandas

# ----------------
# INTERNAL IMPORTS
# ----------------
from ..utils import IOTools


def export(obj, desired_type=None, parent_dir=None):
    """Export all attributes of an object that are DataFrames as csv files.
    """
    desired_type = desired_type or pandas.core.frame.DataFrame

    if parent_dir == None:
        parent_dir = gptr.filedialog.askdirectory()
        parent_dir = os.path.abspath(parent_dir)

    for i in obj._introspect().items():
        if i[-1] == desired_type:
            path = os.path.join(parent_dir, "{}.csv".format(i[0]))
            obj.__dict__[i[0]].to_csv(path)

        if issubclass(i[-1], Handle):
            dirn = os.path.join(parent_dir, i[0])
            # dirn = "/".join([parent_dir,i[0]])
            IOTools.mkdir(dirn)
            export(obj.__dict__[i[0]], desired_type, dirn)

# Omin's core handle
# ---------------------------
# Essentially container for DataFrames.

class Handle(object):
    """The core omin handle base class."""

    def __init__(self):
        """Initalize the core handle."""
        self.numbers = dict()
        self.metadata = dict()
        self.type = type(self)

    def _introspect(self):
        """Return list of object attributes and their types.
        """
        obj_ids = dict()
        # Try to make a list of the types of things inside obj.
        try:
            # Make list all things inside of an object
            for name, thing in self.__dict__.items():
                obj_ids[name] = type(thing)
            return obj_ids
        except Exception:
            pass

    def export_all(self, **kwargs):
        export(self, **kwargs)

    def __repr__(self):
        """Show all attributes."""
        return "Attributes: "+", ".join(list(self.__dict__.keys()))
