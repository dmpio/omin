# -*- coding: utf-8 -*-
# Copyright 2017 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach,
# Blair Chesnut, and Elizabeth Hauser.
import re
import os
import pandas as pd
import guipyter as gptr
from ..utils import IOTools


def export(obj, desired_type=None, parent_dir=None):
    """Export all attributes of an object that are DataFrames as csv files.
    """
    desired_type = desired_type or pd.DataFrame

    if parent_dir == None:
        parent_dir = gptr.filedialog.askdirectory()
        parent_dir = os.path.abspath(parent_dir)

    for i in obj._introspect().items():
        if i[-1] == desired_type:
            path = os.path.join(parent_dir, "{}.csv".format(i[0]))


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

    numbers = dict()
    metadata = dict()
    # type = type(self)

    # def __init__(self):
    #     """Initalize the core handle."""
    #     self.numbers = dict()
    #     self.metadata = dict()
    #     self.type = type(self)

    @classmethod
    def _introspect(cls):
        """Return list of object attributes and their types.
        """
        obj_ids = dict()
        # Try to make a list of the types of things inside obj.
        try:
            # Make list all things inside of an object
            for name, thing in cls.__dict__.items():
                obj_ids[name] = type(thing)

            return obj_ids

        except Exception:
            pass

    @classmethod
    def export_all(cls, **kwargs):
        export(cls, **kwargs)

    @classmethod
    def __repr__(cls):
        """Show all attributes."""
        return "Attributes: "+", ".join(list(cls.__dict__.keys()))
