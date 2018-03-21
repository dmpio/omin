# -*- coding: utf-8 -*-
"""
omin.core.base
--------------

Provides, base class used in the Container classes.
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


# --------------
# REPR DECORATOR
# --------------

def repr_dec(cls):
    """Replace the __repr__ function of a given class with the repr_wrapper.

    When a class is decorated with this function each instance of that class
    will print a list of it's attributes along with that attributes type.
    Hopefully this will make navigating through class instances easier.

    Parameters
    ----------
    cls : class

    Returns
    -------
    cls : class
        With modified __repr__ function.
    """
    def repr_wrapper(self):
        """Show all attributes.
        """
        # Create a list of attributes.
        keys = sorted(list(self.__dict__.keys()))
        # Create a list of the attribute values types.
        value_types = list(map(lambda x: type(x).__name__, [self.__dict__[i] for i in keys]))
        # Zip the two lists above into a dict.
        att_dict = dict(zip(keys, value_types))
        # Create the template string for all the reaults.
        template = "{}: {}\n"
        # Create the results string.
        result = type(self).__name__+" Attributes\n"
        # Create a fancy border for the results.
        result += (len(result)-1)*"-" + "\n"
        # For every key and value in att_dict...
        for k,v in att_dict.items():
            # Format them and concatenate to result.
            result += template.format(k, v)
        return result

    # Replace the classes __repr__ function.
    cls.__repr__ = repr_wrapper
    return cls


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


# ------------
# HANDLE CLASS
# ------------


@repr_dec
class Handle(object):
    """The core omin handle base class.

    Attributes
    ----------
    numbers : dict
        DEPRECATE

    metadata : data

    type : type
    """

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
