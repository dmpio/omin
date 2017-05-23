# -*- coding: utf-8 -*-
"""Tools for handling objects."""
# LICENSE
# -------

# Copyright 2017 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach,
# Blair Chesnut, and Elizabeth Hauser.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files, (the software)), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions: The above copyright
# notice and this permission notice shall be included in all copies or
# substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS",
# WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
# TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM. OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import re


def inspectObject(obj):
    """Return list of object attributes and their types.

    Parameters
    ----------
    obj: :obj:
        Any kind of object that is not a builtin.

    Returns
    -------
    obj_ids: list
        A list of lists to be specific.
    """
    obj_ids = []
    # Make list of builtins
    # bi_list = dir(__builtins__)
    # Try to make a list of the types of things inside obj.
    try:
        # Make list all things inside of an object
        for name, thing in obj.__dict__.items():
            obj_ids.append([name, type(thing).__name__, thing])

        return obj_ids

    except Exception:
        pass


def objectWalker(obj, desired_type=None, att_list=None):
    """Recursively walks through an object an genrates a list of its attributes.

    Parameters
    ----------
    obj: :obj:

    desired_type: str

    Returns
    -------
    att_list: list

    """
    # Create an attributes list if none has been asigned.
    if att_list is None:
        att_list = []

    # Create list of types that we do not want to search through.
    no_list = dir(__builtins__)
    no_list.append("DataFrame")
    desired_type = desired_type or "DataFrame"

    obj_inspect = inspectObject(obj)

    for i in obj_inspect:

        # Test if object attribute is desired_type.
        if re.search(desired_type, i[1], re.I) is not None:

            att_list.append(i)

        # If the object attribute is not the desired_type then...
        elif i[1] not in no_list:
            # Try to use recursively objectWalker on the object.
            try:
                objectWalker(obj.__dict__[i[0]], desired_type, att_list)
            except Exception:
                pass
    return att_list