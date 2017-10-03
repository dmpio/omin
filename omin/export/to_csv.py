# -*- coding: utf-8 -*-
"""Tools for generating csv output."""

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


def inspect_object(obj):
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
    obj_ids = dict()
    # Try to make a list of the types of things inside obj.
    try:
        # Make list all things inside of an object
        for name, thing in obj.__dict__.items():
            obj_ids[name] = type(thing).__name__

        return obj_ids

    except Exception:
        pass

def object_walker(obj, desired_type=None, parent_dir=None):

    desired_type = desired_type or "DataFrame"

    if parent_dir == None:
        parent_dir = omin.gui.select_dir()

    no_list = dir(__builtins__)
    no_list.append("DataFrame")

    inspected = inspect_object(obj)

    for i in inspected.items():
        if i[-1] == desired_type:
            obj.__dict__[i[0]].to_csv(parent_dir+"/{}.csv".format(i[0]), index=False)

        if i[-1] not in no_list:
            dirn = "/".join([parent_dir,i[0]])
            # print(dirn)
            omin.IOTools.mkdir(dirn)
            # object_walker(obj.__dict__[i[0]], desired_type, parent_dir)
            object_walker(obj.__dict__[i[0]], desired_type, dirn)
