# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import string
# from omin.norm import *


def excellor(compob, venn_list, parent_file, file_name):
    """Takes objects containing DataFrames and outputs those as Excel files.

    Parameters
    ----------
    compob : (:obj)
    venn_list : list
    parent_file : str
    file_name : str

    Returns
    -------
    excel file
    """
    from pandas import ExcelWriter
    import re

    # Make sure the string has the file extension added on
    if not file_name.endswith(".xlsx"):
        file_name = file_name + '.xlsx'

    fn = re.sub('[^0-9a-zA-Z]+',  # regex pattern
                '_',  # replacement
                file_name.split(".")[0])  # input string

    file_name = fn + "." + file_name.split(".")[1]

    writer = ExcelWriter(parent_file + '/' + file_name)
    for i in venn_list:
        compob.__dict__[i].to_excel(writer, i)
    writer.save()
    print("Your file has been saved:", parent_file, "/", file_name)
