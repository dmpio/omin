# -*- coding: utf-8 -*-
"""
Copyright 2017 James Draper, Paul Grimsrud, Deborah Muoio

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files, Omics Modeling Integrating
Normalization (OMIN), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom
the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM.
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

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
