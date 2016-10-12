# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import string
from omin.norm import *

# def excellor(koabun, wtabun):
#     abun = pd.concat([koabun,wtabun],axis=1)
#     logabun = pd.concat([logger(koabun).com,logger(wtabun).com],axis=1)
#     logoave = pd.concat([logdiv(logger(koabun)),logdiv(logger(wtabun))],axis=1)
#     logopool = pd.concat([overPooler(koabun),overPooler(wtabun)],axis=1)
#     return pd.concat([abun,logabun,logabun,logoave,logopool],axis=1)
#
def excellor(compob, venn_list, parent_file, file_name):
    """

    Parameters
    ----------
    compob
    venn_list
    parent_file
    file_name

    Returns
    -------

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


def xLc(let):
    """
    xLc() short for Excel Columnizer
    
    Converts column letter labels from excel into a number that can be used as dataframe index.
    e.g.;
    xLc("A") returns 0, xLc("BA") returns 51
    """
    all_let = string.ascii_letters
    alph = all_let[all_let.index('A'):]
    if len(let) == 1:
        out = alph.index(let)
    else:
        out = (alph.index(let[0])+1)*26 + alph.index(let[1])
    return out

def cLx(i):
    alpha = string.ascii_letters
    alpha = alpha[alpha.index("A"):]
    if i > 25:
        isel = int(i/25)-1
        return alpha[isel]+alpha[(i-1)%25]
    else:
        return alpha[i]
    
#cLx(25)+":"+cLx(89)

def rgb2hex(rgb_tuple):
    """
    Takes a tuple of RGB values e.g. (253,253,245) and returns a 
    hexidecimal value for that color.
    """
    hexcolor = '#%02x%02x%02x' % rgb_tuple
    g = np.array(rgb_tuple)>256
    if g.any() == True:
        print("Please enter only values between 0 and 256")
    else:
        hexcolor = '#%02x%02x%02x' % rgb_tuple
        return hexcolor.upper()
