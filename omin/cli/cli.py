# -*- coding: utf-8 -*-

# LICENSE
# -------

# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach,
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



import argparse
import sys
import time

import guipyter
import guipyter as gptr
# import omin
from ..core import handles


class FileHandleCLI(object):

    def __init__(self):
        self.proteins_file_name = None
        self.peptide_groups_file_name = None

    def set_proteins_file_name(self):
        print("Select proteins file:\n")
        time.sleep(1)
        # self.proteins = gptr.DataLoader()
        self.proteins_file_name = gptr.filedialog.askopenfilename()

    def set_peptide_groups_file_name(self):
        time.sleep(1)
        print("Select peptide groups file:\n")
        # self.peptide_groups= gptr.DataLoader()
        self.peptide_groups_file_name = gptr.filedialog.askopenfilename()

    def get_proteins_file_name(self):
        return self.proteins_file_name

    def get_peptide_groups_file_name(self):
        return self.peptide_groups_file_name

    def show_selected(self):
        self.get_proteins_file_name()
        self.get_peptide_groups_file_name()



def main():
    title = "omin normalize now:\n"
    print(len(title)*"-")
    print(title)
    print(len(title)*"-",'\n')

    fh = FileHandleCLI()
    input("Press Enter to continue...\n")
    fh.set_peptide_groups_file_name()
    fh.set_proteins_file_name()
    fh.show_selected()

    om = handles.Process(peptides_file=fh.get_peptide_groups_file_name(),
                         proteins_file=fh.get_proteins_file_name())

    print(om.raw_proteins.shape)
