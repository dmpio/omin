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

import re
# All handles from core need to be imported.
from ..core.handles import *
from .select_process import select_process
from tkinter import Tk, filedialog

def select():
    """Generate instance of tkinter.filedialog.

    Parameters
    ----------
    None

    Returns
    -------
    selected_files: tuple

    """

    # Create Tk root
    root = Tk()

    # Hide the main window
    root.withdraw()

    # Raise the root to the top of all windows.
    root.call('wm', 'attributes', '.', '-topmost', True)

    selected_files = filedialog.askopenfilename(multiple=True)
    print("Omin will attempt to process the following files:")
    # If selected_files is a list then print then to show the user.
    if type(selected_files) == tuple:
        for i in selected_files:
            print(i)

    return selected_files


def selectInputPD(process_type = None):
    """
    Parameters
    ----------
    process_type : str

    Returns
    -------
    output : obj
        Whatever type of object that the user defines.
    """

    if process_type == None:
        process_type = select_process()

    pd_result_files = select()

    rx = re.compile("[Pp]eptide")

    peptide_file = list(filter(rx.findall, pd_result_files))[0]

    rx = re.compile("[Pp]roteins")
    protein_file = list(filter(rx.findall, pd_result_files))[0]

    function = eval(process_type)
    output = function(peptide_file, protein_file)

    return output
