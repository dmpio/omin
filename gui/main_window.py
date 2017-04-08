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

import tkinter
from tkinter import filedialog

class omin_app(tkinter.Tk):
    def __init__(self):
        self.selected_files = None
        self.omin_obj = None

    def run(self):
        root = tkinter.Tk()
        root.title(('Omin GUI Demo'))
        root.resizable(width=True, height=True)
        root.minsize(width=600, height=200)
        # root.maxsize(width=666, height=666)
        root.call('wm', 'attributes', '.', '-topmost', True)
        tkinter.Button(root, text=('Select Input'), command=self.select_files, width=10).pack(side=tkinter.TOP)
        tkinter.Listbox(root).pack(side=tkinter.BOTTOM)
        root.mainloop()


    def select_files(self):
        selected_files = filedialog.askopenfilename(multiple=True)
        self.selected_files = selected_files
        return


if __name__ == "__main__":
    oa = omin_app()
    oa.run()
