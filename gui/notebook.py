# -*- coding: utf-8 -*-
"""Tools ipywidget based hybrid gui."""

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

import traitlets
# import warnings
from IPython.display import display
from ipywidgets import widgets
from tkinter import Tk, filedialog
from ..core.handles import Process
# warnings.filterwarnings("ignore")


class SelectFilesButton(widgets.Button):
    """A file widget that leverages tkinter.filedialog."""

    def __init__(self, *args, **kwargs):
        """Initialize the SelectFilesButton class."""
        super(SelectFilesButton, self).__init__(*args, **kwargs)
        # Add the selected_files trait
        self.add_traits(files=traitlets.traitlets.List())
        # Create the button.
        self.description = "Select Files"
        self.icon = "square-o"
        self.style.button_color = "orange"
        # Set on click behavior.
        self.on_click(self.select_files)

    @staticmethod
    def select_files(b):
        """Generate instance of tkinter.filedialog.

        Parameters
        ----------
        b : obj:
            An instance of ipywidgets.widgets.Button
        """
        # Create Tk root
        root = Tk()
        # Hide the main window
        root.withdraw()
        # Raise the root to the top of all windows.
        root.call('wm', 'attributes', '.', '-topmost', True)
        # List of selected fileswill be set to b.value
        b.files = filedialog.askopenfilename(multiple=True)

        b.description = "Files Selected"
        b.icon = "check-square-o"
        b.style.button_color = "lightgreen"

class RunButton(widgets.Button):
    """Button that begins processing the selected files."""

    def __init__(self, *args, **kwargs):
        """Initialize the SelectFilesButton class."""
        super(RunButton, self).__init__(*args, **kwargs)
        self.add_traits(files=traitlets.traitlets.List())
        self.add_traits(data=traitlets.traitlets.Instance(object))
        self.description = "Run"
        self.button_style = "info"
        self.icon = "play"
        # Enter the number
        self.process = "Process"
        self.on_click(self.run_process)

    @staticmethod
    def run_process(b):
        process = eval(b.process)
        b.data = process(b.files)


class OminNotebook(object):
    """A hybrid GUI for Omin."""

    def __init__(self):
        """Initialize the dashboard."""
        self.select_files_button = SelectFilesButton()
        self.run_button = RunButton()

    # @staticmethod
    def on_value_change(self,change):
        if type(change["old"]) != list:
            self.run_button.files = self.select_files_button.files
            display(self.run_button)

    @property
    def files(self):
        """Getter method for the selected files."""
        return self.select_files_button.files

    @property
    def data(self):
        return self.run_button.data

    @property
    def explore(self):
        view = lambda Attribute:display(self.data.__dict__[Attribute])
        vw = widgets.interactive(view,Attribute=list(self.data.__dict__.keys()))
        display(vw)

    def __repr__(self):
        """Show the dashboard on call."""

        display(self.select_files_button)
        self.select_files_button.observe(self.on_value_change, names="files")
        return ""
