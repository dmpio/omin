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
import numpy as np
import pandas as pd
from datetime import datetime
from dominate import tags
# import warnings
from IPython.display import display
# from IPython.display import HTML
from ipywidgets import widgets
from tkinter import Tk, filedialog
from ..core.handles import Process
from ..stats.tools import Compare
from ..visualize import Volcano
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
        """Run the process defined in the __init__ def."""
        process = eval(b.process)
        b.data = process(b.files)


class OminNotebook(object):
    """A hybrid GUI for Omin."""

    def __init__(self):
        """Initialize the dashboard."""
        self._header = "Omin Notebook"
        self._title = ""
        self._time_stamp = "{:%I:%M %p %A %B %d %Y}".format(datetime.now())
        self.select_files_button = SelectFilesButton()
        self.select_files_button.observe(self.on_value_change, names="files")
        self.run_button = RunButton()
        self.start_panel = widgets.VBox()
        self.start_panel.children = [self.top[0],
                                     self.top[1],
                                     self.select_files_button]

    def on_value_change(self, change):
        """Display the Run Button upon value change."""
        if type(change["old"]) != list:
            self.run_button.files = self.select_files_button.files
            # display(self.run_button)
            updated_panel = self.start_panel.children + (self.run_button,)
            self.start_panel.children = updated_panel

    @property
    def header(self):
        return self._header

    @property
    def title(self):
        return self._title

    @property
    def time_stamp(self):
        return self._time_stamp

    @property
    def top(self):
        header = widgets.HTML(tags.h1(self.header).render())
        time_stamp = widgets.HTML(tags.h4(self.time_stamp).render())
        return header, time_stamp

    @property
    def files(self):
        """Getter method for the selected files."""
        return self.select_files_button.files

    @property
    def data(self):
        """Getter method for the data held in button."""
        return self.run_button.data

    @property
    def explore(self):
        """Interactively explore the 'data' object."""
        view = lambda Attribute: display(self.data.__dict__[Attribute])
        vw = widgets.interactive(view,
                                 Attribute=list(self.data.__dict__.keys()))
        display(vw)

    def __repr__(self):
        """Show the dashboard on call."""
        display(self.start_panel)
        return ""

class MakeComparison(object):
    def __init__(self, dataframe=None):
        self.data = dataframe
        self.numerator = widgets.SelectMultiple(description="Numerator",
                                               options=list(dataframe.columns),
                                               layout=widgets.Layout(width="90%"))

        self.denominator = widgets.SelectMultiple(description="Denominator",
                                                 options=list(dataframe.columns),
                                                 layout=widgets.Layout(width="90%"))

        self.compute = widgets.Button(description="Compute")
        self.compare = widgets.VBox([self.numerator,
                                     self.denominator])
        self.compute.on_click(self.on_compute)

    def on_compute(self, b):
        num = self.data[list(self.numerator.value)]
        dem = self.data[list(self.denominator.value)]
        pvl = Compare.ttester(num.apply(np.log2), dem.apply(np.log2))
        lfc = Compare.log2FC(num.apply(np.log2), dem.apply(np.log2))
        # qvl = omin.Compare.bh_fdr(pvl)
        Volcano.simple(lfc, pvl)
        # plt.show()
        # TODO: Make this update the same graph seamlessly.

    def __repr__(self):
        display(self.compare)
        display(self.compute)

        return ""
