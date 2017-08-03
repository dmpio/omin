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


import numpy as np
from dominate import tags
from ipywidgets import widgets
from IPython.display import display


from ..stats.tools import Compare
from ..visualize import Volcano

from ..utils.string_tools import StringTools
from .widget_utils import SelectFilesButton
from .widget_utils import RunButton


class OminNotebook(object):
    """A hybrid GUI for Omin."""

    def __init__(self):
        """Initialize the dashboard."""
        self._header = "Omin Notebook"
        self._title = ""
        # self._time_stamp = "{:%I:%M %p %A %B %d %Y}".format(datetime.now())
        self._time_stamp = StringTools.time_stamp()
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
