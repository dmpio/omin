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

from dominate import tags
from ipywidgets import widgets
from IPython.display import display
from functools import partial

from ..utils.string_tools import StringTools

from .super_selection_container import SuperAccordion

from .widget_utils import widget_append
from .widget_utils import focus_new_child
from .widget_utils import SelectFilesButton
from .widget_utils import RunButton
from .widget_utils import SelectFilesPanel


class OminNotebookController(object):
    """A hybrid GUI for Omin."""

    def __init__(self):
        """Initialize the dashboard."""
        self._banner = "Omin Notebook"
        self._title = ""
        self._time_stamp = StringTools.time_stamp()
        # self._sections = ['Select Files',
        #                   'Process',
        #                   'Quick Stats',
        #                   'Results',
        #                   'Export']

        # Create new_analysis panel.
        self.new_analysis = SelectFilesPanel(selector_description='New Analysis')
        # Create Quick Stats
        self.quick_stats = widgets.HTML("")
        # Create the Run button.
        self.run_button = RunButton()
        # Link the files trait between new Analysis and run
        widgets.jslink((self.new_analysis.select_files, 'files'), (self.run_button, 'files'))
        # Create the Accordian widget
        # self.accordion = widgets.Accordion()
        self.accordion = SuperAccordion()
        # [self.accordion.set_title(i, title) for i, title in enumerate(self._sections)]
        self.accordion['Select Files'] = self.new_analysis
        # # Load the partial function for focussing on the new child widget.
        # loaded_focus_new_child = partial(focus_new_child,
        #                                  target=self.accordion)
        # # On the addition of a new child widget the accordion will focus on it.
        # self.accordion.observe(loaded_focus_new_child,
        #                        names='children',
        #                        type='change')

        # Create the header panel
        self.header_panel = widgets.VBox()
        self.header_panel.children = [self.header[0], self.header[1]]
        # Compose the main widget box.
        self.main = widgets.VBox([self.header_panel, self.accordion])

        # Observers

        # # Quick Stats observer
        # add_qs = partial(widget_append,
        #                  new_widget=self.quick_stats,
        #                  target=self.accordion)
        #
        # self.new_analysis.select_files.observe(add_qs,
        #                                        names='files',
        #                                        type='change')

        # # Set up for the Run Button.
        # add_rb = partial(widget_append,
        #                  new_widget=self.run_button,
        #                  target=self.accordion)
        #
        # self.new_analysis.select_files.observe(add_rb,
        #                                        names='files',
        #                                        type='change')

        # Set up for the Run Button.
        add_rb = partial(widget_append,
                         new_widget=self.run_button,
                         target=self.accordion,
                         title='Run Button')

        self.new_analysis.select_files.observe(add_rb,
                                               names='files',
                                               type='change')


    @property
    def banner(self):
        return self._banner

    @property
    def title(self):
        return self._title

    @property
    def time_stamp(self):
        return self._time_stamp

    @property
    def header(self):
        banner = widgets.HTML(tags.h1(self.banner).render())
        time_stamp = widgets.HTML(tags.h4(self.time_stamp).render())
        return banner, time_stamp

    @property
    def files(self):
        """Getter method for the selected files."""
        return self.select_files_button.files

    @property
    def data(self):
        """Getter method for the data held in button."""
        return self.run_button.data

    def __repr__(self):
        """Show the dashboard on call."""
        display(self.main)
        return ""
