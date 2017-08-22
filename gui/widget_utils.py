# -*- coding: utf-8 -*-
"""Customized ipywidgets."""

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
from IPython.display import display
from IPython.display import HTML
from ipywidgets import widgets
# from traitlets import traitlets
from functools import partial
from datetime import datetime
from tkinter import Tk, filedialog
from dominate import tags
from .super_selection_container import SuperAccordion
from ..core.handles import Process
from warnings import warn


def timestamp():
    """Return a ipywidget timestamp."""
    warn("deprecated: omin.gui.timestamp, use omin.StringTools.timestamp")
    ts = tags.h4("{:%I:%M:%S %p %A %B %d %Y}".format(datetime.now()))
    # ts_widget = widgets.HTML(ts.render())
    ts_widget = HTML(ts.render())
    return ts_widget


def widget_append(change, target, new_widget, title=None):
    """Append new widget to a taget layout."""
    if change['new'] is not change['old']:
        # Only for SuperSelectionContainer
        if isinstance(new_widget, SuperAccordion):
            if title is not None:
                if new_widget not in target.children:
                    print('success')
                    target[title] = new_widget
                    return
        else:
            if new_widget not in target.children:
                target.children += (new_widget,)
                return


def toggle_append(change, target, new_widget):
    """Append new widget to a taget layout on a toggle change."""
    if change['new'] is True:
        if new_widget not in target.children:
            target.children += (new_widget,)
            return
    # Remove the widget if toggle is set to false
    if change['new'] is not False:
        if new_widget in target.children:
            target.children = list(filter(lambda x: x is not new_widget,
                                          target.children))


def shift_focus(layout, index):
    """Change the selected_index of a layout widget on the fly."""
    if hasattr(layout, 'selected_index'):
        layout.selected_index = index


def focus_new_child(change, target=None):
    """On the addition of a child widget a target layout will focus on it."""
    if len(change['new']) is not len(change['old']):
        if target is not None:
            assert hasattr(target, 'selected_index')
            target.selected_index = len(target.children) - 1


class LoadedButton(widgets.Button):
    """A button that can holds a value as a attribute."""

    def __init__(self, value=None, *args, **kwargs):
        """Initalize the loaded button class."""
        super(LoadedButton, self).__init__(*args, **kwargs)
        # Create the value attribute.
        self.add_traits(value=traitlets.Any(value))


class SelectFilesButton(widgets.Button):
    """A file widget that leverages tkinter.filedialog."""
    files = traitlets.List([], help="List of file paths").tag(sync=True)

    def __init__(self, *args, **kwargs):
        """Initialize the SelectFilesButton class."""
        super(SelectFilesButton, self).__init__(*args, **kwargs)
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

    files = traitlets.List([], help="List of file paths").tag(sync=True)

    def __init__(self, *args, **kwargs):
        """Initialize the SelectFilesButton class."""
        super(RunButton, self).__init__(*args, **kwargs)
        self.add_traits(data=traitlets.Instance(object))
        self.description = "Run"
        self.button_style = "info"
        self.icon = "play"
        self.process = "Process"
        self.on_click(self.run_process)

    @staticmethod
    def run_process(b):
        """Run the process defined in the __init__ def."""
        process = eval(b.process)
        b.data = process(b.files)


class SelectFilesPanel(widgets.VBox):
    """General File Selection Panel."""

    def __init__(self, selector_description="", *args, **kwargs):
        """Initialize the SelectFilesPanel."""
        super(SelectFilesPanel, self).__init__(*args, **kwargs)
        # Create selector button.
        self.files = []
        self.selector = widgets.ToggleButton(description=selector_description,
                                             value=False)
        self.select_files = SelectFilesButton()
        self.children = [self.selector]
        # self.panel = widgets.VBox([self.selector])
        # Load the kwargs into trigger.
        loaded_toggle_append = partial(toggle_append,
                                       target=self,
                                       new_widget=self.select_files)
        # If selector is then selecte_files is added to panel.
        self.selector.observe(loaded_toggle_append,
                              names="value", type='change')

    # def __repr__(self):
    #     """Show the panel."""
    #     display(self.panel)
    #     return ""
