# -*- coding: utf-8 -*-
"""Visualization tools leveraging bokeh."""

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

# Load Bokeh tools
from bokeh.plotting import figure

from bokeh.models.widgets import (DataTable,
                                  TableColumn,
                                  HTMLTemplateFormatter)

from bokeh.models import (ColumnDataSource,
                          HoverTool,
                          TapTool,
                          OpenURL,
                          ToolbarBox,
                          )

from bokeh.io import output_notebook, push_notebook, show


# class DataSource(ColumnDataSource):
#     """Subclass of ColumnDataSource."""
#
#     def __init__(self, *args, **kwargs):
#         super(DataSource, self).__init__(*args, **kwargs)


class Bokomin(object):
    """A class for integration of Bokeh with Omin."""

    def __init__(self):
        """Initialize the base class."""
        output_notebook(hide_banner=True)
        pass

    def __repr__(self):
        """Show all attributes."""
        return "Attributes: "+", ".join(list(self.__dict__.keys()))