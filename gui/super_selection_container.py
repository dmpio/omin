# -*- coding: utf-8 -*-
"""Response ipywidgets conatiner."""

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


from traitlets import Tuple, Dict, CInt, Unicode
from ipywidgets import Box, CoreWidget
from ipython_genutils.py3compat import unicode_type


class SuperSelectionContainer(Box, CoreWidget):
    """Base class used to display multiple child widgets."""
    _titles = Dict(help="Titles of the pages").tag(sync=True)
    selected_index = CInt(help="""The index of the selected page.""", allow_none=True).tag(sync=True)

    def __init__(self, *args, **kwargs):
        # super(SuperSelectionContainer, self).__init__(*args, **kwargs)
        super().__init__(*args, **kwargs)
        # self.children = Tuple(allow_none=True).tag(sync=True)
        self.observe(self.focus_new_child, names='children', type='change')

        # Public methods
    def set_title(self, index, title):
        """Sets the title of a container page.

        Parameters
        ----------
        index : int
            Index of the container page
        title : unicode
            New title
        """
        # JSON dictionaries have string keys, so we convert index to a string
        index = unicode_type(int(index))
        self._titles[index] = title
        self.send_state('_titles')

    def get_title(self, index):
        """Gets the title of a container pages.

        Parameters
        ----------
        index : int
            Index of the container page
        """
        # JSON dictionaries have string keys, so we convert index to a string
        index = unicode_type(int(index))
        if index in self._titles:
            return self._titles[index]
        else:
            return None


    @property
    def reverse_titles(self):
        return dict((v, int(k)) for k, v in self._titles.items())

    def add_child(self, new):
        """Create a page for the added widget.
        """
        if hasattr(new, 'title'):
            if new.title not in self._titles.values():
                self.set_title(len(self.children), new.title)

        if new not in self.children:
            self.children += (new,)

    def remove(self, index):
        # If the number of children is less than t set the tuple to ().
        if len(self.children) is 1:
            self.children = ()

        if len(self.children) > 1:
            self.children = list(filter(lambda x: x is not self.children[index],
                                         self.children))

    def focus_new_child(self, change):
        if len(change['new']) is not len(change['old']):
            self.selected_index = len(self.children) - 1

    def __add__(self, other):
        # self.children += (other,)
        self.add_child(other)

    def __getitem__(self, value):
        try:
            return self.children[self.reverse_titles[value]]
        except:
            return self.children[value]

    def __setitem__(self, key, value):
        child_list = list(self.children)
        if type(key) is str:
            if key in self.reverse_titles:
                child_list[self.reverse_titles[key]] = value
            else:
                self.set_title(len(self.children), key)
                # child_list[self.reverse_titles[key]] = value
                child_list.append(value)

        if type(key) is int:
            child_list[key] = value

        self.children = (child_list)

class SuperAccordion(SuperSelectionContainer):
    """Displays children each on a separate accordion page."""
    _view_name = Unicode('AccordionView').tag(sync=True)
    _model_name = Unicode('AccordionModel').tag(sync=True)

class SuperTab(SuperSelectionContainer):
    """Displays children each on a separate tab page."""
    _view_name = Unicode('TabView').tag(sync=True)
    _model_name = Unicode('TabModel').tag(sync=True)
