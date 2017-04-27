"""jupyter notebook gui."""

from ipywidgets import widgets
from IPython.display import display
from .select_files import select


class SelectFiles(object):

    @staticmethod
    def on_select_files_click(b):
        return select()

    def __init__(self):
        button = widgets.Button(description="Select Files")
        button.on_click(self.on_select_files_click)
        display(button)
