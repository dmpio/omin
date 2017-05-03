"""jupyter notebook gui."""

from ipywidgets import widgets
from IPython.display import display
from tkinter import Tk, filedialog

class OminNotebook(object):

    def __init__(self):
        # Create the select_files_button.
        self.select_files_button = widgets.Button(description="Select Files", icon="file-text-o")
        self.select_files_button.style.button_color = "orange"
        # Create value attribute of select_files_button.
        setattr(self.select_files_button, "value", None)
        # Set on click behavior.
        self.select_files_button.on_click(self.select_files)

        # Select processing style.
        self.process_select = widgets.Select(
            options=['RawData', 'PreProcess', 'Process'],
            value='Process',
            description='Select Process:',
            disabled=False)

        # Run button
        self.run_button = widgets.Button(description="RUN!",icon="rocket")
        self.run_button.style.button_color = "lightgreen"

    @staticmethod
    def select_files(b):
        """Generate instance of tkinter.filedialog.
        """
        # Create Tk root
        root = Tk()
        # Hide the main window
        root.withdraw()
        # Raise the root to the top of all windows.
        root.call('wm', 'attributes', '.', '-topmost', True)
        files = filedialog.askopenfilename(multiple=True)
        b.value = files

    def __repr__(self):
        display(self.select_files_button,
                self.process_select,
                self.run_button)
        return ""
