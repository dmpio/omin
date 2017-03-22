# -*- coding: utf-8 -*-
from tkinter import Tk, filedialog

def select():
    """Generate instance of tkinter.filedialog.

    Parameters
    ----------
    None

    Returns
    -------
    selected_files: tuple

    """

    # Create Tk root
    root = Tk()

    # Hide the main window
    root.withdraw()

    # Raise the root to the top of all windows.
    root.call('wm', 'attributes', '.', '-topmost', True)

    save_as_file_name = filedialog.asksaveasfilename()
    return save_as_file_name

if __name__ == "__main__":
    fn = select()
    print(fn)
