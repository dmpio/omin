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

    selected_files = filedialog.askopenfilename(multiple=True)
    return selected_files
