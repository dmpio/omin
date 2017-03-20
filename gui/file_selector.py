import re
# All handles from core need to be imported.
from ..core.handles import *
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


def selectInputPD(process_type):
    """
    Parameters
    ----------
    process_type : str

    Returns
    -------
    output : obj
        Whatever type of object that the user defines.
    """

    pd_result_files = select()

    rx = re.compile("[Pp]eptide")

    peptide_file = list(filter(rx.findall, pd_result_files))[0]

    rx = re.compile("[Pp]roteins")
    protein_file = list(filter(rx.findall, pd_result_files))[0]

    function = eval(process_type)
    output = function(peptide_file, protein_file)

    return output
