
import os
from glob import glob

from ..core import containers
from ..core import handles

here = os.getcwd()

example_data_dir = os.path.join(os.path.split(here)[0], "example_data", "_data")


# print(glob(example_data_dir+"/*"))

mcd_def_malonylation_file_list = os.path.join(example_data_dir, "mcd_def_malonylation", "*")

mcd_def_malonylation_file_list = glob(mcd_def_malonylation_file_list)


def test_containers_peptide_groups():
    peptide_groups = containers.PeptideGroups(mcd_def_malonylation_file_list[0])
    assert type(peptide_groups) != None

def test_containers_proteins():
    proteins = containers.Proteins(mcd_def_malonylation_file_list[-1])
    assert type(proteins) != None


def test_handles_process():
    process = handles.Process(file_list=mcd_def_malonylation_file_list)
    assert type(process) != None
