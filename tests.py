"""Omin tests.

FIXME: ADD NOSETEST
"""
if __name__ == "__main__":
    import os
    import itertools
    import omin
    if "ExampleData" in list(os.walk("."))[0][1]:
        for n in list(os.walk("ExampleData")):
            # print(n)
            if len(n[-1]) > 1:
                # print(n[-1])
                files = list(itertools.product([n[0]], n[-1]))
                files = ["\\".join([j[0], j[-1]]) for j in files]
                print(files)
    # raw_data = omin.core.handles.RawData(peptide_groups_file, proteins_file)
