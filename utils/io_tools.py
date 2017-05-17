

class IOTools(object):
    """Tools for file handling and dir building."""

    @staticmethod
    def mkdir(directory):
        """Create a directory.

        Parameters
        ----------
        directory : str
        """
        assert type(directory) == str
        # directory = StringTools.remove_punctuation(directory)
        directory = directory.replace(" ", "_")
        if not os.path.exists(directory):
            os.makedirs(directory)
            # print(directory, "is ready.")
            return directory
        else:
            # print(directory, "already exists.")
            return directory
# Testing

if __name__ == "__main__":
    print("Testing utils.py ...")
    from omin.core.handles import RawData

    try:
        data = RawData(
            "ExampleData\crat_ex\_E749_4154_010716_PeptideGroups.txt",
            "ExampleData\crat_ex\_E749_4154_010716_Proteins.txt"
            )
        if data.raw_peptides.shape == (8712, 277):
            print("ExampleData has loaded correctly.")
    except Exception:
        print("Loading ExampleData failed.")
    try:
        if type(inspectObject(data)) == list:
            print("inspectObject works!")
    except Exception:
        print("inspectObject has failed.")
