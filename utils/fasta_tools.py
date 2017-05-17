"""Tools for investigating with fasta and uniprot."""

class FastaTools(object):
    """A collection of tools for handling FASTA formatted strings."""

    @classmethod
    def fasta2Seq(cls, fasta):
        """Join strings in list into one string.

        Parameters
        ----------
        fasta : str

        Returns
        -------
        one_string : str
        """
        one_string = "".join(fasta[1:])
        return one_string

    @staticmethod
    def seqNum(cls, seq):
        """Take a sequence and returns the number for each amino acid.

        Notes
        -----
            Review this functions useage and rewrite annotation.

        Parameters
        ----------
        seq : str

        Returns
        -------
        seq_num : str

        """
        seq_num = [seq[i]+str(i+1) for i in range(len(seq))]
        return seq_num

# === UNIPROT TOOLS ===


class UniprotTools(object):

    @classmethod
    def getFasta(cls, upID):
        """Takes a uniprot ID and returns a string.

        Parameters
        ----------
        upID : str

        Returns
        -------
        fasta : str

        """
        try:
            fasta = str
            link = "http://www.uniprot.org/uniprot/"+upID+".fasta"
            try:
                urlout = urlopen(link)
                # print("hey")

            except HTTPError:
                print(
                      "Try again, the formatting may be off.",
                      "Try removing dashes."
                      )

            # fasta = urlout.read().decode("utf-8")
            fasta = urlout.read().decode("utf-8").split("\n")
            return fasta
        except TypeError as err:
            print("Please enter UniProt ID as string.")
