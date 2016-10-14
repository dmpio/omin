from urllib.request import urlopen

def upFastaGet(upID):
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
        urlout = urlopen(link)
        #fasta = urlout.read().decode("utf-8")
        fasta = urlout.read().decode("utf-8").split("\n")
    except TypeError as err:
        print("Please enter UniProt ID as string.")
    return fasta

def fasta2Seq(fasta):
    """ Joins strings in list into one string. Neglecting the first which is generally the annotation.

    Parameters
    ----------
    fasta : str

    Returns
    -------

    """
    one_string = "".join(fasta[1:])
    return one_string

def seqNum(seq):
    """Takes a sequence and returns the number for each amino acid (I think).

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
