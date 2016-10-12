from urllib.request import urlopen

def upFastaGet(upID):
    """
    Takes a string as input and returns a fasta sequence as a list.
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
    """
    Joins strings in list into one string. Neglecting the first which is generally the annotation.
    """
    return "".join(fasta[1:])

def seqNum(seq):
    return [seq[i]+str(i+1) for i in range(len(seq))]