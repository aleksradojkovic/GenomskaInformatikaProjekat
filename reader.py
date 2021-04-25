from Bio import SeqIO
import gzip

def GetSequencesFromFile(file):
    with gzip.open(file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            yield str(record.seq).upper()