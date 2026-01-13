from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
c_s=set()
for l in open(sys.argv[1]):
    c_s.add(l.rstrip())
with open(sys.argv[2]) as handle:
    for seq in SeqIO.parse(handle, "fasta"):
        if(seq.id not in c_s):
            sys.stdout.write(seq.format("fasta"))
            sys.stdout.flush()
handle.close()
