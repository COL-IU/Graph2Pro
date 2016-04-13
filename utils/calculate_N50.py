from sys import argv
fileR = argv[1]
from Bio import SeqIO

tmpLength = 0.0
lengthArray = []
totalLength = 0.0
for seq_record in SeqIO.parse(fileR,'fasta') :
        if len(seq_record) >= 500 :
		lengthArray.append(len(seq_record.seq))
		totalLength += len(seq_record.seq)

N50L = totalLength/2.0
lengthArray.sort(reverse=True)
N50 = 0
for entry in lengthArray :
	tmpLength += entry
	if tmpLength >= N50L:
		N50 = entry
		break
print N50

