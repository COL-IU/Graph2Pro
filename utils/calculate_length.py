from sys import argv
input_file = argv[1]

count = 0.0
length = 0.0

from Bio import SeqIO
for seq_record in SeqIO.parse(input_file,'fasta') :
	count += 1
	length += len(seq_record)
	#print count,"\t",len(seq_record)
print count,"\t",length/count
