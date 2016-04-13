from sys import argv
input_file = argv[1]

from Bio import SeqIO
for seq_record in SeqIO.parse(input_file,'fasta') :
	if len(seq_record) >= 500 :
		print '>'+seq_record.id
		print seq_record.seq

