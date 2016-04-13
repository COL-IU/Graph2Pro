#from Bio.Seq import MutableSeq
#from Bio.Alphabet import IUPAC
from sys import argv
filename = argv[1]
from Bio import SeqIO
fileHash = dict();
write_filename = argv[1].replace('.fasta','.fixedKR.fasta')
print write_filename
f = open(write_filename,'w')

for seq_record in SeqIO.parse(filename, "fasta"):
        refID = seq_record.id
	des = str(seq_record.description)
	#print des
        seq = str(seq_record.seq)
	seq_rev = ''
	if seq[len(seq)-1] == 'K' or seq[len(seq) -1] == 'R' :
		seq_2 = seq[0:len(seq)-1]
		seq_rev = seq_2[::-1]+seq[len(seq)-1]
	else:
		seq_rev = seq[::-1]
	f.write(">")
	f.write(des)
	f.write("\n")
	f.write(seq)
	f.write("\n")
	f.write(">XXX_")
	f.write(des)
	f.write("\n")
	f.write(seq_rev)
	f.write("\n")
f.close()

