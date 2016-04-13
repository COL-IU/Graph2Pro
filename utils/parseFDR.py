from sys import argv
filename = argv[1]
FDRLevel = argv[2]
filename_FDR = argv[1]+"."+ FDRLevel +".tsv"
fw = open(filename_FDR,'w')
forward = 0.0
reverse = 0.0
titleHash = dict()
with open(filename) as file:
	for line in file:
		if line.startswith('#'):
			fw.write(line)
			continue
		lines = line.split('\t')
		title = lines[1]
		if title in titleHash :
			continue;
		else:
			titleHash[title] = 1
		evalue = float(lines[13])
		protein = lines[10]
		peptide = lines[9]
		peptide_new = peptide.replace('I','L')
		#if evalue <= 2.2E-11 :
		if protein.startswith('XXX_') :
			reverse += 1
		else:
			forward += 1
		spectrumQvalue = reverse/forward
		fw.write(line)
		if spectrumQvalue > float(FDRLevel) :
			break
fw.close()
print forward,"\t",reverse,"\t",spectrumQvalue,"\t",evalue

