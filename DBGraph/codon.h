int  tot_codon = 61;
char codon[61][3] = {
"ATT", "ATC", "ATA",		/* ILE */
"CTT", "CTC", "CTA", "CTG", "TTA", "TTG",  /* LEU */
"GTT", "GTC", "GTA", "GTG",	/* VAL	*/
"TTT", "TTC",			/* PHE	*/
"ATG",				/* MET */
"TGT", "TGC",			/* CYS */
"GCT", "GCC", "GCA", "GCG",	/* ALA */
"GGT", "GGC", "GGA", "GGG",	/* GLY */
"CCT", "CCC", "CCA", "CCG",	/* PRO */
"ACT", "ACC", "ACA", "ACG",	/* THR */
"TCT", "TCC", "TCA", "TCG", "AGT", "AGC", 	/* SER */
"TAT", "TAC",			/* TYR */
"TGG", 				/* TRP */
"CAA", "CAG",			/* GLN */
"AAT", "AAC",			/* ASN */
"CAT", "CAC",			/* HIS */
"GAA", "GAG",			/* GLU */
"GAT", "GAC",			/* ASP */
"AAA", "AAG",			/* LYS */
"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"	/* ARG */
};

char stop_codon[3][3] = {"TAG", "TAA", "TGA"};

char codon2aa[61] = "IIILLLLLLVVVVFFMCCAAAAGGGGPPPPTTTTSSSSSSYYWQQNNHHEEDDKKRRRRRR";

int total_aa = 21;
char aalist[21] = "ACDEFGHIKLMNPQRSTVWY*"

char codonlist[64] = "KNNKTTTTIIIMRSSRQHHQPPPPLLLLRRRR*YY*SSSSLFFL*CCWEDDEAAAAVVVVGGGG";

/*
AAA AAC AAT AAG
ACA ACC ACT ACG
ATA ATC ATT ATG
AGA AGC AGT AGG
CAA CAC CAT CAG
CCA CCC CCT CCG
CTA CTC CTT CTG
CGA CGC CGT CGG
TAA TAC TAT TAG
TCA TCC TCT TCG
TTA TTC TTT TTG
TGA TGC TGT TGG
GAA GAC GAT GAG
GCA GCC GCT GCG
GTA GTC GTT GTG
GGA GGC GGT GGG
*/
