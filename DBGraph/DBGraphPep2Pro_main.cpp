//Package: DBGraphPep2Pro (by H.T., Jan 2016)
//Program Name: DBGraphPep2Pro_main (the main program) 
//Latest modification on Jan 27th, 2016
//Yuzhen Ye <yye@indiana.edu> and Haixu Tang <hatang@indiana.edu>
//Affiliation: School of Informatics and Computing, Indiana University, Bloomington

#include "GraphPepTrav.h"
#include <time.h>

void printusage(char *error);

int main(int argc, char **argv)
{
	int copt, inp;
	int min_contig_leg = 500; //shortest contig length to be explored; default length 500
	int max_length = 5000;
	int max_pep=100;
	int kmer = 31;
	int max_d = 20;
	bool SOAP2 = true;
	bool FastG = false;
	bool fastq = false;
	char tmp_codon2aanum[64] = {8, 11, 11, 8, 16, 16, 16, 16, 7, 7, 7, 10, 14, 15, 15, 14, 13, 6, 6, 13, 12,
		 12, 12, 12, 9, 9, 9, 9, 14, 14, 14, 14, 20, 19, 19,
		 20, 15, 15, 15, 15, 9, 4, 4, 9, 20, 1, 1, 18, 3, 2, 2, 3, 0, 0, 0, 0, 17, 17, 17, 17, 5, 5, 5, 5};
	int  mis_cleavage = 0;
	int  tryp_num = 2;
//	int  tryp_num = 0;	// only consider stop codon
	char tryp_sites[2] = {8, 14};

	time_t t0 = time(NULL);
	GraphPepTrav *eng = new GraphPepTrav();
	char edgefile[1000], edgeseqfile[1000], pepseqfile[1000];
	edgefile[0] = edgeseqfile[0] = pepseqfile[0] = 0; 
	char outfile[1000], outnucfile[1000];
	sprintf(outfile, "%s.out", pepseqfile);
	outnucfile[0] = 0;

	// set defaults
	eng -> set_max_length(max_length);
	eng -> set_kmersize(kmer);
	eng -> set_codon2aanum(tmp_codon2aanum);
	eng -> set_cutsite(tryp_num, tryp_sites);
	eng -> set_mis_cleavage(mis_cleavage);
	eng -> set_max_depth(max_d);
	eng -> set_max_pep(max_pep);

	//get options
	while ((copt=getopt(argc,argv,"e:o:n:p:s:uL:m:k:c:d:f")) != EOF)	{
		switch(copt) {
			case 'e':
			  sscanf(optarg,"%s", edgefile);
			  continue;
			case 's':
			  sscanf(optarg,"%s", edgeseqfile);
			  continue;
			case 'p':
			  sscanf(optarg,"%s", pepseqfile);
			  continue;
			case 'o':
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 'n':
			  sscanf(optarg,"%s", outnucfile);
			  continue;
			case 'm':
			  sscanf(optarg,"%d", &max_pep);
			  eng -> set_max_pep(max_pep);
			  continue;
			case 'L':
			  sscanf(optarg,"%d", &max_length);
			  eng -> set_max_length(max_length);
			  continue;
			case 'd':
			  sscanf(optarg,"%d", &max_d);
			  eng -> set_max_depth(max_d);
			  continue;
			case 'k':
			  sscanf(optarg,"%d", &kmer);
			  eng -> set_kmersize(kmer);
			  continue;
			case 'c':
			  sscanf(optarg,"%d", &mis_cleavage);
			  eng -> set_mis_cleavage(mis_cleavage);
			case 'u':
			  SOAP2 = false;
			  continue;
			case 'f':
			  FastG = true;
			  continue;
			default:
			  printusage("unknown input");
			}
	
		optind--;
	}

	if(edgefile[0] == 0 && !FastG) {
		printusage("Graph file not specified");
	}
	if(edgeseqfile[0] == 0) {
		printusage("Contig sequence file not specified");
	}
	if(pepseqfile[0] == 0) {
		printusage("Peptide sequence file not specified");
	}
	if(outnucfile[0] == 0)	{
		sprintf(outnucfile, "%s.out", outfile);
	}

	//load assembly graph: inputs, edgefile & edgeseqfile
	if(FastG)	{
		eng -> loadFastG(edgeseqfile);
	} else	{
		eng -> loadsoap(SOAP2, edgefile, edgeseqfile);
	}

printf("input peptides\n");
	//load identified peptides onto edges
	eng -> loadpepmap(pepseqfile);

printf("Connect peptides...\n");
        time_t t1 = time(NULL);
	eng -> ConnectMapppedPeptides();
        time_t t2 = time(NULL);
printf("Traverse the graph...\n");
	eng -> GraphPep2Pro();
        time_t t3 = time(NULL);

printf("Output proteins...\n");
	//write transcript assemblies
	eng -> writefile(outfile);

	delete eng;

        time_t tn = time(NULL);
	if(int(difftime(t3, t2)/60) == 0) {
	        printf("Total time used: %d sec (%d sec for computing proteins)\n",int(difftime(tn, t0)), int(difftime(t3, t2)));
	}
	else {
	        printf("Total time used: %d min (%d min for computing proteins)\n",int(difftime(tn, t0)/60), int(difftime(t3, t2)/60));
	}
}

void printusage(char *error)
{
	
	cout<<"Error: "<<error<<endl;
	cout<<"DBGraphPep2Pro version 0.1"<<endl;
	cout<<"Usage: DBGraphPep2Pro -e edgefile -s edgeseqfile -p PeptideSeqFile -n TranscriptOutFile -o ProteinOutFile -k kmersize -u -c #-mis-cleavage -L Max_Seq_len -d Max_Depth -m Max_Pep_per_Edge\n";
	cout<<"-e edgeFile: The input edge file name"<<endl;
	cout<<"-s edgeSeqFile: The input edge sequence (contig) file name"<<endl;
	cout<<"-p PepSeqFile: The input sequence (identified peptides) file name"<<endl;
	cout<<"-o ProteinOutFile(base name only): The output protein Sequences file name"<<endl;
	cout<<"-n TranscriptSeqFile: The output transcript sequences file name"<<endl;
	cout<<"-L Max_Seq_len: maximum protein sequence length (for memory allocation, default 3000)"<<endl;
	cout<<"-k kmersize: default 31"<<endl;
	cout<<"-c mis-cleavage: default 0"<<endl;
	cout<<"-d max_depth: default 10"<<endl;
	cout<<"-m Max_Pep_per_Edge: default 100"<<endl;
	cout<<"-u (SOAP when set; default off for SOAP2)"<<endl;
	cout<<"-f (FastG when set; default off for SOAP2)"<<endl;
	exit(-1);
}
