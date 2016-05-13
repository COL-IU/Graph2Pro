//Package: DBGraph2Pro (by H.T., Jul 2015)
//Program Name: DBGraph2Pro_main (the main program) 
//Latest modification on Jan 22th, 2014
//Yuzhen Ye <yye@indiana.edu> and Haixu Tang <hatang@indiana.edu>
//Affiliation: School of Informatics and Computing, Indiana University, Bloomington

#include "GraphTrav.h"
#include <time.h>

void printusage(char *error);

int main(int argc, char **argv)
{
	int copt, inp;
	int min_contig_leg = 500; //shortest contig length to be explored; default length 500
	int peptide_min = 8;
//	int peptide_min = 25;
	int peptide_max = 40;
//	int peptide_max = 15;
	int max_length = 3000;
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
	GraphTrav *eng = new GraphTrav();
	char edgefile[1000], edgeseqfile[1000];
	edgefile[0] = edgeseqfile[0] = 0; 
	char outfile[1000], singleoutfile[1000];
	strcpy(outfile, "tmp.out");

	// set defaults
	eng -> set_min_contig_leg(min_contig_leg);
	eng -> set_peptide_min(peptide_min);
	eng -> set_peptide_max(peptide_max);
	eng -> set_max_length(max_length);
	eng -> set_kmersize(kmer);
	eng -> set_codon2aanum(tmp_codon2aanum);
	eng -> set_cutsite(tryp_num, tryp_sites);
	eng -> set_mis_cleavage(mis_cleavage);
	eng -> set_max_depth(max_d);

	//get options
	while ((copt=getopt(argc,argv,"e:o:s:ul:p:m:k:c:d:f")) != EOF)	{
		switch(copt) {
			case 'e':
			  sscanf(optarg,"%s", edgefile);
			  continue;
			case 's':
			  sscanf(optarg,"%s", edgeseqfile);
			  continue;
			case 'o':
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 'l':
			  sscanf(optarg,"%d", &min_contig_leg);
			  eng -> set_min_contig_leg(min_contig_leg);
			  continue;
			case 'p':
			  sscanf(optarg,"%d", &peptide_min);
			  eng -> set_peptide_min(peptide_min);
			  continue;
			case 'm':
			  sscanf(optarg,"%d", &peptide_max);
			  eng -> set_peptide_max(peptide_max);
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
		printusage("Sequence file not specified");
	}

	//load assembly graph: inputs, edgefile & edgeseqfile
	if(FastG)	{	 // assembly graph in FastG format, only the contig file is needed.
		eng -> loadFastG(edgeseqfile);
	} else	{		//assembly graph in SOAP format, both the contig and graph files are needed.
		eng -> loadsoap(SOAP2, edgefile, edgeseqfile);
	}

        time_t t1 = time(NULL);
	eng -> Graph2Pro();
        time_t t2 = time(NULL);

	//write transcript assemblies
	eng -> writefile(outfile);

	delete eng;

        time_t tn = time(NULL);
	if(int(difftime(t2, t1)/60) == 0) {
	        printf("Total time used: %d sec (%d sec for computing proteins)\n",int(difftime(tn, t0)), int(difftime(t2, t1)));
	}
	else {
	        printf("Total time used: %d min (%d min for computing proteins)\n",int(difftime(tn, t0)/60), int(difftime(t2, t1)/60));
	}
}

void printusage(char *error)
{
	
	cout<<"Error: "<<error<<endl;
	cout<<"DBGraph2Pro version 0.1"<<endl;
	cout<<"Usage: DBGraph2Pro -e edgefile -s edgeseqfile -o outfile -p min_peptide_len -m max_peptide_len -l min_contig_len -k kmersize -u -c #-mis-cleavage -L Max_Seq_len -d Max_Depth\n";
	cout<<"-e edgeFile: The input edge file name"<<endl;
	cout<<"-s edgeSeqFile: The input edge sequence (contig) file name"<<endl;
	cout<<"-o OutFile(base name only): Protein Sequences files"<<endl;
	cout<<"-p min_peptide_len: minimum peptide length to be output (default 6)"<<endl;
	cout<<"-m max_peptide_len: maximum peptide length to be output (default 50)"<<endl;
	cout<<"-l min_contig_len: minimum contig length to be explored (longer contigs will be executed by FGS)"<<endl;
	cout<<"-L Max_Seq_len: maximum sequence length (for memory allocation, default 3000)"<<endl;
	cout<<"-k kmersize: default 31"<<endl;
	cout<<"-c mis-cleavage: default 0"<<endl;
	cout<<"-d max_depth: default 10"<<endl;
	cout<<"-u (SOAP when set; default off for SOAP2)"<<endl;
	cout<<"-f (FastG when set; default off for SOAP2)"<<endl;
	exit(-1);
}
