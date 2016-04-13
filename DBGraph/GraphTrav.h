#ifndef __GRAPH2TRAV_H_
#define __GRAPH2TRAV_H_

#include "seqgraph.h"
//#include "codon.h"

class GraphTrav:public seqgraph
{
	private:
		int	num_cutsite;
		char	cutsite[20];
//		int	*peptidelen; 
		int	peptide_min; //only output peptides of this length or longer
		int	peptide_max; // only output peptides of this length or shorter 
		int	max_length; //length of the longest peptide
		int	num_peptides;  //number of peptides
		int 	min_contig_leg; // minimum length of contigs to be explored by FGS
		int	max_miscleavage; //maximum number of allowed miscleavage 
		int	max_pep;
		int	depth;
		int	max_depth;
		PEPTIDESEQ	*pepseq;	// peptide sequences sorted in a BST
		char	codon2aanum[64];
	public:	
		void	set_min_contig_leg(int inp) {min_contig_leg = inp; }
		void	set_peptide_min(int inp) {peptide_min = inp; }
		void	set_peptide_max(int inp) {peptide_max = inp; }
		void	set_max_depth(int inp) {max_depth = inp; }
		void	set_max_length(int inp) {max_length = inp; }
		void	set_codon2aanum(char *inp) {int i; for(i = 0; i < 64; i ++) codon2aanum[i] = inp[i]; }
		void	set_cutsite(int num, char *inp_cut)	{int i; num_cutsite = num; for(i = 0; i < num; i ++) cutsite[i] = inp_cut[i]; }
		void	set_mis_cleavage(int inp)	{max_miscleavage = inp; }
		int	appendpep(char *lastpep, int len_lastpep, char *tmppep, int length);
		char 	check_pepseq(char *seq1, int length1, char *seq2, int length2);
		void	writefile(char *outfile);
		PEPTIDESEQ* insert_pepseq(PEPTIDESEQ *pepseq, char *seq, int length);
		int	outputseq(ofstream &fs, PEPTIDESEQ *pepseq, int numpep);
		void 	output_peptides(ofstream &fp, int index, char *seq, int length, int beg_edge_index, int beg_edge_offset, int end_edge_index, int end_edge_offset, char flag_stopcodon);
		int 	copypepseq(char *lastpep, int len_lastpep, char begtag, int beg_edge_index, int beg_edge_offset, int end_edge_index, int end_edge_offset, char flag_stopcodon);
//		int	copypepseq(char **pepseq, int *peptidelen, int num_peptide, char *lastpep, int len_lastpep);
		int	trans_codon2aa(char *tmpcodon, int len_codonseq, char *pepseq);
		void	traverse_graph(EDGE *inedge, int offset, char *lastpep, int len_lastpep, int init_edge_index, int init_edge_offset);
		void	Graph2Pro(void);
		int	Loc_Stopcodon(char *seq, int length);
		char	extract_peptides(char *seq, int length, char label, char *lastpep, int *len_lastpep, int *leftover, int *ncut, char begtag, int init_edge_index, int *init_edge_offset, int end_edge_index, int *end_edge_offset);
		void	continue_extract_peptides(EDGE *inedge, int offset, char *lastpep, int *len_lastpep, int *leftover,
				 int *ncut, char begtag, int beg_edge_index, int beg_edge_offset, int end_edge_index, int end_edge_offset);
		char	ckcutsite(char seq);
};


#endif
