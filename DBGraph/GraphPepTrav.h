#ifndef __GRAPH2TRAV_H_
#define __GRAPH2TRAV_H_

#include "seqgraph.h"
//#include "codon.h"

class GraphPepTrav:public seqgraph
{
	private:
		int	num_cutsite;
		char	cutsite[20];
		int	num_id_peptides;
		ID_PEPTIDE *id_peptides; // identified peptides

		EDGEMAP *pep_edgemap;	// one for each edge

		int	num_proteins;  //number of output proteins
		int	*len_proteins;
		char	**seq_proteins; // output protein sequences
		char	**mark_proteins; // label identified peptides
		int	max_length;

		int	max_miscleavage; //maximum number of allowed miscleavage 
		int	depth;
		int	max_depth;

		int	max_pep;

		char	codon2aanum[64];
	public:	
		void	set_max_pep(int inp) {max_pep = inp; }
		void	set_max_depth(int inp) {max_depth = inp; }
		void	set_max_length(int inp) {max_length = inp; }
		void	set_codon2aanum(char *inp) {int i; for(i = 0; i < 64; i ++) codon2aanum[i] = inp[i]; }
		void	set_cutsite(int num, char *inp_cut)	{int i; num_cutsite = num; for(i = 0; i < num; i ++) cutsite[i] = inp_cut[i]; }
		void	set_mis_cleavage(int inp)	{max_miscleavage = inp; }

		void	loadpepmap(char *pepseqfile);
		void	input_pepnum(char *pepseqfile);
		void	input_pepmap(char *pepseqfile);

		void	pep2edgemap(void );
		void	sort_edgemap(void );

		void	markpeptide(char *markseq, int i1, int i2);
		int	extract_next_peptides(EDGE *inedge, int offset, char *lastseq, int *len_lastseq);
		int	extract_next_peptides_same_edge(EDGE *inedge, int offset, char *lastseq, int len_lastseq);
		int	extract_prev_peptides_same_edge(EDGE *inedge, int offset, char *lastseq, int len_lastseq);
		int	appendpep(char *lastpep, int len_lastpep, char *tmppep, int length);
		void	writefile(char *outputfile);
		void	outputseq(int index, ofstream &fp, char *seq, char *mark_proteins, int length);
		void	outnucseq(int index, ofstream &fp, char *seq, int length);
		int	trans_codon2aa(char *tmpcodon, int len_codonseq, char *pepseq);

		void	ConnectMapppedPeptides(void );
		void	GraphPep2Pro(void);
		void	TraverseGraph(void);
		char	ckcutsite(char seq);
		int	parse_pep(char *seq, int length);
};


#endif
