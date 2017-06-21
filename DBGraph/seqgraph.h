#ifndef __SEQGRAPH_H_
#define __SEQGRAPH_H_

#include "lib.h"

class seqgraph {
	protected: 
		int kmersize;
		__int128_t maskv;
		int num_vertex;
		int *vindex;
		VERTEX *vertex; //??

		int num_edge;
		EDGE *edge; //contigs -- contain seq etc

		int hashw;
		int hashtablesize;

	public:
		seqgraph(void);
		~seqgraph(void) { cleanup(); }
		void init();

		void cleanup(void);
		void set_kmersize(int inp) { kmersize = inp; set_maskv(); }
		void set_maskv(void);
		void loadsoap(bool SOAP2, char *edgefile, char *edgeseqfile);
		void loadFastG(char *edgeseqfile);
		void index_vertex(void);
		void lindex_vertex(void);
		void print_vertex_degree(void);
};

#endif
