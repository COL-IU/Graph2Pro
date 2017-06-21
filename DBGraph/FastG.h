#ifndef __FASTG_H_
#define __FASTG_H_

#include "lib.h"

namespace FastG
{
	int countedge(char *edgefile);
	int readedge(char *edgeseqfile, EDGE *edge, VERTEX *vertex, int *edgeindex, int num_edge, int kmersize);
	int get_next_edge(char *tmpstr, int *start, int end, int *edgeindex);
	uint256_t str2code(char *seq, int kmersize);
};

#endif
