#ifndef __LIB_H_
#define __LIB_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

#define pow1(x) ((x) * (x))
#define rev(x) ((2 + (x)) % 4)
/*
#define max(x,y) ((x) > (y) ? (x) : (y))
#define min(x,y) ((x) < (y) ? (x) : (y))
*/
#define numc(i, j) ((i) < (j)? (j)*((j)+1)/2+(i): (i)*((i)+1)/2+(j))
#define LARGENUMBER 100000000
#define MAX_CON_LEN 5000000
#define MIN_CON_LEN 2000000
#define MIN_OVERLAP 5
#define MIN_LEG 1000
#define MBIG 1000000000
#define MSEED 16183398
#define MZ 0
#define FAC (1.0/MBIG)
#define MAXD 5
#define MAX_PRO_LEN 1000

#define total_nuc 16
#define na_name "actgnrywsmkhbvdx"
#define code_name "ACTG"
#define	total_aa 21
#define aalist "ACDEFGHIKLMNPQRSTVWY*"
#define codonlist "KNNKTTTTIIIMRSSRQHHQPPPPLLLLRRRR*YY*SSSSLFFL*CCWEDDEAAAAVVVVGGGG"
#define siteaa {8, 14} // Cleavage site: after Lys or Arg
#define nonsiteaa {12} // not cut before Pro

typedef struct edge {
	unsigned int	tag : 4 ; 	//4 bits for tag edges to be removed for various reasons; bit 1-3 indicates the edge is visited for offset 0-2, respectively
	bool	multitag;
	int	index;		/* index of the edge	*/
	char	*seq;
	int	length;
	int	bal_edge;
	int	nodeindex[2];
} EDGE;

typedef struct edgelist {
        EDGE    *edge;
        struct edgelist *next;
} EDGELIST;

typedef struct peptideseq {
	int	length;
	char	*seq;
	int	beg_edge_index;
	int	end_edge_index;
	int	beg_edge_offset;
	int	end_edge_offset;
	char	flag_stopcodon;
	struct peptideseq *parent;
	struct peptideseq *left;
	struct peptideseq *right;
} PEPTIDESEQ;

typedef struct peptidelink {
	int	index;
	char	*pepseq;
	int	length;
	char	*upper_seq, *down_seq;
	int	nextpep_index;
	int	prevpep_index;
} PEPTIDELINK;

typedef struct id_peptide {
	char	*seq;
	int	length;
	char	label;
	char	*prev_seq, *next_seq;
	int	len_prev_seq, len_next_seq;
	int	prev_id;
	int	next_id;
	int	beg_edge_index;
	int	end_edge_index;
	int	beg_edge_offset;
	int	end_edge_offset;
	char	flag_stopcodon;
} ID_PEPTIDE;

typedef struct edgemap {
	int	**id_pep_index; 	// sorted indices of identified peptides in frame 1, 2 and 3, respectively
	int	**id_pep_pos;
	int	*num_id_pep;
} EDGEMAP;

typedef struct u256bit {
	__int128_t	bits[2];
} uint256_t;

typedef struct vertex {
	unsigned int	tag : 1; // tagging visited vertices
	__int128_t	index;
	uint256_t	lindex;
	int	indegree, outdegree;
	EDGELIST	*lastedge;
	EDGELIST	*nextedge;
} VERTEX;

typedef struct hashpos {
	int readindex, pos;
	struct hashpos *next;
} HASH;

using namespace std;

#endif
