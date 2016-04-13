#ifndef __SOAP_H_
#define __SOAP_H_

#include "lib.h"

namespace soap
{
	int countedge(char *edgefile);
	int readedge(bool SOAP2, char *edgefile, char *edgeseqfile, EDGE *edge, VERTEX *vertex, int num_vertex, int kmersize);
	int insert_edge(EDGE *edge, __int128_t code1, __int128_t code2, VERTEX *vertex, int num_vertex);
	int link_vertex(EDGE *edge, VERTEX *vertex, int num_vertex);
	int loc_vertex(__int128_t code, VERTEX *vertex, int num_vertex);
	EDGELIST* add_edge_to_vertex(EDGE *edge, EDGELIST *edgelist0);
	void* delete_edgelist(EDGELIST *edgelist);

	//from readvertex.c
	__int128_t code2kmer(char *str);
	int index2code(__int128_t code, char *str);
	int readvertex(bool SOAP2, char *vertexfile, VERTEX *vertex);
	int encodekmer(char *printstr, int length, char *str);
	int printcode(char *str, char *printstr);
	char encodechar(char c);

};

#endif
