#ifndef __READAPP_H_
#define __READAPP_H_

#include "lib.h"

namespace seq
{
	int readseq_num(char *filename);
	int readseq(char **src_seq, char **src_name, int *len_seq, char *filename);
	char aa2int(char c);
	char char2int(char c);
	char char2intgen(char c);
	int outputseq(ofstream &fp, char *seq, char *name, int pbeg, int pend);
	int writeseq(ofstream &fp, char *seq, char *name, int len_seq);
};

#endif
