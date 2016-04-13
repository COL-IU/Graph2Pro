#ifndef __SMALLAPP_H_
#define __SMALLAPP_H_

#include "lib.h"

namespace smallapp
{
	void* 	ckalloc(int amount);
	FILE*   ckopen(char *name, char *mode);
	int  	ran_number(int n, int *idum);
	double 	random1(int *idum);
	int 	hashvalue(int hashw);
	int 	comparcontig(const void * p1, const void * p2);
	int 	comparint(const void * p1, const void * p2);
	int	comparvertex(const void * p1, const void * p2);
	char 	strand2char(char *strandtmp);
};

#endif
