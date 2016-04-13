#include "smallapp.h"


FILE* smallapp::ckopen(char *name, char *mode)
{
	FILE *fp;

	if ((fp = fopen(name, mode)) == NULL)	{
		printf("Cannot open %s.\n", name);
		exit(-1);
	}
	return(fp);
}


/* ckalloc - allocate space; check for success */

void* smallapp::ckalloc(int amount)
{
	void *p;

	if(amount == 0)	{
		amount = (unsigned) 100;
	}
	if ((p = (void *) calloc( (unsigned) amount, 1)) == NULL)	{
		printf("Ran out of memory.\n");
                printf("There may be errors as follows:\n");
                printf("1) Not enough memory.\n");
                printf("2) The ARRAY may be overrode.\n");
                printf("3) The wild pointers.\n");
                exit(-1);
	}
	return(p);
}

double smallapp::random1(int *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj, mk;
	int i, ii,k;

	if(*idum < 0 || iff == 0)	{ /* initialization */
		iff = 1;
		mj = MSEED - (*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;
		for(i = 1; i < 54; i ++)	{
			ii = (21 + i) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if(mk < MZ)	mk += MBIG;
			mj = ma[ii];
		}
		for(k = 1; k <= 4; k ++)	{
			for(i = 1; i <= 55; i ++)	{
				ma[i] -= ma[1 + (i + 30) % 55];
				if(ma[i] < MZ)	ma[i] += MBIG;
			}
		}
		inext = 0;
		inextp = 31;
		*idum = 1;
	}
	if( ++inext == 56) inext = 1;
	if( ++inextp == 56) inextp = 1;
	mj = ma[inext] - ma[inextp];
	if(mj < MZ) mj += MBIG;
	ma[inext] = mj;
	return(mj * FAC);
}

int  smallapp::ran_number(int n, int *idum)
{
	double	     t;
	int          p;

	t = random1(idum);
	if(t == 1.0)	{
		t = t - 1.0e-10;
	}
	p = ( int ) (n * t);
	return( p );
}

int smallapp::hashvalue(int hashw)
{
        int     i, k;
        k = 1;
        for(i = 0; i < hashw; i ++)     {
                k *= 4;
        }
        return(k);
}

//sort vertex based on their index (for fast retrival later)
int smallapp::comparvertex(const void * p1, const void * p2)
{
        VERTEX *v1 = (VERTEX *)p1;
        VERTEX *v2 = (VERTEX *)p2;
        if(v1 -> index > v2 -> index)   {
                return(1);
        } else if(v1 -> index < v2 -> index)    {
                return(-1);
        } else  {
                return(0);
        }
}

char smallapp::strand2char(char *strandtmp)
{
        if(strandtmp[0] == '+') {
                return(0);
        } else  {
                return(1);
        }
}
