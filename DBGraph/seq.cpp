#include "seq.h"
#include "smallapp.h"

int seq::readseq_num(char *filename)
{
	std::ifstream fp(filename);
	if(!fp) { cout<<"open file error "<<filename<<endl; exit(0); }
	int n = 0;
	std::string s;
	while(!fp.eof()) {
		std::getline(fp, s);
		const char* str = s.c_str();
		if(str[0] == '>')	n ++;
	}
	fp.close();
	return(n);
}
				
int seq::readseq(char **src_seq, char **src_name, int *len_seq, char *filename)
{
	int	i, k, n, len;
        char   *buf;

	std::ifstream fp(filename);
	if(!fp) { cout<<"open file error "<<filename; exit(0); }
        buf = (char *)smallapp::ckalloc(sizeof(char) * 1000000);

	n = 0;
	k = -1;
	std::string s;
	while(!fp.eof()) {
		std::getline(fp, s);
		const char *str = s.c_str();
		if(str[0] == '#')	continue;
		if(str[0] == '>')	{
			if(k >= 0)
                        {
                          src_seq[k] = (char *)smallapp::ckalloc(sizeof(char) * (n + 1));
                          memcpy((void *)src_seq[k], (void *)buf, sizeof(char) * (n + 1));
                          len_seq[k] = n;
                        }
			n = 0;
			k ++;
			sscanf(&str[1], "%s", buf);
                        src_name[k] = (char *)smallapp::ckalloc(sizeof(char) * (strlen(buf) + 1));
                        strcpy(src_name[k], buf);
		} else {
			for(i = 0; i < s.length(); i ++)	{
				if(str[i] >= 'a' && str[i] <= 'z') {
					buf[++ n] = char2int(str[i]);
				} else if(str[i] >= 'A' && str[i] <= 'Z') {
					buf[++ n] = char2int(str[i]
						 - 'A' + 'a');
				}
			}
		}
	}
	fp.close();
        
        if (k >= 0)
        {
          src_seq[k] = (char *)smallapp::ckalloc(sizeof(char) * (n + 1));
          memcpy((void *)src_seq[k], (void *)buf, sizeof(char) * (n + 1));
          len_seq[k] = n;
        }
	k ++;

        free((void *)buf);
	return(k);
}

char seq::aa2int(char c)
{
	int	i, k;

	for(i = 0; i < total_aa; i ++)	{
		if(c == aalist[i] || c == aalist[i] - 'A' + 'a')	{
			break;
		}
	}
	int	idum = 1234; //to be checked YY
	if(i == total_aa)	{
		printf("Not found %c\n", c);
		i = smallapp::ran_number(20, &idum);
	}
	return(i);
}

char seq::char2int(char c)
{
	int	i, k;

	for(i = 0; i < total_nuc; i ++)		{
		if(c == na_name[i])	{
			k = i;
			break;
		}
	}

	int	idum = 1234; //to be checked YY
	if(i == total_nuc)	{
		printf("Not found %c\n", c);
		k = smallapp::ran_number(4, &idum);
	}

	if(k > 3)	{
		k = smallapp::ran_number(4, &idum);
	}

	return(k);
}

char seq::char2intgen(char c)
{
	int	i, k;

	for(i = 0; i < total_nuc; i ++)		{
		if(c == na_name[i])	{
			k = i;
			break;
		}
	}

	int idum = 1234; //to-be-checked by YY
	if(i == total_nuc)	{
		printf("Not found %c\n", c);
		k = smallapp::ran_number(4, &idum);
	}

	return(k);
}


int seq::outputseq(ofstream &fp, char *seq, char *name, int pbeg, int pend)
{
	int i, j, l;
	char str[1000];

	sprintf(str, ">%s-%d-%d %d\n", name, pbeg, pend, pend-pbeg+1);
	fp<<str;
	l = 0; 
	for(i = pbeg; i <= pend; i ++)	{
		fp<<na_name[seq[i]];
		l ++;
		if(l == 60)	{
			l = 0;
			fp<<endl;
		}
	}
	if(l != 60)	fp<<endl;
	return(1);
}

int seq::writeseq(ofstream &fp, char *seq, char *name, int len_seq)
{
	int i, j, l;

	fp<<">"<<name<<endl;
	l = 0; 
	for(i = 0; i < len_seq; i ++)	{
		fp<<na_name[seq[i]];
		l ++;
		if(l == 60)	{
			l = 0;
			fp<<endl;
		}
	}
	if(l != 60)	fp<<endl;
	return(1);
}
