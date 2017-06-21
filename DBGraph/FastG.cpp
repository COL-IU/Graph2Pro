//Package: tag
//Program Name: FastG (for loading input file in FastG format)
//Latest modification on April 22th, 2016
//Yuzhen Ye <yye@indiana.edu> and Haixu Tang <hatang@indiana.edu>
//Affiliation: School of Informatics and Computing, Indiana University, Bloomington
#include "FastG.h"
#include "soap.h"
#include "smallapp.h"
#include "seq.h"

int FastG::countedge(char *edgeseqfile)
{
	int	num_edge = 0;
	std::string s;
	std::ifstream fp_edge(edgeseqfile);

        if(!fp_edge) {
                cout<<"cannot open file "<<edgeseqfile<<endl;
                exit(0);
        }
	while(!fp_edge.eof()) {
		std::getline(fp_edge, s);
		const char* str = s.c_str();
		if(str[0] == '>')	num_edge ++;
	}
	fp_edge.close();
	return(num_edge);
}

int FastG::readedge(char *edgeseqfile, EDGE *edge, VERTEX *vertex, int *edgeindex, int num_edge, int kmersize)
{
	int	i, j, k, l, n, m;
	int	num_vertex;
	uint256_t	c1, c2;
	char	*seq;
	int	seqlen, max_len;
	int	len, strand;
	char	node1[100], node2[100];
	char	tmpstr[5001], edgename[501], prevedgename[501];
	char	*tmpseq;
	std::string s;

	cout<<"edgeseqfile "<<edgeseqfile<<endl;
	edgename[0] = '\0';
	max_len = 0;

	tmpseq = (char *) smallapp::ckalloc(MAX_CON_LEN * sizeof(char));
	seqlen = 0;

	std::ifstream fp_edge(edgeseqfile);
	if(!fp_edge) {
		cout<<"cannot open file "<<edgeseqfile<<endl;
		exit(0);
	}
	n = -1;
	while(!fp_edge.eof()) {
		std::getline(fp_edge, s);
		if(s.empty()) continue;
		const char* str = s.c_str();
/*
printf("str %s\n", str);
getchar();
*/
		if(str[0] == '>')	{
			if(n >= 0)	{
				edge[n].seq = (char *) smallapp::ckalloc(seqlen * sizeof(char));
				for(j = 0; j < seqlen; j ++)	{
					edge[n].seq[j] = tmpseq[j];
				}
				edge[n].length = seqlen;
				if(edge[n].length > max_len)	max_len = edge[n].length;
				seqlen = 0;
			}
			replace(s.begin(), s.end(), ':', ' ');
			strcpy(prevedgename, edgename);
			sscanf(&str[1], "%s%s", edgename, tmpstr);
			sscanf(&edgename[5], "%d", &j);
			j --;
			l = strlen(edgename);
			n ++;
			if(edgename[l - 1] == ';')	{
				edgename[l - 1] = '\0';
				l --;
			}
			if(edgename[l - 1] == '\'')	{
				edgeindex[2 * j + 1] = n;
			} else	{
				edgeindex[2 * j] = n;
			}
/*
if(j == 241)	{
printf("str %s\n", str);
printf("edgename %s j %d n %d\n", edgename, j, n);
getchar();
}
j = 241;
printf("j %d edgeindex %d %d\n", j, edgeindex[2 * j], edgeindex[2 * j + 1]);
*/
			// verify palindromic edges have no reverse complementary edges in fastG file
			if(n > 0)	{
				if(!strncmp(edgename, prevedgename, l - 1) && edgename[l - 1] == '\'')	{		// reverse complement of the previous edge
					edge[n].bal_edge = n - 1;
					edge[n - 1].bal_edge = n;
				} else if(edgename[l - 1] != '\'' && n > 0 && prevedgename[l - 1] != '\'')	{
					edge[n - 1].bal_edge = n - 1;	// palindromic edge
				}
			}
		} else	{
			// get edge sequnece
			l = strlen(str);
			for(i = 0; i < l; i ++)	{
				if(str[i] >= 'a' && str[i] <= 'z')	{
					tmpseq[seqlen ++] = seq::char2int(str[i]);
				} else if(str[i] >= 'A' && str[i] <= 'Z')	{
					tmpseq[seqlen ++] = seq::char2int(str[i] - 'A' + 'a');
				}
			}
		}
	}
	if(seqlen > 0)	{
		edge[n].seq = (char *) smallapp::ckalloc(seqlen * sizeof(char));
		for(j = 0; j < seqlen; j ++)	{
			edge[n].seq[j] = tmpseq[j];
		}
		edge[n].length = seqlen;
		if(edge[n].length > max_len)	max_len = edge[n].length;
		n ++;
	}

	free((void *) tmpseq);
	fp_edge.close();
	if(n != num_edge)	{
		printf("Num-of-edges not equal: %d %d\n", n, num_edge);
		exit(0);
	}
	printf("%d edges input\n", n);

	int startedge, outedge, outvertex;
	num_vertex = 0;
	std::ifstream fp_edge2(edgeseqfile);
	if(!fp_edge2) {
		cout<<"cannot open file "<<edgeseqfile<<endl;
		exit(0);
	}
	n = 0;
	while(!fp_edge2.eof()) {
		std::getline(fp_edge2, s);
		if(s.empty()) continue;
		const char* str = s.c_str();
		if(str[0] == '>')	{
			replace(s.begin(), s.end(), ':', ' ');
			sscanf(&str[1], "%s%s", edgename, tmpstr);
			sscanf(&edgename[5], "%d", &j);
			j --;
			l = strlen(edgename);
			if(edgename[l - 1] == ';')	{
				tmpstr[0] = '\0';	// reset tmpstr if there is no outgoing edges
				edgename[l - 1] = '\0';
				l --;
			}
			if(edgename[l - 1] == '\'')	{
				startedge = edgeindex[2 * j + 1];
			} else	{
				startedge = edgeindex[2 * j];
			}
			outvertex = -1;
			l = strlen(tmpstr);
			k = 0;
			while(k < l - 6)	{
				outedge = get_next_edge(tmpstr, &k, l, edgeindex);
				k ++;
				if(edge[outedge].nodeindex[0] >= 0)	{
					outvertex = edge[outedge].nodeindex[0];
					break;
				}
			}
			if(outvertex < 0)	{
				outvertex = num_vertex;
				vertex[num_vertex].lindex = str2code(&(edge[startedge].seq[edge[startedge].length - kmersize]), kmersize);
				num_vertex ++;
			}
			edge[startedge].nodeindex[1] = outvertex;
			vertex[outvertex].lastedge = soap::add_edge_to_vertex(&edge[startedge], vertex[outvertex].lastedge);
			vertex[outvertex].indegree ++;
			l = strlen(tmpstr);
			k = 0;
			while(k < l - 6)	{
				outedge = get_next_edge(tmpstr, &k, l, edgeindex);
				k ++;
				if(edge[outedge].nodeindex[0] < 0)	{
					edge[outedge].nodeindex[0] = outvertex;
					vertex[outvertex].nextedge = soap::add_edge_to_vertex(&edge[outedge], vertex[outvertex].nextedge);
					vertex[outvertex].outdegree ++;
				}
			}
		}
	}
	fp_edge2.close();
	printf("num_vertex %d\n", num_vertex);

	std::ifstream fp_edge3(edgeseqfile);
	if(!fp_edge3) {
		cout<<"cannot open file "<<edgeseqfile<<endl;
		exit(0);
	}
	n = 0;
	while(!fp_edge3.eof()) {
		std::getline(fp_edge3, s);
		if(s.empty()) continue;
		const char* str = s.c_str();
		if(str[0] == '>')	{
			replace(s.begin(), s.end(), ':', ' ');
			sscanf(&str[1], "%s%s", edgename, tmpstr);
			sscanf(&edgename[5], "%d", &j);
			j --;
			l = strlen(edgename);
			if(edgename[l - 1] == ';')	{
				tmpstr[0] = '\0';	// reset tmpstr if there is no outgoing edges
				edgename[l - 1] = '\0';
				l --;
			}
			if(edgename[l - 1] == '\'')	{
				startedge = edgeindex[2 * j + 1];
			} else	{
				startedge = edgeindex[2 * j];
			}
			if(edge[startedge].nodeindex[0] < 0)	{
				vertex[num_vertex].lindex = str2code(edge[startedge].seq, kmersize);
				edge[startedge].nodeindex[0] = num_vertex;
				vertex[num_vertex].nextedge = soap::add_edge_to_vertex(&edge[startedge], vertex[num_vertex].nextedge);
				vertex[num_vertex].outdegree ++;
				num_vertex ++;
			}
		}
	}
	fp_edge3.close();
	printf("num_vertex %d\n", num_vertex);

/*
	for(i = 0; i < num_vertex; i ++)	{
		if(vertex[i].outdegree == 0)	{
			printf("i %d degrees %d %d\n", i, vertex[i].indegree, vertex[i].outdegree);
			getchar();
		}
	}

for(i = 0; i < num_edge; i ++)	{
	printf("edge %d length %d start %d outdegree %d end %d indegree %d\n", i, edge[i].length, edge[i].nodeindex[0], vertex[edge[i].nodeindex[0]].outdegree,edge[i].nodeindex[1], vertex[edge[i].nodeindex[1]].indegree);
getchar();
}
*/

	return(num_vertex);
}

int FastG::get_next_edge(char *tmpstr, int *start, int end, int *edgeindex)
{
	int	i, j, k;

	i = *start;
	if(i < end - 6)	{
		sscanf(&tmpstr[i + 5], "%d", &k); 
		k --;
	}
	for(; *start < end; (*start) ++)	{
		if(tmpstr[*start] == ',' || tmpstr[*start] == ';')	{
			break;
		}
	}
	if(tmpstr[*start - 1] == '\'')	{
//printf("k %d edgeindex %d\n", k, edgeindex[2 * k + 1]);
		return(edgeindex[2 * k + 1]);
	} else	{
		return(edgeindex[2 * k]);
	}
}

uint256_t FastG::str2code(char *seq, int kmersize)
{
	int	i, k, l;
	uint256_t	code;

	code.bits[0] = code.bits[1] = 0;
	k = 0;
	for(i = 0; i < kmersize; i ++)	{
		code.bits[k] *= 4;
		code.bits[k] += seq[i];
		if(i == 64)	{
			k ++;
		} 
	}
	return(code);
}
