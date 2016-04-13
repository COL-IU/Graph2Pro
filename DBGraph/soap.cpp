#include "soap.h"
#include "smallapp.h"
#include "seq.h"

int soap::countedge(char *edgefile)
{
	int	num_edge;
	std::ifstream fp_edge(edgefile);
        if(!fp_edge) {
                cout<<"cannot open file "<<edgefile<<endl;
                exit(0);
        }
	std::string s;
	std::getline(fp_edge, s);
	sscanf(s.c_str(), "%*s%d", &num_edge);
	fp_edge.close();
	return(num_edge);
}

int soap::readedge(bool SOAP2, char *edgefile, char *edgeseqfile, EDGE *edge, VERTEX *vertex, int num_vertex, int kmersize)
{
	int	i, j, k, l, n, m;
	__int128_t	c1, c2;
	int	num_edge;
	char	*seq;
	int	seqlen, max_len;
	int	len, strand;
	int	edgeindex;
	char	node1[100], node2[100];
	char	tmpstr[501];

	cout<<"edgefile "<<edgefile<<endl;
	std::ifstream fp_edge(edgefile);
	if(!fp_edge) {
		cout<<"cannot open file "<<edgefile<<endl;
		exit(0);
	}
	std::string s;
	std::getline(fp_edge, s);
	sscanf(s.c_str(), "%*s%d", &num_edge);

	n = 0;
	max_len = 0;
	while(!fp_edge.eof()) {
		std::getline(fp_edge, s);
		if(s.empty()) continue;
		if(SOAP2)	{	/* SOAP-denovo2 format	*/
			replace(s.begin(), s.end(), ',', ' ');
			sscanf(s.c_str(), "%*s%d%d%*s%*s%s%*s%s", &len, &strand, node1, node2);
		} else	{	/* SOAP-denovo format */
			sscanf(s.c_str(), "%*s%s", tmpstr);
			for(i = 0; i < strlen(tmpstr) - 1; i ++)	{
				if(tmpstr[i] == ',')	{
					tmpstr[i] = ' ';
				}
			}
			sscanf(tmpstr, "%d%s%s%d%*d", &len, node1, node2, &strand);
			len += kmersize;	/*	SOAP-denovo format	*/
		}
		c1 = code2kmer(node1);
		c2 = code2kmer(node2);
		insert_edge(&edge[n], c1, c2, vertex, num_vertex);
		edge[n].length = len;
		edge[n].seq = new char[edge[n].length];
		if(edge[n].length > max_len)	max_len = edge[n].length;
		if(strand == 0)	{	// 1: forward edge; 0: palindromic edge (its revserse complement is itself)
			edge[n].bal_edge = n;
		} else if(strand == -1) {	// -1: reverse complement of the previous edge
			edge[n].bal_edge = n - 1;
			edge[n - 1].bal_edge = n;
		}
		n ++;
	}
	fp_edge.close();
	if(n != num_edge)	{
		printf("Num-of-edges not equal: %d %d\n", n, num_edge);
		exit(0);
	}

	seq = new char[max_len + 1000];
	std::ifstream fp_edgeseq(edgeseqfile);
	if(!fp_edgeseq) {
		cout<<"open file error: "<<edgeseqfile<<endl;
		exit(1);
	}
	n = 0; 
	seqlen = 0;
	while(!fp_edgeseq.eof()) {
		std::getline(fp_edgeseq, s);
		const char* str = s.c_str();
		if(str[0] == '>')	{
			if(n > 0)	{
				if(seqlen != edge[edgeindex].length)	{
					printf("Length not match for edge %d: %d vs %d\n", edgeindex, seqlen, edge[edgeindex].length);
					exit(0);
				}
				for(i = 0; i < seqlen; i ++)	{
					edge[edgeindex].seq[i] = seq[i];
					//cout<<na_name[seq[i]];
				}
				edge[edgeindex].tag = 1;
				if(edge[edgeindex].bal_edge != edgeindex)	{
					for(i = 0; i < seqlen; i ++)	{
						edge[edge[edgeindex].bal_edge].seq[i] = rev(seq[seqlen - 1 - i]);
					}
					edge[edge[edgeindex].bal_edge].tag = 1;
					n ++;
				}
			}
			sscanf(&str[1], "%d", &edgeindex); //edgeindex transformed from read name????
			edgeindex = edgeindex - 1;
			seqlen = 0;
			n ++;
		} else	{
			l = strlen(str);
			for(i = 0; i < l; i ++)    {
				if(str[i] >= 'a' && str[i] <= 'z')	{
					seq[seqlen ++] = seq::char2int(str[i]);
				} else if(str[i] >= 'A' && str[i] <= 'Z')	{
					seq[seqlen ++] = seq::char2int(str[i] - 'A' + 'a');
				}
			}
		}
	}
	if(seqlen != edge[edgeindex].length)	{
		printf("Length not match for edge %d: %d vs %d\n", edgeindex, seqlen, edge[edgeindex].length);
		exit(0);
	}
	for(i = 0; i < seqlen; i ++)	{
		edge[edgeindex].seq[i] = seq[i];
	}
	edge[edgeindex].tag = 1;
	if(edge[edgeindex].bal_edge != edgeindex)	{
		for(i = 0; i < seqlen; i ++)	{
			edge[edge[edgeindex].bal_edge].seq[i] = rev(seq[seqlen - 1 - i]);
		}
		edge[edge[edgeindex].bal_edge].tag = 1;
		n ++;
	}
	fp_edgeseq.close();
	delete[] seq;
	printf("%d edge sequences included.\n", n);
	for(i = 0; i < num_edge; i ++)	{
		if(edge[i].tag == 1)	{
			link_vertex(&edge[i], vertex, num_vertex);
		}
	}
	return(num_edge);
}

int soap::insert_edge(EDGE *edge, __int128_t code1, __int128_t code2, VERTEX *vertex, int num_vertex)
{
	int	i, j, k;

	edge -> nodeindex[0] = loc_vertex(code1, vertex, num_vertex);
	edge -> nodeindex[1] = loc_vertex(code2, vertex, num_vertex);
}

int soap::link_vertex(EDGE *edge, VERTEX *vertex, int num_vertex)
{
	vertex[edge -> nodeindex[0]].nextedge = add_edge_to_vertex(edge, vertex[edge -> nodeindex[0]].nextedge);
	vertex[edge -> nodeindex[0]].outdegree ++;
	vertex[edge -> nodeindex[1]].lastedge = add_edge_to_vertex(edge, vertex[edge -> nodeindex[1]].lastedge);
	vertex[edge -> nodeindex[1]].indegree ++;
}

//binary search to locate vertex for given code
int soap::loc_vertex(__int128_t code, VERTEX *vertex, int num_vertex)
{
	int	i, j, k, l, n;
	char	str[200];

	if(num_vertex == 0)	{
		index2code(code, str);
		printf("Vertex not found: code %s\n", str);
		exit(0);
	}
	n = num_vertex / 2;
	if(vertex[n].index == code)	{
		return(n);
	} else if(vertex[n].index < code)	{
		return(n + 1 + loc_vertex(code, &vertex[n + 1], num_vertex - n - 1));
	} else	{
		return(loc_vertex(code, vertex, n));
	}
}

EDGELIST* soap::add_edge_to_vertex(EDGE *edge, EDGELIST *edgelist0)
{
	int	i, j, k, l;
	EDGELIST	*edgelist;

	edgelist = new EDGELIST[1];
	edgelist -> edge = edge;
	edgelist -> next = edgelist0;
	return(edgelist);
}

void* soap::delete_edgelist(EDGELIST *edgelist)
{
	EDGELIST *edgelist0;

	if(edgelist)	{
		edgelist0 = edgelist -> next;
		delete edgelist;
		delete_edgelist(edgelist0);
	}
}

int soap::readvertex(bool SOAP2_tag, char *edgefile, VERTEX *vertex)
{
	int	num_vertex;
	char	node1[501], node2[501], node[501];
	int	i, j, k, l;
	char	tmpstr[501];
	std::string s;

	num_vertex = 0;
	std::ifstream fp(edgefile);
	if(!fp) {
		cout<<"open file error "<<edgefile<<endl;
		exit(1);
	}
	std::getline(fp, s);
	while(!fp.eof()) {
		std::getline(fp, s);	
		if(SOAP2_tag)	{	/* SOAP-denovo2 format */
			replace(s.begin(), s.end(), ',', ' ');
			sscanf(s.c_str(), "%*s%*s%*s%*s%*s%s%*s%s", node1, node2);
			vertex[num_vertex].index = code2kmer(node1);
			num_vertex ++;
			vertex[num_vertex].index = code2kmer(node2);
			num_vertex ++;
		} else	{	/* SOAP-denovo format */
			sscanf(s.c_str(), "%*s%s", tmpstr);
			l = strlen(tmpstr);
			for(i = 0; i < l - 1; i ++)	{
				if(tmpstr[i] == ',')	{
					tmpstr[i] = ' ';
				}
			}
			sscanf(tmpstr, "%*d%s%s%*d%*d", node1, node2);
			vertex[num_vertex].index = code2kmer(node1);
			num_vertex ++;
			vertex[num_vertex].index = code2kmer(node2);
			num_vertex ++;
		}
	}
	fp.close();

	std::qsort(vertex, num_vertex, sizeof(VERTEX), smallapp::comparvertex);

/*	Remove duplicated nodes	*/

	k = 1;
	for(i = 1; i < num_vertex; i ++)	{
		if(vertex[i].index != vertex[i - 1].index)	{
			if(i > k)	{
				vertex[k ++] = vertex[i];
			} else	{
				k ++;
			}
		}
	}
	num_vertex = k;

	return(num_vertex);
}

int soap::printcode(char *str, char *printstr)
{
	int	i, j, k, l, c1, c2, code;

	l = strlen(str);

	for(i = l - 1; i >= 0; i --)	{
		if(str[i] >= '0' && str[i] <= '9')	{
			code = str[i] - '0';
		} else if(str[i] >= 'a' && str[i] <= 'f')	{
			code = 10 + str[i] - 'a';
		}
		c1 = code / 4;
		c2 = code - c1 * 4;
		printstr[2 * i + 1] = code_name[c2];
		printstr[2 * i] = code_name[c1];
	}
	printstr[2 * l] = '\0';
	return(2 * l);
}

__int128_t soap::code2kmer(char *str)
{
	int	i, k, l, code;;
	__int128_t base;
	__int128_t n;

	l = strlen(str);
	k = 0;
	n = 0;
	for(i = l - 1; i >= 0; i --)	{
		if(str[i] >= '0' && str[i] <= '9')	{
			code = str[i] - '0';
		} else if(str[i] >= 'a' && str[i] <= 'f')	{
			code = 10 + str[i] - 'a';
		}
		base = code;
		base = base << k;
		n += base;
		k += 4;
	}
	return(n);
}

int soap::index2code(__int128_t code, char *str)
{
	int	i, k, l, n;
	char	tmpstr[200];

	i = 0;
	do	{
		k = code & 15; 
		if(k >= 0 && k <= 9)	{
			tmpstr[i] = k + '0';
		} else {
			tmpstr[i] = k - 10 + 'a';
		}
		i ++;
		code = code >> 4;
	} while(code != 0);
	for(n = 0; n < i; n ++)	{
		str[n] = tmpstr[i - 1 - n];
	}
	str[i] = '\0';

	return(i);
}

int soap::encodekmer(char *printstr, int length, char *str)
{
	int	i, k, l, n, code, c1, c2;
	char	*tmpstr;

	tmpstr = new char[length];
	tmpstr[0] = encodechar(printstr[0]);
	n = 1;
	for(i = 1; i < length; i += 2)	{
		c1 = encodechar(printstr[i]);
		c2 = encodechar(printstr[i + 1]);
		code = c1 * 4 + c2;
		if(code >= 0 && code <= 9)	{
			tmpstr[n] = code + '0';
		} else {
			tmpstr[n] = code - 10 + 'a';
		}
		n ++;
	}
	for(i = 0; i < n; i ++)	{
		if(tmpstr[i] != '0')	{
			break;
		}
	}
	k = 0;
	for(; i < n; i ++)	{
		str[k ++] = tmpstr[i];
	}
	str[k] = '\0';
	delete[] tmpstr;
	return(k);
}

char soap::encodechar(char c)
{
	/* codes: A -> 0, C -> 1, T -> 2, G -> 3	*/
	if(c == 'A' || c == 'a')	{
		return(0);
	} else if(c == 'C' || c == 'c')	{
		return(1);
	} else if(c == 'T' || c == 't')	{
		return(2);
	} else if(c == 'G' || c == 'g')	{
		return(3);
	} else	{
		printf("Character unfound: %c\n", c);
		exit(0);
	}
}
