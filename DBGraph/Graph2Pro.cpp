#include "GraphTrav.h"
#include "smallapp.h"
#include "soap.h"

void GraphTrav::Graph2Pro(void)
{
	int	i, j, k, l, n;
	int	num_part;
	int	pos[2];
	char	*start_edge;
	int	*edge_offset;
	int	offset;
	int	initedge_index, initedge_offset;
	char	*label;
	char	*initpep;
	int	len_initpep;

	cout<<"\n>>>Now Traverse the graph..."<<endl;
	//i = 610685;
	printf("Check at the beginning of graph traversal\n");
	printf("Peptides of lengths between %d and %d will be output\n", peptide_min, peptide_max);
        //printf("Edge i = %d tag %d, begtag %d endtag %d length %d pathseg-length %d\n", i, edge[i].tag, edge[i].begtag, edge[i].endtag, edge[i].length, pathlist[i] -> length);

	max_pep = 100 * num_edge;

	pepseq = NULL;

	start_edge = new char[num_edge];
	edge_offset = new int[num_edge * 3];

	for(i = 0; i < num_edge; i ++)	{
		edge[i].tag = edge[i].multitag = 0; //inilize the tag -- set to 0 
		start_edge[i] = 0;
	}

	for(i = 0; i < num_edge; i ++)	{
		if(vertex[edge[i].nodeindex[0]].indegree == 0 && edge[i].length > 0)	{
			start_edge[i] = 1;
		} else	{
			for(offset = 0; offset < 3; offset ++)	{
				edge_offset[3 * i + offset] = offset + Loc_Stopcodon(&(edge[i].seq[kmersize + offset]), edge[i].length - kmersize - offset);
			}
		}
	}
	num_peptides = 0;

	initpep = new char[max_length];

	for(i = 0; i < num_edge; i ++)	{
/*
printf("edge %d %d %d start_edge %d num_peptides %d nodes %d %d\n", i, &edge[i], edge[i].length, start_edge[i], num_peptides,
	 edge[i].nodeindex[0], edge[i].nodeindex[1]);
*/
/*
if(i > 1613)	{
	getchar();
}
*/
		len_initpep = 0;
		initedge_index = i;
		if(start_edge[i] == 1 && edge[i].tag == 0)	{
			for(offset = 0; offset < 3; offset ++)	{
/*
if(i == 1581451)	{
printf("pos1 i %d offset %d num_peptides %d\n", i, offset, num_peptides);
getchar();
}
*/
				initedge_offset = offset;
				traverse_graph(&edge[i], offset, initpep, len_initpep, initedge_index, initedge_offset);
			}
		} else	{
			for(offset = 0; offset < 3; offset ++)	{
				if(edge_offset[3 * i + offset] > offset)	{
					initedge_offset = kmersize + edge_offset[3 * i + offset];
/*
if(i == 1581451)	{
printf("pos2 offset %d %d num_peptides %dl length %d initedge_offset %d\n", offset, edge_offset[3 * i + offset], num_peptides, edge[i].length, initedge_offset);
getchar();
}
*/
					traverse_graph(&edge[i], initedge_offset, initpep, len_initpep, initedge_index, initedge_offset);
				}
/*
if(i == 747261)	{
printf("after offset %d\n", offset);
printf("continue\n");
}
*/
			}
		}
	}

	delete[] start_edge;
	delete[] edge_offset;
	delete[] initpep;
}

int GraphTrav::Loc_Stopcodon(char *seq, int length)
{
	int	i, l;
	char	*pepseq2;

	pepseq2 = new char[12];

	for(i = 0; i < length - 2; i += 3)	{
		l = trans_codon2aa(&seq[i], 3, pepseq2);
		if(pepseq2[0] == 20 || ckcutsite(pepseq2[0]) == 1)	{	// stop codon or cleavage site
			if(i < length - 3)	{
				return(i + 3);
			}
		}
	}

	delete[] pepseq2;
	return(0);
}

char GraphTrav::ckcutsite(char seq)
{
	int	i;
	for(i = 0; i < num_cutsite; i ++)	{
		if(seq == cutsite[i])	{
			return(1);
		}
	}
	return(0);
}

void GraphTrav::traverse_graph(EDGE *inedge, int offset, char *lastpep, int len_lastpep, int initedge_index, int initedge_offset)
{
	int	i, j, k, l;
	int	leftover;
	int	len_tmplast=0, len_tmplast2 = 0;
	char	tmpcodon[10], tmppep[3];
	char	midtag;
	char	flag_stopcodon;
	char	*tmplastpep, *tmplastpep2;
	int	tmp_initedge_index, tmp_initedge_offset;
	int	tmp_initedge_index2, tmp_initedge_offset2;
	int	end_edge_index, end_edge_offset;
	int	ncut, ncut_cpy, ncut_cpy2;
	char	begtag, tmpbegtag, tmpbegtag2;
	EDGELIST	*nextedge, *nextnextedge;

	tmplastpep = new char[max_length];
	tmplastpep2 = new char[max_length];

	k = 1;
	k = k << (offset % 3);
	if(inedge -> tag & k == 0)	{
		return;			// The edge has been traversed.
	}

	inedge -> tag |= k;	// mark that the edge is visited at offset

	depth = 0;
	flag_stopcodon = 0;

	len_lastpep = 0;
	ncut = 0;
	begtag = 1;
	end_edge_index = initedge_index;
	end_edge_offset = initedge_offset;
	if(inedge -> length - offset >= 3)	{
		if(inedge -> length >= min_contig_leg)	{
			midtag = extract_peptides(&(inedge -> seq[offset]), inedge -> length - offset, 1, lastpep, &len_lastpep,
				 &leftover, &ncut, begtag, initedge_index, &initedge_offset, end_edge_index, &end_edge_offset);
		} else	{
			midtag = extract_peptides(&(inedge -> seq[offset]), inedge -> length - offset, 0, lastpep, &len_lastpep,
				 &leftover, &ncut, begtag, initedge_index, &initedge_offset, end_edge_index, &end_edge_offset);
		}
	} else	{
		leftover = inedge -> length - offset;
	}

	len_tmplast = 0;	// temporarily store the lastpep sequence
	len_tmplast = appendpep(tmplastpep, len_tmplast, lastpep, len_lastpep);
	tmp_initedge_index = initedge_index;
	tmp_initedge_offset = initedge_offset;
	ncut_cpy2 = ncut;
	tmpbegtag = begtag;
	nextedge = vertex[inedge -> nodeindex[1]].nextedge;
	while(nextedge)	{
/*
printf("nextedge %d length %d nodes %d %d multitag %d\n", nextedge -> edge, nextedge -> edge -> length, nextedge -> edge -> nodeindex[0], nextedge -> edge -> nodeindex[1], nextedge -> edge -> multitag);
printf("leftover %d offset %d\n", leftover, offset);
*/
		if(nextedge -> edge -> multitag == 1)	{
			nextedge = nextedge -> next;
			continue;
		}
		len_lastpep = 0;	// restore the lastpep sequence
		len_lastpep = appendpep(lastpep, len_lastpep, tmplastpep, len_tmplast);
		initedge_index = tmp_initedge_index;
		initedge_offset = tmp_initedge_offset;
		ncut = ncut_cpy;
		begtag = tmpbegtag;
		if(leftover == 0)	{
			nextedge -> edge -> multitag = 1;
			if(len_lastpep <= peptide_max && depth <= max_depth)	{
/*
printf("A enter continue nextedge %d length %d\n", nextedge -> edge, nextedge -> edge -> length);
*/
				depth ++;
				if(len_lastpep == 0)	{
					initedge_index = nextedge -> edge -> index;
					initedge_offset = kmersize;
				}
				end_edge_index = nextedge -> edge -> index;
				end_edge_offset = kmersize;
				continue_extract_peptides(nextedge -> edge, 0, lastpep, &len_lastpep, &leftover, &ncut, begtag,
					 initedge_index, initedge_offset, end_edge_index, end_edge_offset);
				depth --;
			}
			nextedge -> edge -> multitag = 0;
		} else	{
			for(i = 0; i < leftover; i ++)	{
				tmpcodon[i] = inedge -> seq[inedge -> length - leftover + i];
			}
			if(nextedge -> edge -> length - kmersize <= 0)	{
				printf("edge length error: %d, kmersize %d\n", nextedge -> edge -> length, kmersize);
				exit(-1);
			}
			if(nextedge -> edge -> length - kmersize + leftover < 3)	{
				tmpcodon[1] = nextedge -> edge -> seq[kmersize];
				len_tmplast2 = 0;	// temporarily store the lastpep sequence
				len_tmplast2 = appendpep(tmplastpep2, len_tmplast2, lastpep, len_lastpep);
				tmp_initedge_index2 = initedge_index;
				tmp_initedge_offset2 = initedge_offset;
				ncut_cpy2 = ncut;
				tmpbegtag2 = begtag;
				nextnextedge = vertex[nextedge -> edge -> nodeindex[1]].nextedge;
				while(nextnextedge)	{
/*
printf("nextnextedge %d length %d nodes %d %d multitag %d\n", nextnextedge -> edge, nextnextedge -> edge -> length, nextnextedge -> edge -> nodeindex[0], nextnextedge -> edge -> nodeindex[1], nextnextedge -> edge -> multitag);
*/
					if(nextnextedge -> edge -> multitag == 1)	{
						nextnextedge = nextnextedge -> next;
						continue;
					}
					len_lastpep = 0;	// restore the lastpep sequence
					len_lastpep = appendpep(lastpep, len_lastpep, tmplastpep2, len_tmplast2);
					initedge_index = tmp_initedge_index2;
					initedge_offset = tmp_initedge_offset2;
					end_edge_index = nextnextedge -> edge -> index;
					end_edge_offset = kmersize + 1;
					ncut = ncut_cpy;
					begtag = tmpbegtag;
					tmpcodon[2] = nextnextedge -> edge -> seq[kmersize];
					l = trans_codon2aa(tmpcodon, 3, tmppep);
					if(ckcutsite(tmppep[0]) == 1)	ncut ++;
					if(tmppep[0] == 20 || ncut > max_miscleavage)	{	// stop codon or cleavage site
						if(tmppep[0] != 20)	{
							len_lastpep = appendpep(lastpep, len_lastpep, tmppep, l);
							flag_stopcodon = 0;
						} else	{
							flag_stopcodon = 1;
						}
/*
for(j = 0; j < length; j ++)
	printf("%c", na_name[seq[j]]);
printf("\n");
printf("pos3 length %d len_pep %d edges %d %d %d %d\n", nextnextedge -> edge -> length, len_lastpep, initedge_index, initedge_offset, end_edge_index, end_edge_offset);
getchar();
*/
						num_peptides = copypepseq(lastpep, len_lastpep, begtag, initedge_index, initedge_offset, 
								end_edge_index, end_edge_offset, flag_stopcodon);
						len_lastpep = 0;
						initedge_index = end_edge_index;
						initedge_offset = end_edge_offset;
						if(tmppep[0] != 20)	{
							begtag = 1;
						} else	{
							begtag = 0;
						}
						ncut = 0;
					} else	{
						len_lastpep = appendpep(lastpep, len_lastpep, tmppep, l);
					}
					nextedge -> edge -> multitag = 1;
					nextnextedge -> edge -> multitag = 1;
					if(len_lastpep <= peptide_max && depth <= max_depth)	{
						depth ++;
						continue_extract_peptides(nextnextedge -> edge, 1, lastpep, &len_lastpep, &leftover,
							 &ncut, begtag, initedge_index, initedge_offset, end_edge_index, end_edge_offset);
						depth --;
					}
					nextedge -> edge -> multitag = 0;
					nextnextedge -> edge -> multitag = 0;
					nextnextedge = nextnextedge -> next;
				}
			} else	{
				for(j = 0; i < 3; i ++, j ++)	{
					tmpcodon[i] = nextedge -> edge -> seq[kmersize + j];
				}
				end_edge_index = nextedge -> edge -> index;
				end_edge_offset = kmersize + j;
				l = trans_codon2aa(tmpcodon, 3, tmppep);
				if(ckcutsite(tmppep[0]) == 1)	ncut ++;
				if(tmppep[0] == 20 || ncut > max_miscleavage)	{	// stop codon or cleavage site
					if(tmppep[0] != 20)	{
						len_lastpep = appendpep(lastpep, len_lastpep, tmppep, l);
						flag_stopcodon = 0;
					} else	{
						flag_stopcodon = 1;
					}
/*
if(nextedge -> edge -> index == 1563520)	{
printf("pos4 edge %d length %d len_pep %d edges %d %d %d %d\n", nextedge -> edge -> index, nextedge -> edge -> length, len_lastpep, initedge_index, initedge_offset, end_edge_index, end_edge_offset);
getchar();
}
*/
					num_peptides = copypepseq(lastpep, len_lastpep, begtag, initedge_index,
							 initedge_offset, end_edge_index, end_edge_offset, flag_stopcodon);
					len_lastpep = 0;
					initedge_index = end_edge_index;
					initedge_offset = end_edge_offset;
					if(tmppep[0] != 20)	{
						begtag = 1;
					} else	{
						begtag = 0;
					}
					ncut = 0;
				} else	{
					len_lastpep = appendpep(lastpep, len_lastpep, tmppep, l);
				}
				nextedge -> edge -> multitag = 1;
				if(len_lastpep <= peptide_max && depth <= max_depth)	{
					depth ++;
					continue_extract_peptides(nextedge -> edge, 3 - leftover, lastpep, &len_lastpep, &leftover,
						 &ncut, begtag, initedge_index, initedge_offset, end_edge_index, end_edge_offset);
					depth --;
				}
				nextedge -> edge -> multitag = 0;
			}
		}
		nextedge = nextedge -> next;
	}
/*
printf("exit traverse %d\n");
*/

	delete[] tmplastpep;
	delete[] tmplastpep2;
}

char GraphTrav::extract_peptides(char *seq, int length, char label, char *lastpep, int *len_lastpep, int *leftover, int *ncut, char begtag, int init_edge_index, int *init_edge_offset, int end_edge_index, int *end_edge_offset)
{
	int	i, j, k, l;
	char	midtag = 0;
	char	tmpcodon[10], tmppep[3];
	char	flag_stopcodon;
	int	offset;

	for(i = 0; i <= length - 3; i += 3)	{
		*end_edge_offset = (*end_edge_offset) + 3;
		l = trans_codon2aa(&seq[i], 3, tmppep);
		if(ckcutsite(tmppep[0]) == 1)	(*ncut) ++;
		if(tmppep[0] == 20 || *ncut > max_miscleavage)	{	// stop codon or cleavage site
			if(tmppep[0] != 20)	{
				*len_lastpep = appendpep(lastpep, *len_lastpep, tmppep, l);
				flag_stopcodon = 0;
			} else	{
				flag_stopcodon = 1;
			}
/*
if(init_edge_index == 1581451 || end_edge_index == 1581451)	{
for(j = 0; j < length; j ++)
	printf("%c", na_name[seq[j]]);
printf("\n");
for(j = 0; j < *len_lastpep; j ++)	{
	printf("%c", aalist[lastpep[j]]);
}
printf("\n");
printf("pos1 length %d len_pep %d edges %d %d %d %d\n", length, *len_lastpep, init_edge_index, *init_edge_offset, end_edge_index, *end_edge_offset);
getchar();
}
*/
			num_peptides = copypepseq(lastpep, *len_lastpep, begtag, init_edge_index, *init_edge_offset, end_edge_index, *end_edge_offset, flag_stopcodon);
			*len_lastpep = 0;
			*init_edge_offset = *end_edge_offset;
			
			*ncut = 0;
			if(label == 1)	{
				midtag = 1;
				break;
			}
		} else	{
			*len_lastpep = appendpep(lastpep, *len_lastpep, tmppep, l);
		}
	}
	if(midtag == 1)	{
		*leftover = length % 3;
		for(i = length - 3 - *leftover; i >= 0; i -= 3)	{
			l = trans_codon2aa(&seq[i], 3, tmppep);
			if(tmppep[0] == 20)	{
				break;
			}
		}
		offset = edge[init_edge_index].length - length;
		*len_lastpep = 0;
		*init_edge_offset = i + 3;
		for(i += 3; i <= length - *leftover - 3; i += 3)	{
			*end_edge_offset = i + 3;
			l = trans_codon2aa(&seq[i], 3, tmppep);
			if(ckcutsite(tmppep[0]) == 1)	(*ncut) ++;
			if(tmppep[0] == 20 || *ncut > max_miscleavage)	{	// stop codon or cleavage site
				if(tmppep[0] != 20)	{
					*len_lastpep = appendpep(lastpep, *len_lastpep, tmppep, l);
					flag_stopcodon = 0;
				} else	{
					flag_stopcodon = 1;
				}
/*
if(init_edge_index == 1581451 || end_edge_index == 1581451)	{
for(j = 0; j < length; j ++)
	printf("%c", na_name[seq[j]]);
printf("\n");
for(j = 0; j < *len_lastpep; j ++)	{
	printf("%c", aalist[lastpep[j]]);
}
printf("\n");
printf("pos5 offset %d length %d len_pep %d edges %d %d %d %d\n", offset, length, *len_lastpep, init_edge_index, *init_edge_offset + offset, end_edge_index, *end_edge_offset + offset);
getchar();
}
*/
				num_peptides = copypepseq(lastpep, *len_lastpep, begtag, init_edge_index, *init_edge_offset + offset, end_edge_index, *end_edge_offset + offset, flag_stopcodon);
				*len_lastpep = 0;
				*init_edge_offset = *end_edge_offset;
				if(tmppep[0] != 20)	{
					begtag = 1;
				} else	{
					begtag = 0;
				}
				*ncut = 0;
			} else	{
				*len_lastpep = appendpep(lastpep, *len_lastpep, tmppep, l);
			}
		}
	} else	{
		if(i == length)	{
			*leftover = 0;	// end at the last position of the edge
		} else {
			*leftover = length - i;
		}
	}
	return(midtag);
}

void GraphTrav::continue_extract_peptides(EDGE *inedge, int offset, char *lastpep, int *len_lastpep, int *leftover, int *ncut, char begtag,
	int beg_edge_index, int beg_edge_offset, int end_edge_index, int end_edge_offset)
{
	int	i, j, k, l;
	char	tmpcodon[10], tmppep[3];
	char	*tmplastpep, *tmplastpep2;
	int	tmp_beg_edge_index, tmp_beg_edge_offset;
	int	tmp_beg_edge_index2, tmp_beg_edge_offset2;
	int	tmp_end_edge_index, tmp_end_edge_offset;
	int	tmp_end_edge_index2, tmp_end_edge_offset2;
	int	len_tmplast=0, len_tmplast2=0;
	char	flag_stopcodon;
	char	tmpbegtag, tmpbegtag2;
	int	ncut_cpy, ncut_cpy2;
	EDGELIST	*nextedge, *nextnextedge;

	for(i = kmersize + offset; i <= inedge -> length - 3; i += 3)	{
		end_edge_offset += 3;
		l = trans_codon2aa(&(inedge -> seq[i]), 3, tmppep);
		if(ckcutsite(tmppep[0]) == 1)	(*ncut) ++;
		if(tmppep[0] == 20 || *ncut > max_miscleavage)	{	// stop codon or cleavage site
			if(tmppep[0] != 20)	{
				*len_lastpep = appendpep(lastpep, *len_lastpep, tmppep, l);
				flag_stopcodon = 0;
			} else	{
				flag_stopcodon = 1;
			}
/*
for(j = 0; j < length; j ++)
	printf("%c", na_name[seq[j]]);
printf("\n");
printf("pos5 length %d len_pep %d edges %d %d %d %d\n", inedge -> length, *len_lastpep, beg_edge_index, beg_edge_offset, end_edge_index, end_edge_offset);
getchar();
*/
			num_peptides = copypepseq(lastpep, *len_lastpep, begtag, beg_edge_index, beg_edge_offset, end_edge_index, end_edge_offset, flag_stopcodon);
			*len_lastpep = 0;
			beg_edge_index = end_edge_index;
			beg_edge_offset = end_edge_offset;
			if(tmppep[0] != 20)	{
				begtag = 1;
			} else	{
				begtag = 0;
			}
			*ncut = 0;
			return;
		} else	{
			*len_lastpep = appendpep(lastpep, *len_lastpep, tmppep, l);
		}
	}

	tmplastpep = new char[(*len_lastpep) + 1];
	len_tmplast = 0;	// temporarily store the lastpep sequence
	len_tmplast = appendpep(tmplastpep, len_tmplast, lastpep, *len_lastpep);
	tmp_beg_edge_index = beg_edge_index;
	tmp_beg_edge_offset = beg_edge_offset;
	tmp_end_edge_index = end_edge_index;
	tmp_end_edge_offset = end_edge_offset;
	ncut_cpy = *ncut;
	tmpbegtag = begtag;
	*leftover = inedge -> length - i;
	nextedge = vertex[inedge -> nodeindex[1]].nextedge;
	while(nextedge)	{
/*
printf("continue nextedge %d length %d nodes %d %d multitag %d depth %d\n", nextedge -> edge, nextedge -> edge -> length, nextedge -> edge -> nodeindex[0], nextedge -> edge -> nodeindex[1], nextedge -> edge -> multitag, depth);
printf("leftover %d\n", *leftover);
*/
		if(nextedge -> edge -> multitag == 1)	{
			nextedge = nextedge -> next;
			continue;
		}
		*len_lastpep = 0;	// restore the lastpep sequence
		*len_lastpep = appendpep(lastpep, *len_lastpep, tmplastpep, len_tmplast);
		beg_edge_index = tmp_beg_edge_index;
		beg_edge_offset = tmp_beg_edge_offset;
		end_edge_index = tmp_end_edge_index;
		end_edge_offset = tmp_end_edge_offset;
		*ncut = ncut_cpy;
		begtag = tmpbegtag;
		if(*leftover == 0)	{
			end_edge_index = nextedge -> edge -> index;
			end_edge_offset = kmersize;
			nextedge -> edge -> multitag = 1;
			if(*len_lastpep <= peptide_max && depth <= max_depth)	{
				depth ++;
				continue_extract_peptides(nextedge -> edge, 0, lastpep, len_lastpep, leftover, ncut, begtag,
					beg_edge_index, beg_edge_offset, end_edge_index, end_edge_offset);
				depth --;
			}
		} else	{
			for(i = 0; i < *leftover; i ++)	{
				tmpcodon[i] = inedge -> seq[inedge -> length - *leftover + i];
			}
			if(nextedge -> edge -> length - kmersize <= 0)	{
				printf("edge length error: %d, kmersize %d\n", nextedge -> edge -> length, kmersize);
				exit(-1);
			}
			if(nextedge -> edge -> length - kmersize + *leftover < 3)	{
				tmpcodon[1] = nextedge -> edge -> seq[kmersize];
				tmplastpep2 = new char[(*len_lastpep) + 1];
				len_tmplast2 = 0;	// temporarily store the lastpep sequence
				len_tmplast2 = appendpep(tmplastpep2, len_tmplast2, lastpep, *len_lastpep);
				tmp_beg_edge_index2 = beg_edge_index;
				tmp_beg_edge_offset2 = beg_edge_offset;
				tmp_end_edge_index2 = end_edge_index;
				tmp_end_edge_offset2 = end_edge_offset;
				ncut_cpy2 = *ncut;
				tmpbegtag2 = begtag;
				nextnextedge = vertex[nextedge -> edge ->  nodeindex[1]].nextedge;
				while(nextnextedge)	{
/*
printf("continue nextnextedge %d length %d nodes %d %d multitag %d depth %d\n", nextnextedge -> edge, nextnextedge -> edge -> length, nextnextedge -> edge -> nodeindex[0], nextnextedge -> edge -> nodeindex[1], nextnextedge -> edge -> multitag, depth);
*/
					if(nextnextedge -> edge -> multitag == 1)	{
						nextnextedge = nextnextedge -> next;
						continue;
					}
					*len_lastpep = 0;	// restore the lastpep sequence
					*len_lastpep = appendpep(lastpep, *len_lastpep, tmplastpep2, len_tmplast2);
					end_edge_index = tmp_end_edge_index2;
					end_edge_offset = tmp_end_edge_offset2;
					*ncut = ncut_cpy2;
					tmpcodon[2] = nextnextedge -> edge -> seq[kmersize];
					l = trans_codon2aa(tmpcodon, 3, tmppep);
					end_edge_index = nextnextedge -> edge -> index;
					end_edge_offset = kmersize;
					if(ckcutsite(tmppep[0]) == 1)	(*ncut) ++;
					if(tmppep[0] == 20 || *ncut > max_miscleavage)	{	// stop codon or cleavage site
						if(tmppep[0] != 20)	{
							*len_lastpep = appendpep(lastpep, *len_lastpep, tmppep, l);
							flag_stopcodon = 0;
						} else	{
							flag_stopcodon = 1;
						}
/*
for(j = 0; j < length; j ++)
	printf("%c", na_name[seq[j]]);
printf("\n");
printf("pos6 length %d len_pep %d edges %d %d %d %d\n", nextnextedge -> edge -> length, *len_lastpep, beg_edge_index, beg_edge_offset, end_edge_index, end_edge_offset);
getchar();
*/
						num_peptides = copypepseq(lastpep, *len_lastpep, begtag, beg_edge_index, beg_edge_offset,
								end_edge_index, end_edge_offset, flag_stopcodon);
						*len_lastpep = 0;
						beg_edge_index = end_edge_index;
						beg_edge_offset = end_edge_offset;
						if(tmppep[0] != 20)	{
							begtag = 1;
						} else	{
							begtag = 0;
						}
						*ncut = 0;
					} else	{
						*len_lastpep = appendpep(lastpep, *len_lastpep, tmppep, l);
					}
					nextedge -> edge -> multitag = 1;
					nextnextedge -> edge -> multitag = 1;
					if(*len_lastpep <= peptide_max && depth <= max_depth)	{
						depth ++;
						continue_extract_peptides(nextnextedge -> edge, 1, lastpep, len_lastpep, leftover,
							 ncut, begtag, beg_edge_index, beg_edge_offset, end_edge_index, end_edge_offset);
						depth --;
					}
					nextnextedge -> edge -> multitag = 0;
					nextnextedge = nextnextedge -> next;
				}
				delete[] tmplastpep2;
				begtag = tmpbegtag2;
			} else	{
				for(j = 0; i < 3; i ++, j ++)	{
					tmpcodon[i] = nextedge -> edge -> seq[kmersize + j];
				}
				l = trans_codon2aa(tmpcodon, 3, tmppep);
				end_edge_index = nextedge -> edge -> index;
				end_edge_offset = kmersize + j;
				if(ckcutsite(tmppep[0]) == 1)	(*ncut) ++;
				if(tmppep[0] == 20 || *ncut > max_miscleavage)	{	// stop codon or cleavage site
					if(tmppep[0] != 20)	{
						*len_lastpep = appendpep(lastpep, *len_lastpep, tmppep, l);
						flag_stopcodon = 0;
					} else	{
						flag_stopcodon = 1;
					}
/*
for(j = 0; j < length; j ++)
	printf("%c", na_name[seq[j]]);
printf("\n");
printf("pos7 length %d len_pep %d edges %d %d %d %d\n", nextedge -> edge -> length, *len_lastpep, beg_edge_index, beg_edge_offset, end_edge_index, end_edge_offset);
getchar();
*/
					num_peptides = copypepseq(lastpep, *len_lastpep, begtag, beg_edge_index, beg_edge_offset,
								 end_edge_index, end_edge_offset, flag_stopcodon);
					*len_lastpep = 0;
					beg_edge_index = end_edge_index;
					beg_edge_offset = end_edge_offset;
					if(tmppep[0] != 20)	{
						begtag = 1;
					} else	{
						begtag = 0;
					}
					*ncut = 0;
				} else	{
					*len_lastpep = appendpep(lastpep, *len_lastpep, tmppep, l);
				}
				nextedge -> edge -> multitag = 1;
				if(*len_lastpep <= peptide_max && depth <= max_depth)	{
					depth ++;
					continue_extract_peptides(nextedge -> edge, 3 - *leftover, lastpep, len_lastpep, leftover,
						 ncut, begtag, beg_edge_index, beg_edge_offset, end_edge_index, end_edge_offset);
					depth --;
				}
			}
		}
		nextedge -> edge -> multitag = 0;
		nextedge = nextedge -> next;
	}

	delete[] tmplastpep;
}

int GraphTrav::appendpep(char *lastpep, int len_lastpep, char *tmppep, int length)
{
	int	i;

	for(i = 0; i < length; i ++)	{
		lastpep[len_lastpep ++] = tmppep[i];
	}
	return(len_lastpep);
}


int GraphTrav::trans_codon2aa(char *tmpcodon, int len_codonseq, char *tmppepseq)
{
	int	i, j, k, q;

	q = 0;
	for(i = 0; i < len_codonseq; i += 3)	{
		k = 0;
		for(j = i; j < i + 3; j ++)	{
			k = k * 4 + tmpcodon[j];
		}
		tmppepseq[q ++] = codon2aanum[k];
	}
	return(q);
}

int GraphTrav::copypepseq(char *lastpep, int len_lastpep, char begtag, int beg_edge_index, int beg_edge_offset, int end_edge_index, int end_edge_offset, char flag_stopcodon)
{
	int i;
	char c;
	PEPTIDESEQ *tmppepseq, *parentpepseq;

	if(len_lastpep >= peptide_min && len_lastpep <= peptide_max)	{
		if(begtag == 1 || lastpep[len_lastpep - 1] != 20)	{
			if(!pepseq)	{
				pepseq = new PEPTIDESEQ[1];
				pepseq -> seq = new char[len_lastpep + 1];
				for(i = 0; i < len_lastpep; i ++)	{
					pepseq -> seq[i] = lastpep[i];
				}
				pepseq -> left = pepseq -> right = pepseq -> parent = NULL;
				pepseq -> length = len_lastpep;
/*
if(beg_edge_index == 1876151 && beg_edge_index == 4590727)	{
printf("Input %d %d %d %d\n", beg_edge_index, beg_edge_offset, end_edge_index, end_edge_offset);
getchar();
}
*/
				pepseq -> beg_edge_index = beg_edge_index;
				pepseq -> beg_edge_offset = beg_edge_offset;
				pepseq -> end_edge_index = end_edge_index;
				pepseq -> end_edge_offset = end_edge_offset;
				pepseq -> flag_stopcodon = flag_stopcodon;
				num_peptides ++;
/*
printf("input pepseq %d seq %d left %d right %d\n", pepseq, pepseq -> length, pepseq -> left, pepseq -> right);
*/
				return(num_peptides);
			}
			tmppepseq = pepseq;
			while(tmppepseq)	{
				c = check_pepseq(lastpep, len_lastpep, tmppepseq -> seq, tmppepseq -> length);
				if(c == 2)	{
					parentpepseq = tmppepseq;
					tmppepseq = tmppepseq -> left;
				} else if(c == 1)	{
					parentpepseq = tmppepseq;
					tmppepseq = tmppepseq -> right;
				} else	{
					return(num_peptides);
				}
			}
			if(len_lastpep == 0)	{
				printf("length error\n");
				exit(-1);
			}
			tmppepseq = new PEPTIDESEQ[1];
			tmppepseq -> seq = new char[len_lastpep + 1];
			for(i = 0; i < len_lastpep; i ++)	{
				tmppepseq -> seq[i] = lastpep[i];
			}
			tmppepseq -> length = len_lastpep;
/*
if(beg_edge_index == 1876151 || beg_edge_index == 4590727)	{
printf("tmppepseq %d Length %d Input %d %d %d %d\n", tmppepseq, len_lastpep, beg_edge_index, beg_edge_offset, end_edge_index, end_edge_offset);
getchar();
}
*/
			tmppepseq -> beg_edge_index = beg_edge_index;
			tmppepseq -> beg_edge_offset = beg_edge_offset;
			tmppepseq -> end_edge_index = end_edge_index;
			tmppepseq -> end_edge_offset = end_edge_offset;
			tmppepseq -> flag_stopcodon = flag_stopcodon;
			tmppepseq -> parent = parentpepseq;
			tmppepseq -> left = tmppepseq -> right = NULL;
			num_peptides ++;
			if(c == 2)	{
				parentpepseq -> left = tmppepseq;
			} else	{
				parentpepseq -> right = tmppepseq;
			}
		}
	}
	return(num_peptides);
}

char GraphTrav::check_pepseq(char *seq1, int length1, char *seq2, int length2)
{
	int	i;
	if(length1 < length2)	{
		return(2);
	} else if(length1 > length2)	{
		return(1);
	} else	{
		for(i = 0; i < length1; i ++)	{
			if(seq1[i] < seq2[i])	{
				return(2);
			} else if(seq1[i] > seq2[i])	{
				return(1);
			}
		}
		return(0);
	}
}

void GraphTrav::writefile(char *outputfile)
{
	int	numpep = 0;
	std::ofstream fs(outputfile);
	if(!fs) { cout<<"cannot open file to write "<< outputfile <<endl; exit(0); }
/*
printf("continue %d\n", num_peptides);
*/
	printf("%d peptides generated.\n", num_peptides);
	numpep = outputseq(fs, pepseq, numpep);
/*
printf("continue10 numpep %d\n", numpep);
*/
}

int GraphTrav::outputseq(ofstream &fs, PEPTIDESEQ *pepseq, int numpep)
{
	if(pepseq)	{
		numpep = outputseq(fs, pepseq -> left, numpep);
/*
if(pepseq -> beg_edge_index == 1876151 || pepseq -> beg_edge_index == 4590727)	{
printf("pepseq %d Length %d Input %d %d %d %d\n", pepseq, pepseq -> length, pepseq -> beg_edge_index, pepseq -> beg_edge_offset, pepseq -> end_edge_index, pepseq -> end_edge_offset);
getchar();
}
*/
		output_peptides(fs, numpep, pepseq -> seq, pepseq -> length, pepseq -> beg_edge_index, pepseq -> beg_edge_offset,
				pepseq -> end_edge_index, pepseq -> end_edge_offset, pepseq -> flag_stopcodon);
		free((void *) pepseq -> seq);
		numpep ++;
		numpep = outputseq(fs, pepseq -> right, numpep);
		free((void *) pepseq);
	}
	return(numpep);
}

void GraphTrav::output_peptides(ofstream &fp, int index, char *seq, int length, int beg_edge_index, int beg_edge_offset, int end_edge_index, int end_edge_offset, char flag_stopcodon)
{
	int	i, j, k, l, m, n;
	char	str[1000];	

	sprintf(str, ">peptide%d %d %d %d %d %d %d", index, length, beg_edge_index, beg_edge_offset, end_edge_index, end_edge_offset, flag_stopcodon);
	fp<<str<<endl;
	for(i = 0; i < length; i ++)	{
		sprintf(str, "%c", aalist[seq[i]]);
		fp<<str;
		if(i % 60 == 59)	{
			fp<<endl;
		}
	}
	if(i % 60 != 0)	{
		fp<<endl;	
	}
}
