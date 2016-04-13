#include "GraphPepTrav.h"
#include "smallapp.h"
#include "soap.h"
#include "seq.h"

void GraphPepTrav::loadpepmap(char *pepseqfile)
{
	int	i, j, k, l;

	input_pepnum(pepseqfile);
	id_peptides = (ID_PEPTIDE *) smallapp::ckalloc(num_id_peptides * sizeof(ID_PEPTIDE));
	input_pepmap(pepseqfile);

printf("generate peptide map ...\n");
	pep2edgemap();

printf("sort peptide map ...\n");
	sort_edgemap();
}

void GraphPepTrav::input_pepnum(char *pepseqfile)
{
	std::ifstream fp_pepseq(pepseqfile);
	if(!fp_pepseq)	{
		cout<<"cannot open file "<<pepseqfile<<endl;
		exit(-1);
	}
	num_id_peptides = 0;
	std::string s;
	while(!fp_pepseq.eof()) {
		std::getline(fp_pepseq, s);
		if(s.empty()) continue;
		num_id_peptides ++;
	}
	fp_pepseq.close();
}

void GraphPepTrav::input_pepmap(char *pepseqfile)
{
	int	i, j, k, l;
	int	count;
	int	num_multi_edge_pep;
	char	temp[10000];
	double	perc;

	std::istringstream sline;

	std::ifstream fp_pepseq(pepseqfile);
	if(!fp_pepseq)	{
		cout<<"cannot open file "<<pepseqfile<<endl;
		exit(-1);
	}
	num_multi_edge_pep = i = 0;
	std::string s;
	string field = "";
	if(!fp_pepseq.eof())	{
		std::getline(fp_pepseq, s);	// skip the first line
	}
	while(!fp_pepseq.eof()) {
		std::getline(fp_pepseq, s);
		if(s.empty()) continue;
		if(s.at(0) == '#')	continue;
		sline.str(s);	// put a line in a line stream
		count = 0;
		while (std::getline(sline, field, '\t' ))	{
//		sscanf(s.c_str(), "%*s%s%d", id_peptides[i].seq, &(id_peptides[i].beg_edge_index), &(id_peptides[i].beg_edge_offset),
//				&(id_peptides[i].end_edge_index), &(id_peptides[i].end_edge_offset), &(id_peptides[i].flag_stopcodon));
//		id_peptides[i].length = strlen(id_peptides[i].seq);
//printf("count %d field %s\n", count, field.c_str());
			if(count == 0)	{
				id_peptides[i].length = field.length();
				id_peptides[i].seq = new char[id_peptides[i].length + 2];
				std:strcpy(id_peptides[i].seq, field.c_str());
				id_peptides[i].length = parse_pep(id_peptides[i].seq, id_peptides[i].length);
			} else if(count == 5)	{
				id_peptides[i].beg_edge_index = std::atoi(field.c_str());
			} else if(count == 6)	{
				id_peptides[i].beg_edge_offset = std::atoi(field.c_str());
			} else if(count == 7)	{
				id_peptides[i].end_edge_index = std::atoi(field.c_str());
			} else if(count == 8)	{
				id_peptides[i].end_edge_offset = std::atoi(field.c_str());
			} else if(count == 9)	{
				id_peptides[i].flag_stopcodon = std::atoi(field.c_str());
			}
/*
			if(!std::strcmp(field.c_str(), "DBGraph"))	count += 2;	// shift columns
			if(count == 9)	{
				id_peptides[i].length = field.length();
				id_peptides[i].seq = new char[id_peptides[i].length + 2];
				std:strcpy(id_peptides[i].seq, field.c_str());
				id_peptides[i].length = parse_pep(id_peptides[i].seq, id_peptides[i].length);
			} else if(count == 18)	{
				id_peptides[i].beg_edge_index = std::atoi(field.c_str());
			} else if(count == 19)	{
				id_peptides[i].beg_edge_offset = std::atoi(field.c_str());
			} else if(count == 20)	{
				id_peptides[i].end_edge_index = std::atoi(field.c_str());
			} else if(count == 21)	{
				id_peptides[i].end_edge_offset = std::atoi(field.c_str());
			} else if(count == 22)	{
				id_peptides[i].flag_stopcodon = std::atoi(field.c_str());
			}
*/
			count ++;
		}
/*
if(id_peptides[i].beg_edge_index == 0 || id_peptides[i].end_edge_index == 0)	{
printf("%s \n", s.c_str());
printf("edge %d seq %s location %d %d %d %d %d\n", i, id_peptides[i].seq, id_peptides[i].beg_edge_index, id_peptides[i].beg_edge_offset, id_peptides[i].end_edge_index, id_peptides[i].end_edge_offset, id_peptides[i].flag_stopcodon);
getchar();
}
*/
		if(id_peptides[i].beg_edge_index != id_peptides[i].end_edge_index)	{
			num_multi_edge_pep ++;
		}
//		i ++;
		if(id_peptides[i].beg_edge_index != id_peptides[i].end_edge_index || (id_peptides[i].end_edge_offset - id_peptides[i].beg_edge_offset) / 3 == 
			id_peptides[i].length)	{
			i ++;
		} else	{
			printf("peptide ");
			for(j = 0; j < id_peptides[i].length; j ++)	{
				printf("%c", aalist[id_peptides[i].seq[j]]);
			}
			printf(" positions %d %d %d %d\n", id_peptides[i].beg_edge_index, id_peptides[i].beg_edge_offset,
				id_peptides[i].end_edge_index, id_peptides[i].end_edge_offset);
		}
		sline.clear();
	}
	num_id_peptides = i;
	perc = ((double) num_multi_edge_pep) * 100 / num_id_peptides;
	printf("%d out of %d (%.2f%%) are peptides spanning multiple edges.\n", num_multi_edge_pep, num_id_peptides, perc);

	fp_pepseq.close();
}

void GraphPepTrav::pep2edgemap(void )
{
	int	i, j, k, l;
	int	frame;
	char	tmppep[10];

	pep_edgemap = new EDGEMAP[num_edge];

/*
i = 1423660;
printf("edge %d length %d\n", i, edge[i].length);
for(j = 158; j < edge[i].length - 2; j += 3)	{
	l = trans_codon2aa(&(edge[i].seq[j]), 3, tmppep);
	printf("%c", aalist[tmppep[0]]);
}
printf("\n");
i = 1563520;
printf("edge %d length %d\n", i, edge[i].length);
for(j = 32; j < edge[i].length - 2; j += 3)	{
	l = trans_codon2aa(&(edge[i].seq[j]), 3, tmppep);
	printf("%c", aalist[tmppep[0]]);
}
printf("\n");
getchar();
*/
 
	for(i = 0; i < num_edge; i ++)	{
		pep_edgemap[i].num_id_pep = new int[3];
		pep_edgemap[i].id_pep_index = (int **) smallapp::ckalloc(3 * sizeof(int *));
		pep_edgemap[i].id_pep_pos = (int **) smallapp::ckalloc(3 * sizeof(int *));
		for(j = 0; j < 3; j ++)	{
			pep_edgemap[i].num_id_pep[j] = 0;
		}
	}

	for(i = 0; i < num_id_peptides; i ++)	{
		frame = id_peptides[i].beg_edge_offset % 3;
/*
if(i == 9581 || i == 8843)	{
if(id_peptides[i].beg_edge_index == 1574368)	{
printf("i %d index %d frame %d num %d\n", i, id_peptides[i].beg_edge_index, frame, pep_edgemap[id_peptides[i].beg_edge_index].num_id_pep[frame]);
if((k = pep_edgemap[id_peptides[i].beg_edge_index].num_id_pep[frame]) > 0)	{
	l = pep_edgemap[id_peptides[i].beg_edge_index].id_pep_index[frame][0];
	printf("id_pep %d index %d %d %d %d\n", l, id_peptides[l].beg_edge_index, id_peptides[l].beg_edge_offset,
		id_peptides[l].end_edge_index, id_peptides[l].end_edge_offset);
}
}
*/
		if(pep_edgemap[id_peptides[i].beg_edge_index].num_id_pep[frame] % max_pep == 0)	{
			if(pep_edgemap[id_peptides[i].beg_edge_index].num_id_pep[frame] > 0)	{
				std::free((void *) pep_edgemap[id_peptides[i].beg_edge_index].id_pep_index[frame]);
				std::free((void *) pep_edgemap[id_peptides[i].beg_edge_index].id_pep_pos[frame]);
			}
			pep_edgemap[id_peptides[i].beg_edge_index].id_pep_index[frame] =
			 (int *) smallapp::ckalloc((pep_edgemap[id_peptides[i].beg_edge_index].num_id_pep[frame] + max_pep) * sizeof(int));
			pep_edgemap[id_peptides[i].beg_edge_index].id_pep_pos[frame] =
			 (int *) smallapp::ckalloc((pep_edgemap[id_peptides[i].beg_edge_index].num_id_pep[frame] + max_pep) * sizeof(int));
		}
		pep_edgemap[id_peptides[i].beg_edge_index].id_pep_index[frame][pep_edgemap[id_peptides[i].beg_edge_index].num_id_pep[frame]] = i;
		pep_edgemap[id_peptides[i].beg_edge_index].id_pep_pos[frame][pep_edgemap[id_peptides[i].beg_edge_index].num_id_pep[frame]] =
				 id_peptides[i].beg_edge_offset;
/*
if(i == 829 || i == 9691)	{
printf("i %d index %d frame %d num %d offset %d end index %d offset %d\n", i, id_peptides[i].beg_edge_index, frame, pep_edgemap[id_peptides[i].beg_edge_index].num_id_pep[frame],
	id_peptides[i].beg_edge_offset, id_peptides[i].end_edge_index, id_peptides[i].end_edge_offset);
getchar();
}
*/
		pep_edgemap[id_peptides[i].beg_edge_index].num_id_pep[frame] ++;
		if(id_peptides[i].end_edge_index != id_peptides[i].beg_edge_index)	{
			frame = id_peptides[i].end_edge_offset % 3;
/*
printf("end i %d index %d num %d\n", i, id_peptides[i].beg_edge_index, pep_edgemap[id_peptides[i].beg_edge_index].num_id_pep[frame]);
*/
			if(pep_edgemap[id_peptides[i].end_edge_index].num_id_pep[frame] % max_pep == 0)	{
				if(pep_edgemap[id_peptides[i].end_edge_index].num_id_pep[frame] > 0)	{
					std::free((void *) pep_edgemap[id_peptides[i].end_edge_index].id_pep_index[frame]);
					std::free((void *) pep_edgemap[id_peptides[i].end_edge_index].id_pep_pos[frame]);
				}
				pep_edgemap[id_peptides[i].end_edge_index].id_pep_index[frame] =
				 (int *) smallapp::ckalloc((pep_edgemap[id_peptides[i].end_edge_index].num_id_pep[frame] + max_pep) * sizeof(int));
				pep_edgemap[id_peptides[i].end_edge_index].id_pep_pos[frame] =
				 (int *) smallapp::ckalloc((pep_edgemap[id_peptides[i].end_edge_index].num_id_pep[frame] + max_pep) * sizeof(int));
			}
			pep_edgemap[id_peptides[i].end_edge_index].id_pep_index[frame][pep_edgemap[id_peptides[i].end_edge_index].num_id_pep[frame]] = i;
			pep_edgemap[id_peptides[i].end_edge_index].id_pep_pos[frame][pep_edgemap[id_peptides[i].end_edge_index].num_id_pep[frame]] =
				id_peptides[i].end_edge_offset;
/*
if(i == 6102 || id_peptides[i].end_edge_index == 46841 && id_peptides[i].end_edge_offset == 37)	{
	printf("index %d offset %d frame %d\n", id_peptides[i].end_edge_index, id_peptides[i].end_edge_offset, frame);
	getchar();
}
*/
			pep_edgemap[id_peptides[i].end_edge_index].num_id_pep[frame] ++;
		}
/*
if(i == 9581 || i == 8843)	{
printf("i %d index %d %d num %d\n", i, id_peptides[i].beg_edge_index, id_peptides[i].end_edge_index, pep_edgemap[id_peptides[i].beg_edge_index].num_id_pep[frame]);
getchar();
}
*/
	}
}

void GraphPepTrav::sort_edgemap(void )
{
	int	i, j, k, l, n;

	for(i = 0; i < num_edge; i ++)	{
		for(j = 0; j < 3; j ++)	{
/*
printf("edge %d length %d frame %d num_id_pep %d\n", i, edge[i].length, j, pep_edgemap[i].num_id_pep[j]);
getchar();
*/
			for(k = 0; k < pep_edgemap[i].num_id_pep[j] - 1; k ++)	{
				for(l = k + 1; l < pep_edgemap[i].num_id_pep[j]; l ++)	{
					if(pep_edgemap[i].id_pep_pos[j][k] > pep_edgemap[i].id_pep_pos[j][l] ||
					   pep_edgemap[i].id_pep_pos[j][k] == pep_edgemap[i].id_pep_pos[j][l] && 
					   id_peptides[pep_edgemap[i].id_pep_index[j][k]].beg_edge_offset == pep_edgemap[i].id_pep_pos[j][k] &&
					   id_peptides[pep_edgemap[i].id_pep_index[j][l]].end_edge_offset == pep_edgemap[i].id_pep_pos[j][l])	{
						n = pep_edgemap[i].id_pep_pos[j][l];
						pep_edgemap[i].id_pep_pos[j][l] = pep_edgemap[i].id_pep_pos[j][k];
						pep_edgemap[i].id_pep_pos[j][k] = n;
						n = pep_edgemap[i].id_pep_index[j][l];
						pep_edgemap[i].id_pep_index[j][l] = pep_edgemap[i].id_pep_index[j][k];
						pep_edgemap[i].id_pep_index[j][k] = n;
					}
				}
			}
			for(k = 0; k < pep_edgemap[i].num_id_pep[j] - 1; k ++)	{
				if(id_peptides[pep_edgemap[i].id_pep_index[j][k]].beg_edge_index != id_peptides[pep_edgemap[i].id_pep_index[j][k]].end_edge_index)
					continue;
				for(l = k + 1; l < pep_edgemap[i].num_id_pep[j]; l ++)	{
					if(id_peptides[pep_edgemap[i].id_pep_index[j][l]].beg_edge_index !=
						 id_peptides[pep_edgemap[i].id_pep_index[j][l]].end_edge_index)	
							continue;
					if(id_peptides[pep_edgemap[i].id_pep_index[j][k]].beg_edge_offset <=
						 id_peptides[pep_edgemap[i].id_pep_index[j][l]].beg_edge_offset &&
					   id_peptides[pep_edgemap[i].id_pep_index[j][k]].end_edge_offset >=
						 id_peptides[pep_edgemap[i].id_pep_index[j][l]].end_edge_offset)	{
							id_peptides[pep_edgemap[i].id_pep_index[j][l]].label = 1;	// peptide l is fully contained in peptide k
					}
					if(id_peptides[pep_edgemap[i].id_pep_index[j][k]].beg_edge_offset >=
						 id_peptides[pep_edgemap[i].id_pep_index[j][l]].beg_edge_offset &&
					   id_peptides[pep_edgemap[i].id_pep_index[j][k]].end_edge_offset <=
						 id_peptides[pep_edgemap[i].id_pep_index[j][l]].end_edge_offset)	{
							id_peptides[pep_edgemap[i].id_pep_index[j][k]].label = 1;	// peptide k is fully contained in peptide l
					}
				}
			}
		}
	}
}

void GraphPepTrav::ConnectMapppedPeptides(void )
{
	int	i, j, k, l, m, n, max_l, n1, k1;
	char	*tmpseq;
	int	len_tmpseq;
	char	tmppep[10];
	int	start_mark = -1;
	ID_PEPTIDE *id_peptide_pointer, *next_id_peptide_pointer;

	tmpseq = new char[max_length];
	len_tmpseq = 0;

	for(i = 0; i < num_edge; i ++)	{
		for(j = 0; j < 3; j ++)		{		// for each frame 1, 2 and 3
/*
if(i == 1574368)	{
printf("edge %d length %d frame %d num_id_pep %d\n", i, edge[i].length, j, pep_edgemap[i].num_id_pep[j]);
}
*/
			for(k = 0; k < pep_edgemap[i].num_id_pep[j] - 1; k ++)	{
				id_peptide_pointer = &(id_peptides[pep_edgemap[i].id_pep_index[j][k]]);
				if(i != id_peptide_pointer -> beg_edge_index && pep_edgemap[i].id_pep_pos[j][k] == id_peptide_pointer -> beg_edge_offset)	{
					continue;	// the segment starts from a different edge
				}
				if(i != id_peptide_pointer -> end_edge_index)	{
					continue;	// the segment starts from a different edge
				}
				if(id_peptide_pointer -> label == 1)	{	//skip identified peptides fully contained in other identified peptides
					continue;
				}
				for(k1 = k + 1; k1 < pep_edgemap[i].num_id_pep[j]; k1 ++)	{
					if(id_peptides[pep_edgemap[i].id_pep_index[j][k1]].label == 0 &&
					   id_peptides[pep_edgemap[i].id_pep_index[j][k1]].beg_edge_index == i)	{
						break;
					}
				}
				if(k1 == pep_edgemap[i].num_id_pep[j])	{
					continue;
				}
				next_id_peptide_pointer = &(id_peptides[pep_edgemap[i].id_pep_index[j][k1]]);
				if(id_peptide_pointer -> flag_stopcodon == 1)	{	//the identified peptide ends with a stop codone
					id_peptide_pointer -> len_next_seq = 0;
					continue;
				}
/*
if(i == 1574368)	{
printf("id_peptide_pointer %d index %d %d %d %d, next_id_peptide_pointer %d index %d %d %d %d\n",
	pep_edgemap[i].id_pep_index[j][k], id_peptide_pointer -> beg_edge_index, id_peptide_pointer -> beg_edge_offset,
	id_peptide_pointer -> end_edge_index, id_peptide_pointer -> end_edge_offset,
	pep_edgemap[i].id_pep_index[j][k1], next_id_peptide_pointer -> beg_edge_index, next_id_peptide_pointer -> beg_edge_offset,
	next_id_peptide_pointer -> end_edge_index, next_id_peptide_pointer -> end_edge_offset);
printf("peptide1 ");
for(m = id_peptide_pointer -> beg_edge_offset; m < id_peptide_pointer -> end_edge_offset; m += 3)	{
	l = trans_codon2aa(&(edge[i].seq[m]), 3, tmppep);
	printf("%c", aalist[tmppep[0]]);
}
printf("\n");
printf("peptide2 ");
for(m = next_id_peptide_pointer -> beg_edge_offset; m < next_id_peptide_pointer -> end_edge_offset; m += 3)	{
	l = trans_codon2aa(&(edge[i].seq[m]), 3, tmppep);
	printf("%c", aalist[tmppep[0]]);
}
printf("\n");
printf("interal peptide ");
for(m = id_peptide_pointer -> end_edge_offset; m < next_id_peptide_pointer -> beg_edge_offset; m += 3)	{
	l = trans_codon2aa(&(edge[i].seq[m]), 3, tmppep);
	printf("%c", aalist[tmppep[0]]);
}
printf("\n");
}
*/
				if(next_id_peptide_pointer -> beg_edge_offset >= id_peptide_pointer -> end_edge_offset)	{
					max_l = min(edge[i].length - 2, next_id_peptide_pointer -> beg_edge_offset);
					for(m = id_peptide_pointer -> end_edge_offset; m < max_l; m += 3)	{
						l = trans_codon2aa(&(edge[i].seq[m]), 3, tmppep);
						if(tmppep[0] == 20)	{	// stop codon
							id_peptide_pointer -> next_id = -1;	// stop codon
							id_peptide_pointer -> next_seq = (char *) smallapp::ckalloc((len_tmpseq + 1) * sizeof(char));
							id_peptide_pointer -> len_next_seq = 0;
							id_peptide_pointer -> len_next_seq = appendpep(id_peptide_pointer -> next_seq,
								id_peptide_pointer -> len_next_seq, tmpseq, len_tmpseq);
/*
if(len_tmpseq > max_length - 3)	{
	printf("pos1 len_tmpseq %d\n", len_tmpseq);
}
*/
							len_tmpseq = 0;
							break;
						} else	{
							len_tmpseq = appendpep(tmpseq, len_tmpseq, tmppep, l);
						}
					}
					if(m < edge[i].length && id_peptide_pointer -> next_id >= 0)	{
						if(m != next_id_peptide_pointer -> beg_edge_offset)	{
							printf("peptide map does not match: edge %d frame %d map %d endpos %d next %d length %d m %d index %d %d\n",
								i, j, k, id_peptide_pointer -> end_edge_offset,
								next_id_peptide_pointer -> beg_edge_offset, edge[i].length, m,
								pep_edgemap[i].id_pep_index[j][k], pep_edgemap[i].id_pep_index[j][k + 1]);
							exit(-1);
						}
//printf("len_tmpseq %d\n", len_tmpseq);
/*
if(i == 1574368)	{
printf("tmpseq ");
for(m = 0; m < len_tmpseq; m ++)	{
	printf("%c", aalist[tmpseq[m]]);
}
printf("\n");
}
*/
						id_peptide_pointer -> next_seq = (char *) smallapp::ckalloc((len_tmpseq + 1) * sizeof(char));
						id_peptide_pointer -> len_next_seq = 0;
						id_peptide_pointer -> len_next_seq = appendpep(id_peptide_pointer -> next_seq,
								id_peptide_pointer -> len_next_seq, tmpseq, len_tmpseq);
						next_id_peptide_pointer -> prev_seq = (char *) smallapp::ckalloc((len_tmpseq + 1) * sizeof(char));
						next_id_peptide_pointer -> len_prev_seq = 0;
						next_id_peptide_pointer -> len_prev_seq = appendpep(next_id_peptide_pointer -> prev_seq,
								next_id_peptide_pointer -> len_prev_seq, tmpseq, len_tmpseq);
/*
if(len_tmpseq > max_length - 3)	{
	printf("pos2 len_tmpseq %d\n", len_tmpseq);
}
printf("continue1\n");
*/
						len_tmpseq = 0;
						id_peptide_pointer -> next_id = pep_edgemap[i].id_pep_index[j][k1] + 1; // increase id by 1
						next_id_peptide_pointer -> prev_id = pep_edgemap[i].id_pep_index[j][k] + 1;
					}
				} else	{
					id_peptide_pointer -> len_next_seq = (next_id_peptide_pointer -> beg_edge_offset - id_peptide_pointer -> end_edge_offset) / 3;
					next_id_peptide_pointer -> len_prev_seq = id_peptide_pointer -> len_next_seq;
					id_peptide_pointer -> next_id = pep_edgemap[i].id_pep_index[j][k1] + 1; // increase id by 1
					next_id_peptide_pointer -> prev_id = pep_edgemap[i].id_pep_index[j][k] + 1;
				}
			}
		}
	}
printf("forward connection done\n");

	len_tmpseq = 0;
	for(i = 0; i < num_edge; i ++)	{
		for(j = 0; j < 3; j ++)		{		// for each frame 1, 2 and 3
			for(k = 1; k < pep_edgemap[i].num_id_pep[j]; k ++)	{ //Note: the first identified peptide in each edge will not be extended backward yet
				id_peptide_pointer = &(id_peptides[pep_edgemap[i].id_pep_index[j][k]]);
				if(id_peptide_pointer -> label == 1)	{
					continue;
				}
				if(id_peptide_pointer -> prev_id == 0)	{	// there is a stop codon between peptides k-1 and k
					for(m = id_peptide_pointer -> beg_edge_offset - 3; m >= 0; m -= 3)	{
						l = trans_codon2aa(&(edge[i].seq[m]), 3, tmppep);
						if(tmppep[0] == 20)	{	// stop codon
							if(start_mark < 0)	{
								start_mark = len_tmpseq;
							}
							n1 = 0;
							id_peptide_pointer -> prev_seq = (char *) smallapp::ckalloc((len_tmpseq + 1) * sizeof(char));
							for(n = len_tmpseq - 1; n >= 0; n --)	{
								id_peptide_pointer -> prev_seq[n1 ++] = tmpseq[n];
							}
							id_peptide_pointer -> len_prev_seq = n1;
							start_mark = -1;
/*
if(len_tmpseq > max_length - 3)	{
	printf("pos3 edge %d len_tmpseq %d\n", i, len_tmpseq);
}
*/
							len_tmpseq = 0;
							break;
						} else if(tmppep[0] == 10)	{ // start codon
							len_tmpseq = appendpep(tmpseq, len_tmpseq, tmppep, l);
							start_mark = len_tmpseq;
							n1 = 0;
							id_peptide_pointer -> prev_seq = (char *) smallapp::ckalloc((len_tmpseq + 1) * sizeof(char));
							for(n = len_tmpseq - 1; n >= 0; n --)	{
								id_peptide_pointer -> prev_seq[n1 ++] = tmpseq[n];
							}
							id_peptide_pointer -> len_prev_seq = n1;
							len_tmpseq = 0;
							break;
						}
						len_tmpseq = appendpep(tmpseq, len_tmpseq, tmppep, l);
					}
					id_peptide_pointer -> prev_id = -1;	// stop codon
				}
			}
		}
	}
printf("backward connection done\n");
	delete[] tmpseq;
}

void GraphPepTrav::GraphPep2Pro(void)
{
	int	i, j, k, l, n, m;
	int	offset;
	char	*tmpseq, *markseq;
	ID_PEPTIDE	*id_peptide_pointer;
	int	len_tmpseq;
	int	tot_len = 0;

	cout<<"\n>>>Now Traverse the graph..."<<endl;
	//i = 610685;
	printf("Check at the beginning of graph traversal\n");
        //printf("Edge i = %d tag %d, begtag %d endtag %d length %d pathseg-length %d\n", i, edge[i].tag, edge[i].begtag, edge[i].endtag, edge[i].length, pathlist[i] -> length);

	TraverseGraph();
	printf("traverse graph done\n");
	printf("Connect peptides into proteins.\n");

	seq_proteins = (char **) smallapp::ckalloc(num_id_peptides * sizeof(char *));
	mark_proteins = (char **) smallapp::ckalloc(num_id_peptides * sizeof(char *));
	len_proteins = (int *) smallapp::ckalloc(num_id_peptides * sizeof(int));
	num_proteins = 0;

	tmpseq = (char *) smallapp::ckalloc(max_length * sizeof(char));
	markseq = (char *) smallapp::ckalloc(max_length * sizeof(char));

	for(i = 0; i < max_length; i ++)	{
		tmpseq[i] = markseq[i] = 0;
	}

	len_tmpseq = 0;
	for(i = 0; i < num_id_peptides; i ++)	{
		if(id_peptides[i].label == 1)	{
			continue;
		}
		id_peptide_pointer = &(id_peptides[i]);
		if(id_peptide_pointer -> prev_id <= 0)	{
			if(id_peptide_pointer -> len_prev_seq > 0 && ckcutsite(id_peptide_pointer -> prev_seq[id_peptides[i].len_prev_seq - 1]) == 1)	{ // cut site
				len_tmpseq = appendpep(tmpseq, len_tmpseq, id_peptide_pointer -> prev_seq, id_peptide_pointer -> len_prev_seq);
			}
			while(id_peptide_pointer -> next_id > 0)	{
				l = len_tmpseq;
				len_tmpseq = appendpep(tmpseq, len_tmpseq, id_peptide_pointer -> seq, id_peptide_pointer -> length);
				markpeptide(markseq, l, len_tmpseq);
				id_peptide_pointer -> label = 1;
				n = id_peptide_pointer -> next_id - 1;
				if(id_peptide_pointer -> len_next_seq > 0)	{
/*
if(i == 0)	{
printf("nextseq ");
for(m = 0; m < id_peptide_pointer -> len_next_seq; m ++)	{
	printf("%c", aalist[id_peptide_pointer -> next_seq[m]]);
}
printf("\n");
}
*/
					len_tmpseq = appendpep(tmpseq, len_tmpseq, id_peptide_pointer -> next_seq, id_peptide_pointer -> len_next_seq);
				} else	{
					len_tmpseq += id_peptide_pointer -> len_next_seq;
				}
				id_peptide_pointer = &(id_peptides[n]);
			}
			l = len_tmpseq;
			len_tmpseq = appendpep(tmpseq, len_tmpseq, id_peptide_pointer -> seq, id_peptide_pointer -> length);
			markpeptide(markseq, l, len_tmpseq);
			if(id_peptides[n].len_next_seq > 0)	{
				len_tmpseq = appendpep(tmpseq, len_tmpseq, id_peptides[n].next_seq, id_peptides[n].len_next_seq);
			}
			id_peptide_pointer -> label = 1;
/*
if(len_tmpseq > 1000)	{
	printf("len_tmpseq %d\n", len_tmpseq);
	getchar();
}
*/
			seq_proteins[num_proteins] = (char *) smallapp::ckalloc((len_tmpseq + 1) * sizeof(char));
			mark_proteins[num_proteins] = (char *) smallapp::ckalloc((len_tmpseq + 1) * sizeof(char));
			len_proteins[num_proteins] = 0;
			len_proteins[num_proteins] = appendpep(seq_proteins[num_proteins], len_proteins[num_proteins], tmpseq, len_tmpseq);
/*
if(i == 0)	{
printf("num_proteins %d length %d proseq ", num_proteins, len_proteins[num_proteins]);
for(m = 0; m < len_tmpseq; m ++)	{
	printf("%c", aalist[tmpseq[m]]);
}
printf("\n");
for(m = 0; m < len_proteins[num_proteins]; m ++)	{
	printf("%c", aalist[seq_proteins[num_proteins][m]]);
}
printf("\n");
}
printf("num_proteins %d\n", num_proteins);
for(j = 0; j < len_proteins[num_proteins]; j ++)	{
	printf("%d ", markseq[j]);
}
printf("\n");
getchar();
*/
			for(j = 0; j < len_proteins[num_proteins]; j ++)	{
/*
if(seq_proteins[num_proteins][j] >= 20)	{
	printf("i %d num_proteins %d pos %d seq %d\n", i, num_proteins, j, seq_proteins[num_proteins][j]);
	getchar();
}
*/
				mark_proteins[num_proteins][j] = markseq[j];
			}
			tot_len += len_proteins[num_proteins];
			num_proteins ++;
//if(num_proteins == 1)	printf("len %d\n", len_proteins[num_proteins - 1]);
			for(j = 0; j < len_tmpseq; j ++)	{
				markseq[j] = 0;
			}
			len_tmpseq = 0;
/*
printf("peptide %d num_protein %d tot_len %d num_id_peps %d\n", i, num_proteins, tot_len, num_id_peptides);
*/
		}
	}

	free((void *) tmpseq);
	free((void *) markseq);
}

void GraphPepTrav::markpeptide(char *markseq, int i1, int i2)
{
	int	i;
	for(i = i1; i < i2; i ++)	{
		markseq[i] = 1;
	}
}

void GraphPepTrav::TraverseGraph(void)
{
	int	i, j, k, l, n;
	EDGE	*start_edge;
	int	offset;
	char	*lastpep;
	int	len_lastpep;

	lastpep = (char *) smallapp::ckalloc(max_length * sizeof(char));

	for(i = 0; i < num_id_peptides; i ++)	{
		if(id_peptides[i].next_id == 0)	{
			len_lastpep = 0;
			start_edge = &edge[id_peptides[i].end_edge_index];
/*
printf("start edge %d number peptides %d %d %d ...\n", id_peptides[i].end_edge_index, pep_edgemap[id_peptides[i].end_edge_index].num_id_pep[0],
		pep_edgemap[id_peptides[i].end_edge_index].num_id_pep[1], pep_edgemap[id_peptides[i].end_edge_index].num_id_pep[2]);
*/

			depth = 0;
			id_peptides[i].next_id = extract_next_peptides(start_edge, id_peptides[i].end_edge_offset, lastpep, &len_lastpep);
/*
printf("edge %d next_id %d\n", i, id_peptides[i].next_id);
*/
			if(id_peptides[i].next_id > 0)	{
				id_peptides[id_peptides[i].next_id - 1].prev_id = i + 1;
				id_peptides[i].len_next_seq = 0;
				id_peptides[i].next_seq = (char *) smallapp::ckalloc((len_lastpep + 1) * sizeof(char));
				id_peptides[i].len_next_seq = appendpep(id_peptides[i].next_seq, id_peptides[i].len_next_seq,
					lastpep, len_lastpep);
				id_peptides[id_peptides[i].next_id - 1].len_prev_seq = 0;
				id_peptides[id_peptides[i].next_id - 1].prev_seq = (char *) smallapp::ckalloc((len_lastpep + 1) * sizeof(char));
				id_peptides[id_peptides[i].next_id - 1].len_prev_seq = appendpep(id_peptides[id_peptides[i].next_id - 1].prev_seq,
					 id_peptides[id_peptides[i].next_id - 1].len_prev_seq, lastpep, len_lastpep);
			} else	{
				len_lastpep = 0;
				id_peptides[i].next_id = -1;
				len_lastpep = extract_next_peptides_same_edge(start_edge, id_peptides[i].end_edge_offset, lastpep, len_lastpep);
				id_peptides[i].len_next_seq = 0;
				id_peptides[i].next_seq = (char *) smallapp::ckalloc((len_lastpep + 1) * sizeof(char));
				id_peptides[i].len_next_seq = appendpep(id_peptides[i].next_seq, id_peptides[i].len_next_seq,
					lastpep, len_lastpep);
			}
		}
	}
	printf("Back-extend peptides...\n");
	for(i = 0; i < num_id_peptides; i ++)	{
		if(id_peptides[i].label == 0 && id_peptides[i].prev_id == 0)	{
			len_lastpep = 0;
			start_edge = &edge[id_peptides[i].beg_edge_index];
			len_lastpep = extract_prev_peptides_same_edge(start_edge, id_peptides[i].beg_edge_offset, lastpep, len_lastpep);
			id_peptides[i].len_prev_seq = 0;
			id_peptides[i].prev_seq = (char *) smallapp::ckalloc((len_lastpep + 1) * sizeof(char));
			id_peptides[i].len_prev_seq = appendpep(id_peptides[i].prev_seq, id_peptides[i].len_prev_seq, lastpep, len_lastpep);
			id_peptides[i].prev_id = -1;
		}
	}

	free((void *) lastpep);
}

int GraphPepTrav::extract_next_peptides(EDGE *inedge, int offset, char *lastseq, int *len_lastseq)
{
	int	i, j, k, l;
	int	next_id = 0;
	int	leftover;
	char	frame;
	int	frameindex;
	int	edgeindex;
	char	*lastseq_cpy, *lastseq_cpy2;
	int	len_lastseq_cpy, len_lastseq_cpy2;
	EDGELIST	*nextedge, *nextnextedge;
	char	tag = 0;
	char	tmppep[10];
	char	tmpcodon[10];

	lastseq_cpy = (char *) smallapp::ckalloc(max_length * sizeof(char));
	lastseq_cpy2 = (char *) smallapp::ckalloc(max_length * sizeof(char));

/*
printf("inedge %d %d len_lastseq %d \n", inedge -> index, inedge -> length, *len_lastseq);
*/
	for(i = offset; i <= inedge -> length - 3; i += 3)	{
		l = trans_codon2aa(&(inedge -> seq[i]), 3, tmppep);
		if(tmppep[0] == 20)	{	// stop codon
			tag = 1;
			next_id = -1;
			break;
		}
		*len_lastseq = appendpep(lastseq, *len_lastseq, tmppep, l);
	}
	if(tag == 0)	{
		leftover = inedge -> length - i;
		frame = (kmersize - leftover) % 3;
/*
printf("edgelength %d offset %d leftover %d frame %d depth %d\n", inedge -> length, offset, leftover, frame, depth);
*/
		nextedge = vertex[inedge -> nodeindex[1]].nextedge;
		len_lastseq_cpy = 0;
		len_lastseq_cpy = appendpep(lastseq_cpy, len_lastseq_cpy, lastseq, *len_lastseq);
		while(nextedge)	{
			*len_lastseq = 0;
			*len_lastseq = appendpep(lastseq, *len_lastseq, lastseq_cpy, len_lastseq_cpy);
			edgeindex = nextedge -> edge -> index;
			if(pep_edgemap[edgeindex].num_id_pep[frame] > 0)	{
				if(id_peptides[pep_edgemap[edgeindex].id_pep_index[frame][0]].beg_edge_index != edgeindex)	{
					nextedge = nextedge -> next;
					continue;
				}
				frameindex = id_peptides[pep_edgemap[edgeindex].id_pep_index[frame][0]].beg_edge_offset;
				if(frameindex % 3 != frame)	{
					printf("Frame does not match: frame %d first frame %d\n", frame, frameindex);
					printf("next %d pep %d %d %d %d\n", edgeindex, id_peptides[pep_edgemap[edgeindex].id_pep_index[frame][0]].beg_edge_index,
						id_peptides[pep_edgemap[edgeindex].id_pep_index[frame][0]].beg_edge_offset,
						id_peptides[pep_edgemap[edgeindex].id_pep_index[frame][0]].end_edge_index,
						id_peptides[pep_edgemap[edgeindex].id_pep_index[frame][0]].end_edge_offset);
					exit(-1);
				}
				for(i = 0; i < frameindex; i += 3)	{
					l = trans_codon2aa(&(nextedge -> edge -> seq[i]), 3, tmppep);
					if(tmppep[0] == 20)	{	// stop codon
						next_id = -1;
						tag = 1;
						break;
					} else	{
						*len_lastseq = appendpep(lastseq, *len_lastseq, tmppep, l);
					}
				}
				if(tag == 0)	{
					next_id = pep_edgemap[edgeindex].id_pep_index[frame][0];
					break;
				}
			} else if(leftover == 0)	{
				if(depth <= max_depth)	{
					depth ++;
					next_id = extract_next_peptides(nextedge -> edge, kmersize, lastseq, len_lastseq);
					depth --;
					if(next_id > 0)	{
						break;
					}
				}
			} else	{
				for(i = 0; i < leftover; i ++)	{
					tmpcodon[i] = inedge -> seq[inedge -> length - leftover + i];
				}
				if(nextedge -> edge -> length - kmersize <= 0)	{
					printf("edge length error: %d, kmersize %d\n", nextedge -> edge -> length, kmersize);
					exit(-1);
				}
/*
printf("nextedge index %d length %d leftover %d\n", nextedge -> edge -> index, nextedge -> edge -> length, leftover);
*/
				if(nextedge -> edge -> length - kmersize + leftover < 3)	{
					tmpcodon[1] = nextedge -> edge -> seq[kmersize];
					nextnextedge = vertex[nextedge -> edge -> nodeindex[1]].nextedge;
					len_lastseq_cpy2 = 0;
					len_lastseq_cpy2 = appendpep(lastseq_cpy2, len_lastseq_cpy2, lastseq, *len_lastseq);
					while(nextnextedge)	{
						*len_lastseq = 0;
						*len_lastseq = appendpep(lastseq, *len_lastseq, lastseq_cpy2, len_lastseq_cpy2);
						tmpcodon[2] = nextnextedge -> edge -> seq[kmersize];
						l = trans_codon2aa(tmpcodon, 3, tmppep);
						if(tmppep[0] == 20)	{	// stop codon or cleavage site
							next_id = -1;
							tag = 1;
							break;
						} else	{
							*len_lastseq = appendpep(lastseq, *len_lastseq, tmppep, l);
						}
						if(depth <= max_depth)	{
							depth ++;
							next_id = extract_next_peptides(nextnextedge -> edge, kmersize + 1, lastseq, len_lastseq);
							depth --;
							if(next_id > 0)	{
								break;
							}
						}
						nextnextedge = nextnextedge -> next;
					}
				} else	{
					for(j = 0; i < 3; i ++, j ++)	{
						tmpcodon[i] = nextedge -> edge -> seq[kmersize + j];
					}
					l = trans_codon2aa(tmpcodon, 3, tmppep);
					if(tmppep[0] == 20)	{	// stop codon or cleavage site
						next_id = -1;
						tag = 1;
						break;
					} else	{
						*len_lastseq = appendpep(lastseq, *len_lastseq, tmppep, l);
					}
					if(depth <= max_depth)	{
						depth ++;
						next_id = extract_next_peptides(nextedge -> edge, kmersize + j, lastseq, len_lastseq);
						depth --;
						if(next_id > 0)	{
							break;
						}
					}
				}
			}
			nextedge = nextedge -> next;
		}
	}

	free((void *) lastseq_cpy);
	free((void *) lastseq_cpy2);
	return(next_id);
}

int GraphPepTrav::extract_next_peptides_same_edge(EDGE *inedge, int offset, char *lastseq, int len_lastseq)
{
	int	i, j, k, l;
	int	next_id = 0;
	int	leftover;
	char	frame;
	EDGELIST	*nextedge, *nextnextedge;
	char	tmppep[10];

	for(i = offset; i <= inedge -> length - 3; i += 3)	{
		l = trans_codon2aa(&(inedge -> seq[i]), 3, tmppep);
		if(tmppep[0] == 20)	{	// stop codon
			break;
		}
		len_lastseq = appendpep(lastseq, len_lastseq, tmppep, l);
	}
	return(len_lastseq);
}

int GraphPepTrav::extract_prev_peptides_same_edge(EDGE *inedge, int offset, char *lastseq, int len_lastseq)
{
	int	i, j, k, l;
	char	tmppep[10], c;

	for(i = offset - 3; i >= 0; i -= 3)	{
		l = trans_codon2aa(&(inedge -> seq[i]), 3, tmppep);
		if(tmppep[0] == 20)	{	// stop codon
			break;
		} else if(tmppep[0] == 10)	{	// start codon
			len_lastseq = appendpep(lastseq, len_lastseq, tmppep, l);
			break;
		}
		len_lastseq = appendpep(lastseq, len_lastseq, tmppep, l);
	}

	for(i = 0; i < len_lastseq / 2; i ++)	{
		c = lastseq[i];
		lastseq[i] = lastseq[len_lastseq - 1 - i];
		lastseq[len_lastseq - 1 - i] = c;
	}
	return(len_lastseq);
}

int GraphPepTrav::appendpep(char *lastpep, int len_lastpep, char *tmppep, int length)
{
	int	i;

	for(i = 0; i < length; i ++)	{
		lastpep[len_lastpep ++] = tmppep[i];
	}
	return(len_lastpep);
}


int GraphPepTrav::trans_codon2aa(char *tmpcodon, int len_codonseq, char *tmppepseq)
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

void GraphPepTrav::writefile(char *outputfile)
{
	int	i, j, num_lp = 0, tot_lpl = 0;
	std::ofstream fs(outputfile);
	if(!fs) { cout<<"cannot open file to write "<< outputfile <<endl; exit(0); }
/*
printf("continue %d\n", num_peptides);
*/
	printf("%d proteins generated.\n", num_proteins);
	for(i = 0; i < num_proteins; i ++)	{
		if(len_proteins[i] >= 150)	{
			tot_lpl += len_proteins[i];
			num_lp ++;
		}
		outputseq(i, fs, seq_proteins[i], mark_proteins[i], len_proteins[i]);
	}
	printf("# long proteins (>150 aa): %d, total length: %d\n", num_lp, tot_lpl);
	fs.close();
/*
printf("continue10 numpep %d\n", numpep);
*/
}

void GraphPepTrav::outputseq(int index, ofstream &fp, char *seq, char *mark, int length)
{
	int	i, j, k, l, m, n;
	char	str[1000];	

	sprintf(str, ">Protein%d %d", index, length);
	fp<<str<<endl;
	for(i = 0; i < length; i ++)	{
		if(mark[i] == 1)	{
			sprintf(str, "%c", aalist[seq[i]]);
		} else	{
			sprintf(str, "%c", aalist[seq[i]] - 'A' + 'a');
		}
		fp<<str;
		if(i % 60 == 59)	{
			fp<<endl;
		}
	}
	if(i % 60 != 0)	{
		fp<<endl;	
	}
}

void GraphPepTrav::outnucseq(int index, ofstream &fp, char *seq, int length)
{
	int	i, j, k, l, m, n;
	char	str[1000];	

	sprintf(str, ">Transcript%d %d", index, length);
	fp<<str<<endl;
	for(i = 0; i < length; i ++)	{
		sprintf(str, "%c", na_name[seq[i]] - 'a' + 'A');
		fp<<str;
		if(i % 60 == 59)	{
			fp<<endl;
		}
	}
	if(i % 60 != 0)	{
		fp<<endl;	
	}
}

char GraphPepTrav::ckcutsite(char seq)
{
	int	i;
	for(i = 0; i < num_cutsite; i ++)	{
		if(seq == cutsite[i])	{
			return(1);
		}
	}
	return(0);
}

int GraphPepTrav::parse_pep(char *seq, int length)
{
	int	i, j, l;

	l = 0;
	for(i = 0; i < length; i ++)	{
		if(seq[i] >= 'A' && seq[i] <= 'Z' || seq[i] >= 'a' && seq[i] <= 'z')	{
			seq[l ++] = seq::aa2int(seq[i]);
		}
	}
	return(l);
}
