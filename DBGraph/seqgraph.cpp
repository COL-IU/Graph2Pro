//Package: DBGraph2Pro
//Program Name: seqgraph (the main program for loading de Bruijn graph) 
//Latest modification on April 18th, 2016
//Yuzhen Ye <yye@indiana.edu> and Haixu Tang <hatang@indiana.edu>
//Affiliation: School of Informatics and Computing, Indiana University, Bloomington

#include "seqgraph.h"
#include "soap.h"
#include "FastG.h"

int ifprint = 0;

seqgraph::seqgraph(void) {
	init();
}

//initialize
void seqgraph::init() {
        kmersize = 31;
        set_maskv();
}

void seqgraph::set_maskv()
{
        maskv = 1;
        maskv = maskv << (2 * kmersize);
        maskv = maskv - 1;
}

void seqgraph::cleanup(void)
{
        delete[] vindex;
        
        int i;
        
        for(i = 0; i < num_edge; i ++)  {
                delete[] edge[i].seq;
        }
        delete[] edge; 
        for(i = 0; i < num_vertex; i ++)        {
                soap::delete_edgelist(vertex[i].lastedge);
                soap::delete_edgelist(vertex[i].nextedge);
        }
        delete[] vertex;
}

//load the graph (from SOAPdenovo)
void seqgraph::loadsoap(bool SOAP2_tag, char *edgefile, char *edgeseqfile) 
{
	int i;
	//get num of edges
	cout<<"\n>>>Now read graph from file "<<edgefile<<endl;
	num_edge = soap::countedge(edgefile);
	cout<<"num edge to be allocated: "<<num_edge<<endl;

	//allocate vertex
	vertex = new VERTEX[2 * (num_edge + 1)];
	//initialize the variables -- it is important!
	cout<<"num_vertex to be allocated: "<<2 * (num_edge + 1)<<endl;
	for(i = 0; i < 2 * (num_edge + 1); i ++) {
		vertex[i].lastedge = vertex[i].nextedge = NULL;
		vertex[i].indegree = vertex[i].outdegree = 0;
	}
	num_vertex = soap::readvertex(SOAP2_tag, edgefile, vertex);
	cout<<"num_vertex actually extracted:  "<<num_vertex<<endl;

	//Initialize class attributes (global) -- edges for contigs
	edge = new EDGE[num_edge]; 
	for(i = 0; i < num_edge; i ++) { 
		edge[i].tag = edge[i].multitag = edge[i].length = 0;
		edge[i].bal_edge = 0;
		edge[i].seq = NULL;
		edge[i].nodeindex[0] = edge[i].nodeindex[1] = 0;
	}

	cout<<"SOAP2_tag "<<SOAP2_tag<<" kmersize "<<kmersize<<endl;
	num_edge = soap::readedge(SOAP2_tag, edgefile, edgeseqfile, edge, vertex, num_vertex, kmersize);
	cout<<"input edge extracted: "<<num_edge<<endl;

	for(i = 0; i < num_edge; i ++)	{
		edge[i].index = i;
	}
	cout<<"print_vertex_degree.."<<endl;
	print_vertex_degree();

	cout<<"hashtablesize "<<hashtablesize<<endl;
	//vindex = (int *) smallapp::ckalloc((hashtablesize + 1) * sizeof(int));
	vindex = new int[hashtablesize + 1];
	for(i = 0; i < hashtablesize + 1; i ++) vindex[i] = 0;
	vindex[hashtablesize] = num_vertex;

	cout<<"index vertex..."<<endl;
	index_vertex();

	for(i = hashtablesize - 1; i >= 0; i --)	{
		if(vindex[i] == 0)	{
			vindex[i] = vindex[i + 1];
		}
	}
}

//load the graph (from FastG format; see: http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf)
void seqgraph::loadFastG(char *edgeseqfile)
{
        int i;
        int *edgeindex;

        //get num of edges
        cout<<"\n>>>Now read edges from file "<<edgeseqfile<<endl;
        num_edge = FastG::countedge(edgeseqfile);
        cout<<"num edge to be allocated: "<<num_edge<<endl;

        //allocate vertex
        vertex = new VERTEX[2 * (num_edge + 1)];
        //initialize the variables -- it is important!
        for(i = 0; i < 2 * (num_edge + 1); i ++) {
                vertex[i].lastedge = vertex[i].nextedge = NULL;
                vertex[i].indegree = vertex[i].outdegree = 0;
        }

        edge = new EDGE[num_edge];
        for(i = 0; i < num_edge; i ++) {
                edge[i].tag = edge[i].multitag = edge[i].length = 0;
                edge[i].bal_edge = 0;
                edge[i].seq = NULL;
                edge[i].nodeindex[0] = edge[i].nodeindex[1] = -1;
        }

        edgeindex = new int[num_edge * 2];
        for(i = 0; i < num_edge * 2; i ++) {
                edgeindex[i] = 0;
        }

        cout<<" kmersize "<<kmersize<<endl;
        num_vertex = FastG::readedge(edgeseqfile, edge, vertex, edgeindex, num_edge, kmersize);
        cout<<"input edge extracted: "<<num_edge<<endl;
        cout<<"vertex created: "<<num_vertex<<endl;

        delete[] edgeindex;

        for(i = 0; i < num_edge; i ++)  {
                edge[i].index = i;
        }
        cout<<"print_vertex_degree.."<<endl;
        print_vertex_degree();

        cout<<"hashtablesize "<<hashtablesize<<endl;
        //vindex = (int *) smallapp::ckalloc((hashtablesize + 1) * sizeof(int));
        vindex = new int[hashtablesize + 1];
        for(i = 0; i < hashtablesize + 1; i ++) vindex[i] = 0;
        vindex[hashtablesize] = num_vertex;

        cout<<"index vertex..."<<endl;
        lindex_vertex();

        for(i = hashtablesize - 1; i >= 0; i --)        {
                if(vindex[i] == 0)      {
                        vindex[i] = vindex[i + 1];
                }
        }
}
                                                                   
void seqgraph::index_vertex(void)
{
	int shift = (kmersize - hashw) * 2;
	printf("kmer %d hashw %d shift %d\n", kmersize, hashw, shift);

	int hashv = vertex[0].index >> shift;
	vindex[hashv] = 1;
	//cout<<"index "<<vertex[0].index<<" hashv "<<hashv<<endl;
	printf("index %ld hashv %d\n", vertex[0].index, hashv);
	for(int i = 1; i < num_vertex; i ++)	{
		hashv = vertex[i].index >> shift;
		if(vindex[hashv] == 0)	{
			vindex[hashv] = i + 1;
		}
	}
}

void seqgraph::lindex_vertex(void)
{
	int hashv;
	int shift = (kmersize - hashw) * 2;
	printf("kmer %d hashw %d shift %d\n", kmersize, hashw, shift);

	if(shift >= 128)	{
		hashv = vertex[0].lindex.bits[1] >> (shift - 128);
	} else	{
		hashv = vertex[0].lindex.bits[0] >> shift;
	}
	vindex[hashv] = 1;
	//cout<<"index "<<vertex[0].index<<" hashv "<<hashv<<endl;
	printf("index %ld %ld hashv %d\n", vertex[0].lindex.bits[0], vertex[0].lindex.bits[1], hashv);
	for(int i = 1; i < num_vertex; i ++)	{
		if(shift >= 128)	{
			int hashv = vertex[0].lindex.bits[1] >> (shift - 128);
		} else	{
			int hashv = vertex[0].lindex.bits[0] >> shift;
		}
		if(vindex[hashv] == 0)	{
			vindex[hashv] = i + 1;
		}
	}
}

void seqgraph::print_vertex_degree(void)
{
	int	i, j, n1, n2;
	int	vertex_in_out[5][5];

	for(i = 0; i < 5; i ++)	{
		for(j = 0; j < 5; j ++)	{
			vertex_in_out[i][j] = 0;
		}
	}
	for(i = 0; i < num_vertex; i ++)	{
		if(vertex[i].indegree > 4)	{
			n1 = 4;
		} else	{
			n1 = vertex[i].indegree;
		}
		if(vertex[i].outdegree > 4)	{
			n2 = 4;
		} else	{
			n2 = vertex[i].outdegree;
		}
		vertex_in_out[n1][n2] ++;
	}
	printf("Degree distribution..\n");
	printf("-----------------------\n");
	for(i = 0; i < 5; i ++)	{
		for(j = 0; j < 5; j ++)	{
			printf("%d ", vertex_in_out[i][j]);
		}
		printf("\n");
	}
	printf("-----------------------\n");
}
