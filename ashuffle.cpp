
/* 
 * File	 : ashuffle.cpp
 * Author: Hasan and Maleq
 * Date	 : August 29, 2015
 */


#include "utility.hpp"
#include "Timer.hpp"
#include "Serial_Multinomial.hpp"

#include <map>
#include <iostream>  
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <climits>

#define  MAX_AGE   128	// maximum age of a person
#define  MAX_SIZE  8	// maximum allowed clique size
#define  MAX_V     ((Vertex)0xFFFFFFFF)
#define  NONE      MAX_V
#define  MAX_ATTEMPT 1000

#define  Min(x,y) ((x)<(y) ? (x) : (y))
#define  Max(x,y) ((x)>(y) ? (x) : (y))
//#define  MAX_SLEVEL 100

#ifndef NULL
	#define  NULL 0
#endif

#define  FNAME_LEN 512	//maximum length of graph file name
#define  MAX_LABEL 6

using namespace std;

typedef  long int	LI;
typedef  long int  	Vertex;
typedef  double  	Weight;
typedef  double  	LD;
typedef  int      	Byte;
typedef  Vertex*	Vlist;

struct Pair {
	Vertex u, v;
};

typedef Pair*  Bin; 

struct Edge {
	Vertex  v;		// the other end vertex of an edge
	//Weight	weight;	// weight of an edge
	//Byte	type;	// label of the other end vertex of an edge
	Vertex 	Key() { return v; }
};

typedef  Edge*  Elist;

class Graph {
public:
    LI   	n;			// # of vertices
	LI   	m;			// # of edges
    Elist	*nlist;		// neighbor list. nlist[i] is the neighbor list of vertex i
	LI		**mlist;
	LI   	*deg;		// deg[i] is the degree of vertex i
	LI		*tdeg;
    Byte  	*label;		// labels of vertices	
    char	gfname[FNAME_LEN];	// graph file name: set in ReadGraph

	Graph() {
		nlist = NULL;
		mlist = NULL;
		deg   = NULL; 
		tdeg  = NULL;
		label = NULL;
		n = m = 0;
	}
	
	~Graph() {
		FreeGraph();
	}
	
	void ReadGraph(const char *fname, const char *demfile);
	void Shuffle(double delta, const char *ofile, int wflag);
	void WriteGraph(const char *ofile, Bin *bin, LI *bsize, LI nbin, bool *eqlLabel);
	void FreeGraph();
};



void Graph::ReadGraph(const char *fname, const char *demfile=NULL)
{
	Vertex	u, v;					// vertices
	LI   	i, k, tmpdeg, lu, lv;
	Vertex  vidx;					// vertex index
	Elist	nbr;					// neighbor
	Weight 	wt;
	LI  	type;

	strcpy(gfname, fname);
	ifstream ifp(fname);
	if(!ifp.is_open()) {
		cout << "Cannot open " << fname << endl;
		exit(1);
	}
	ifp >> n;		// read first line -- the number of vertex
	
	// read dem file and create labels for the vertices
	cout << "Reading the demographic file: " << demfile << " ..."; cout.flush();
	label = new Byte[n];
	for (i=0; i<n; i++)
    	label[i] = INT_MAX;
	
	ifstream dfp(demfile);
	int person_age;
	if(!dfp.is_open()) {
        cout << "Cannot open demographic file: " << fname << endl;
        exit(1);	
	}
	//dfp.ignore(1000,'\n');				// ignore first line: instructions
	while(!dfp.eof()) {
		//dfp.ignore(1000,'\t');		// ignore first column
		dfp >> u >> person_age;			// read id and age
		//dfp.ignore(1000,'\n');			// ignore the rest of the line
		if(u >=0 && u < n)
			label[u] = person_age;  	// label using ages only
	}
	dfp.close();
	cout << "done." << endl;
	
	// Find if any age label is missing, if so, print the first 10 of them
	int error_count = 0;
	int error_limit = 10;
	for (i=0; i<n && error_count <= error_limit; i++) {
		if (label[i] == INT_MAX) {
			cout << "Invalid age: Id " << i << ", Age " << label[i] << endl;
			error_count++;
		}
	}

    nlist = new Elist[n];
	mlist = new LI*[n];
	deg   = new LI[n];
	tdeg  = new LI[n];
	vidx  = 0;
	cout << "Reading graph file: " << fname << " ... "; cout.flush();
 	for (i=0; i<n; i++) {
    	ifp >> u >> tmpdeg;					// read nodes and its degree
		lu = label[u];
		tdeg[u] = 0;
		nbr = nlist[u] = new Edge[tmpdeg];	//neighbor list
		mlist[u] = new LI[tmpdeg];
		for(k=0; k<tmpdeg; k++) {
      		ifp >> v;		// read the adjacent nodes
			lv = label[v];
			if(lu > lv || (lu == lv && u > v))
				nbr[tdeg[u]++].v = v;		// add the edge (u, v) to the graph
			//nbr[k].weight = wt;
			//nbr[k].type = (Byte) type;
		}
		deg[u] = tdeg[u];
	}
	ifp.close();
	cout << " done." << endl;
}


void Graph::FreeGraph()
{
	Timer tf(1);
	FreeMem(deg);
	FreeMem(label);
	FreeAll(nlist, n);
	FreeAll(mlist, n);
	cout << "Freeing the allocated memories inside the destructor function took: "; tf.printsec(); cout << endl;
}

void Graph::Shuffle(double delta, const char *ofile, int wflag)
{
	LI  i, j, k, idx, x, y, repeatedAttempt = 0;
	Bin   *bin;				// edge bins
	LI *bsize, *bidx;		// bin size, bin index (beginning index) 

	Timer total(1);
	Timer t(1);
	Timer tlab(1);
	
	// Find the maximum label; runtime: O(n)
	LI lmax = FindMax(label, n);
	bidx 	= new LI[lmax+1];// additional 1 for degree zero
	bidx[0] = 0;
	for (i=1; i<= lmax; i++)// bidx is used to save time for computing index 
		bidx[i] = bidx[i-1] + i;

	LI nbin = bidx[lmax] + lmax + 1;	// number of bins
	bsize = new LI[nbin];
	for (i=0; i<nbin; i++)
		bsize[i] = 0;
	cout << "Finding lmax, creating and initializing \"bidx\" and \"bsize\" array took: "; tlab.printsec(); cout << endl;
	
	Vertex v;
	LI li, lv;
	Elist nbr;
	
	Timer tcount(1);
	// count the size of each bin; runtime: O(m)
	for (i=0; i<n; i++) {
		nbr = nlist[i];
		li  = label[i];
		for (k=0; k<tdeg[i]; k++) {
			v  = nbr[k].v;
			lv = label[v];
			bsize[bidx[li]+lv]++;
			nlist[v][deg[v]].v = i;
			mlist[i][k] = deg[v];
			mlist[v][deg[v]] = k;
			deg[v]++;
		}
	}
	FreeMem(tdeg);
	cout << "Counting the size of each bin and mapping the counterpart of each edge took: "; tcount.printsec(); cout << endl;
	
	Timer tprob(1);
	LD *prob 	= new LD[nbin];	// probability of an edge switch operation at each bin
	LI *nswitch = new LI[nbin];	// #Shuffle performed at each bin
	
	// count the total number of edges
	m = 0;
	for (i=0; i<nbin; i++)
		m += bsize[i];	// each edge is counted once

	// calculate the probability of an edge switch operation at each bin
	for (i=0; i<nbin; i++)
		prob[i] = LD(bsize[i])/m;
	
	// determine the number of samples
	LI s;		// number of samples
	if (delta < 0.9999)
		s = LI( - double(m)/2.0 * log(1.0 - delta) );
	else
		s = LI( double(m)/2.0 * log(double(m)) );
	cout << "Calculating probability values and #edges took: "; tprob.printsec(); cout << endl;
	
	cout << "|V| = " << n << ", |E| = " << m << ", fraction = " << delta << ", #shuffle = " << s << endl;
	cout << "Max label = " << lmax << ", No. of bins = " << nbin << endl << endl;
	
	// determine which bin performs how many edge switch operations by generating multinomial random variables
	Timer tmult(1);
	Multinomial(s, prob, nbin, nswitch);
	cout << "Computing multinomial distribution took: "; tmult.printsec(); cout << endl;
	FreeMem(prob);	// prob is no longer needed ... so free it
	
	/*
	// checking whether multinomial random variables were generated properly
	LI dummys = 0;
	for (i=0; i<nbin; i++)
		dummys += nswitch[i];
	if(s != dummys) {
		cout << "Error in multinomial random variable generation ..." << endl;
		exit(1);
	}
	*/
	
	// create the bins
	Timer tbincreate(1);
	bin = new Bin[nbin];
	bool *eqlLabel = new bool[nbin];
	LI   **mapEdge = new LI *[nbin];
	cout << "\nPreparing the bins ..." << endl << endl; cout.flush();
	for(i=k=0; i<=lmax; i++) {
		for(j=0; j<i; j++) {
			bin[k] = new Pair[bsize[k]];	// will keep 1 copy for each edge for equal labels: lu != lv
			eqlLabel[k] = false;
			mapEdge[k] = NULL;
			bsize[k++] = 0; // reinitialize the size to count the actual number of pairs in the bin
		}
		
		bin[k] = new Pair[2*bsize[k]];	// will keep 2 copy for each edge for equal labels: lu  = lv
		mapEdge[k] = new LI[2*bsize[k]];
		eqlLabel[k] = true;
		bsize[k++] = 0; // reinitialize the size to count the actual number of pairs in the bin
	}
	
	cout << "Allocating/creating the bins took: "; tbincreate.printsec(); cout << endl;
	LI bi;
	// Now populate the bins with edges
	Timer tbinpop(1);
	for (i=0; i<n; i++) {
		nbr = nlist[i];
		li  = label[i];
		for (k=0; k<deg[i]; k++) {
			v  = nbr[k].v;
			lv = label[v];
			if (li > lv) {		// insert each edge only once for li > lv and twice for li = lv
				bi = bidx[li] + lv;
				bin[bi][bsize[bi]].u = i;
				bin[bi][bsize[bi]].v = v;
				bsize[bi]++;
			}
			else if (li==lv) {
				bi = bidx[li] + lv;
				bin[bi][bsize[bi]].u = i;
				bin[bi][bsize[bi]].v = v;
				mapEdge[bi][bsize[bi]] = mlist[i][k];
				mlist[i][k] = bsize[bi];
				bsize[bi]++;
			}
		}
	}
	cout << "Populating the bins took: "; tbinpop.printsec(); cout << endl;
	
	Timer tmap(1);
	// maps the counterpart edge (v,u) for each edge (u,v) for equal label bins lu=lv so that the edges can be updated later in O(1) time	
	
	for(i=0; i<=lmax; i++) {
		bi = bidx[i] + i;
		for(j=0; j<bsize[bi]; j++)
			mapEdge[bi][j] = mlist[bin[bi][j].v][mapEdge[bi][j]];
	}
	cout << "Creating and mapping the edges for the bins lu=lv took: "; tmap.printsec(); cout << endl << endl;
	
	/*
	// check whether mapping is correct
	for(i=0; i<=lmax; i++) {
		bi = bidx[i] + i;
		for(j=0; j<bsize[bi]; j++) {
			if((bin[bi][j].u != bin[bi][mapEdge[bi][j]].v) || (bin[bi][j].v != bin[bi][mapEdge[bi][j]].u))
				cout << "i=" << i << ", bi=" << bi << ", j=" << j << ", bin[bi][j].u=" << bin[bi][j].u << ", bin[bi][j].v=" << bin[bi][j].v << ", bin[bi][mapEdge[bi][j]].v=" << bin[bi][mapEdge[bi][j]].v << ", bin[bi][mapEdge[bi][j]].u=" << bin[bi][mapEdge[bi][j]].v << endl;
		}
	}
	*/
	
	cout << "In total, bin preparation took: "; t.printsec(); cout << endl << endl;
	cout << "Shuffling edges ... "; cout.flush();
	
	Vertex u, u2, v2;
	LI vidx, v2idx;
	Elist nbr2;
	srand(time(NULL));
	Timer shuf(1);
	
	//pick s pairs of random edges and swap end points if switchable
	for(bi=0; bi<nbin; bi++) {
		for (i=0; i<nswitch[bi]; ) {
			for(k=0; k<MAX_ATTEMPT; k++) {
					
				// pick the first edge uniformly at random
				LI randnum1 = (LI)rand() % bsize[bi];
				u = bin[bi][randnum1].u;
				v = bin[bi][randnum1].v;
					
				// pick the second edge uniformly at random
				LI randnum2 = (LI)rand() % bsize[bi];
				u2 = bin[bi][randnum2].u;
				v2 = bin[bi][randnum2].v;
				// check for loop and edge switch that does not change the edges
				if (u==u2 || u==v2 || v==u2 || v==v2) {
					repeatedAttempt++;
					continue;
				}
				
				// check for parallel edge (u,v2)
				x = u;
				y = v2;
				bool flag = false;
				for(idx = randnum1-1; idx >= 0 && bin[bi][idx].u == x; idx--) {
					if(bin[bi][idx].v == y) {
						flag = true;
						break;
					}
				}
				if(flag) { repeatedAttempt++; continue; }
				
				for(idx = randnum1+1; idx < bsize[bi] && bin[bi][idx].u == x; idx++) {
					if(bin[bi][idx].v == y) {
						flag = true;
						break;
					}
				}
				if(flag) { repeatedAttempt++; continue; }
				
				// check for parallel edge (u2,v)
				x = u2;
				y = v;
				
				for(idx = randnum2-1; idx >= 0 && bin[bi][idx].u == x; idx--) {
					if(bin[bi][idx].v == y) {
						flag = true;
						break;
					}
				}
				if(flag) { repeatedAttempt++; continue; }
				
				for(idx = randnum2+1; idx < bsize[bi] && bin[bi][idx].u == x; idx++) {
					if(bin[bi][idx].v == y) {
						flag = true;
						break;
					}
				}
				if(flag) { repeatedAttempt++; continue; }
				
				Swap(bin[bi][randnum1].v, bin[bi][randnum2].v);	// swap v and v2
				
				if(eqlLabel[bi]) {
					Swap(bin[bi][mapEdge[bi][randnum1]].v, bin[bi][mapEdge[bi][randnum2]].v);	// swap u and u2
					
					LI temp = mapEdge[bi][randnum1];		// update the mapping info of the edges 
					mapEdge[bi][randnum1] = mapEdge[bi][randnum2];
					mapEdge[bi][mapEdge[bi][randnum2]] = randnum1;
					
					mapEdge[bi][randnum2] = temp;
					mapEdge[bi][temp] = randnum2;
				}
					
				i++;	// increase the number of edge switch
				break;
			}
				
			if(k == MAX_ATTEMPT) {
				//printf("Maximum attempt taken for one edge switch at bin %ld: #edge-at-bin = %ld, #switch-should-be-performed = %ld, #switch-performed = %ld ... Continuing to next bin ...\n", bi, bsize[bi], nswitch[bi], i);
				break;
			}
		}
	}
	
	cout << "done.\n\nDue to loops and parallel edges, repeatedAttempt = " << repeatedAttempt << " (" << double(repeatedAttempt)/s*100 << "%)\n\n";
	cout << "Shuffling edges took: "; shuf.printsec(); cout << endl << endl;
	cout << "Total time to shuffle and bin preparation took: "; total.printsec(); cout << endl << endl; cout.flush();
	
	// write the output graph in edge list format
	if(wflag)
		WriteGraph(ofile, bin, bsize, nbin, eqlLabel);
	
	// free allocated memory
	Timer tf(1);
	FreeAll(bin, nbin);
	FreeAll(mapEdge, nbin);
	FreeMem(bidx);
	FreeMem(bsize);
	FreeMem(eqlLabel);
	FreeMem(nswitch);
	cout << "Freeing the allocated memories in \"Shuffle\" function took: "; tf.printsec(); cout << endl;
}

// write the output graph in edge list format: each edge (u,v) is written as <u v> in one line
void Graph::WriteGraph(const char *ofile, Bin *bin, LI *bsize, LI nbin, bool *eqlLabel)
{
	LI i, k;
	cout << "Writing the shuffled graph: " << ofile << " ... "; cout.flush();
	Timer w(1);
	ofstream ofp(ofile);
	for(i=0; i<nbin; i++)
	{
		if(!eqlLabel[i]) {
			for(k=0; k<bsize[i]; k++)
				ofp << bin[i][k].u << "\t" << bin[i][k].v << endl;
		}
		else {
			for(k=0; k<bsize[i]; k++)
				if(bin[i][k].u < bin[i][k].v)
					ofp << bin[i][k].u << "\t" << bin[i][k].v << endl;
		}
	}
	ofp.close();
	cout << "done.\nWriting the output graph took: "; w.printsec(); cout << endl << endl;
}

/*------------------------------------------------------------------------------------*/
int main(const int argc, const char *argv[])
{    
	Graph g;
	time_t start, sec;
	int i, k;

	if (argc < 6) {
		cout << "USAGE: " << argv[0] << "  <ingraph> <demfile> <outgraph> <Shuffle_Ratio> <Write-output-graph (0/1)>" << endl;
		exit(0);
	}
	
	start = time(NULL);
 	g.ReadGraph(argv[1], argv[2]);
	sec = time(NULL) - start;
	cout << "Time to read graphs: " << sec / 60 << " min. " << sec % 60 << " sec." << endl << endl;
	start = time(NULL);
	g.Shuffle(atof(argv[4]), argv[3], atoi(argv[5]));	// includes shuffle and write the output graph
	sec = time(NULL) - start;
	cout << "Total time including preparing the bins, shuffling and writing the graph: " << sec / 60 << " min. " << sec % 60 << " sec." << endl << endl;
	cout << "Program terminating successfully ..." << endl;
	return 0;
}
