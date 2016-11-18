/* 
 * File	  : Par-Lab-Shuffle.hpp  
 * Purpose: Header file for parallel labeled edge shuffle (MPI program)
 * Author : Hasanuzzaman Bhuiyan
 * Email  : mhb@vbi.vt.edu
 * Date	  : November 14, 2014
 */

#ifndef LABEL_SHUFFLE_H_
#define LABEL_SHUFFLE_H_

#include <map>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <vector>

#include "Serial_Multinomial.hpp"
#include "parsort.hpp"
#include "Timer.hpp"
#include "MinHeap.hpp"
#include "utility.hpp"
#include "TSet.hpp"
#include "mpi.h"

#ifndef NULL
	#define  NULL 0
#endif

#define MSG_LEN 	 		512
#define MAX_LOCL_ATTEMPT 	50
#define MAX_GLOB_ATTEMPT 	20
#define MAX_SERV_ATTEMPT 	20
#define PEDGERSIZE			1024
#define PEDGECSIZE			8
//#define DEBUG_PRINT

#define MANUAL_PARTITION_INPUT

using namespace std;

typedef  long int LI;
//typedef  long int Vertex;
typedef  int 	  Vertex;
typedef  int 	  DT;
typedef  int 	  VSize;
typedef  long int SSize;
typedef  double   LD;


#define MPI_TYPE MPI_INT
// #define MPI_TYPE MPI_LONG


typedef struct {
	Vertex u;
	Vertex v;
} Edge;


typedef struct {
	DT index;	 	// index  of bins [0,noBins-1]
	LI noEdge;	 	// number of bin in this binMetaData
	LI noSwitch; 	// number of bin switches in this binMetaData
	LI Key(){		// sort based on the number of bin switches taking place in this binMetaData
		return noSwitch;
	}
} BinMetaData;


typedef struct{
	BinMetaData binMD;
	int binNo, startProc, endProc;
} GlobalBinMetaData;

typedef struct{
	int sProc, eProc;
} Bin2Proc;


typedef struct{
	DT x, y, indx;
} pEdge;

struct recordTime {
    double 	value;
    int   	index;
} in, in1, out;;


// ==================== EdgeSwitch class =========================

class EdgeSwitch {
	public:
	
		int 		  myRank, noProcs;			// rank of this processor, #total processors
		int			  noGlobProc, noLocProc;	// #processors performing global and local edge switch
		VSize     	  noNodes, noEdges;			// #vertices, #edges in the entire graph
		VSize		  noLocEdges, noLocSwitch, noLocSwitchP;	// #local-edges, #local-switch (should be performed), #local-switch (performed) in this processor respectively
		SSize		  noShuffle, myNoGlobShuffle;	// # of total and my (in case of global shuffle) shuffles in the entire graph
		double		  visitRate, partFactor;	// shuffle fraction, partitioning-factor is used in partitioning
		LI			  noBins, noLocBins, noTotGlobBins, maxLocBinSize;	// total number of bins in the entire graph, number of bins local to this processor
		int 		  maxLabel;					// maximum label
		
		LI			  iledge, fledge, nonZeroLocBin;
		Edge		**bin;						// contains the edges of this partition for efficient random pickup
		vector<Edge>  globBin;
		vector<bool>  isBusy;
		vector<pEdge> pendEdge[PEDGERSIZE];
		double		 *cumProb, noGlobSwitch;	// cumulative edge picking probability array
		LI			 *cumNode;					// only for the case lu == lv: cumulative no. of nodes array containing no. of nodes belonging to processors having nodes with same labels as of partial bins
		
		BinMetaData	 *binMetaData; 				// array holding binMetaData infos. such as index, #edges and #switches in a bin
		GlobalBinMetaData globalBinMetaData;	// a structure holding information of global bin belonging to this processor
		LI 			 *binSize, globBinSize;					// array holding size (#edges) of the edge bins
		vector<Bin2Proc> bin2Proc;				// holds the mapping of which bin belongs to which processor
		vector<int>   binContainer;
		
		DT 			*label, *bsize, *tempbsize, *bidx, *recvBuffer, *sendBufIndex, **sendBuffer, *nodeCount;
		bool 		*nodeInArr;
		LI			 maxSendBufSize, maxRecvBufSize, maxDeg, maxLocDeg, ackCounter;
		LI			 fileSize, byteQuota, byteStart, byteEnd, nextBytePos;
		LI			 compSwitch;
		double		 timeMult, timePart, timeSort, timeShuffle;
		map <DT, DT> binMap;
		MPI_Status   status;
		MPI_Request  request;
		
		// BEGIN: these variables are used for computation in member functions only
		int			 p1, q1, p2, q2, rn;
		int			 Pj, tag, terminateCount, send[2], recv[2];
		// END: these variables are used for computation in member functions only
		
		EdgeSwitch() {
			bin 			= NULL;
			binMetaData 	= NULL;
			binSize 		= NULL;
			cumProb			= NULL;
			noBins = noLocBins = 0;
			noLocEdges = noLocSwitch = noLocSwitchP = globBinSize = 0;
			terminateCount = 0;
			myNoGlobShuffle = 0;
		}
		
		~EdgeSwitch() { FreeGraph(); }
		
		void Initialize(int myrank, int no_of_processors, double v, double f) { 
			myRank 		= myrank;	
			noProcs 	= no_of_processors; 
			visitRate 	= v; 
			partFactor 	= f;
			in.index = in1.index = myRank;
			in.value = in1.value = 0.0L;
		}
		
		void Initialize(int myrank, int no_of_processors, double v, double f, LI max_loc_bsize, double f_no_glob_procs) { 
			myRank 		= myrank;	
			noProcs 	= no_of_processors; 
			visitRate 	= v; 
			partFactor 	= f;
			maxLocBinSize = max_loc_bsize;
			//noGlobProc 	= LI(ceil(f_no_glob_procs * noProcs));
			noGlobProc 	= LI(f_no_glob_procs);
			noLocProc 	= noProcs - noGlobProc;
			in.index = in1.index = myRank;
			in.value = in1.value = 0.0L;
		}
		
		void GraphPreProcessing(const char *graphFileName, const char *demFileName);
			void ReadNoNodes(const char *graphFileName);
			void ReadDemographicFile(const char *demFileName);
			void AllocateDataStructure();
			void ReadGraphToComputeBinSize(const char *graphFileName);
			void MultinomialAndSort();
			void Partition();
				void AssignLocalBins();
				void AssignGlobalBins();
			void ReadGraphToPopulateBins(const char *graphFileName);
				void AllocateDataStructure3();
				void ServeGraphReadMsgsTillSend();
				void ServeGraphReadMsgs();
				void ServeAReadMsg();
			
		void Shuffle();
			void ShufflePreProcessing();
			void ShuffleAtLocalBin();
			void ShuffleAtPartialBin();
				void LocalShuffleAtPartialBin();
				void GlobalShuffleAtPartialBin();
			bool ServeRecvdMessages();
			void SendTerminateSignal();
		
		bool CreateParallelEdgeInGBin(DT index, DT y);
		bool EdgeIsPending(DT x2, DT y1);
		DT   RemoveFromPending(DT x2, DT y1);
		void WriteGraph(const char *oname);
		
		int  PrintAvgMaxMinTime(double t);		
		void FreeGraph();
};
// ======================= END: EdgeSwitch class ===========================


/* =========================== GraphPreProcessing function =============================
 *
 * Reads the demographic-file, input graph, partition the graph and populate the bins
 *
======================================================================================== */
void EdgeSwitch::GraphPreProcessing(const char *graphFileName, const char *demFileName) { // 1st param: original graph name, 2nd param: demographic file name
	
	ReadNoNodes(graphFileName);
	ReadDemographicFile(demFileName);
	AllocateDataStructure();
	ReadGraphToComputeBinSize(graphFileName);
	MultinomialAndSort();
	Partition();
	ReadGraphToPopulateBins(graphFileName);
}


// Read number of vertices from the graph file
void EdgeSwitch::ReadNoNodes(const char *graphFileName) { // 1st param: original graph name
	if(myRank == 0) {
		ifstream ifp(graphFileName);
		if(!ifp.is_open()) {
			cout << "Can not open file: " << graphFileName << endl; cout.flush();
			exit(1);
		}
		ifp >> noNodes;
		ifp.close();
	}
	MPI_Bcast(&noNodes, 1, MPI_LONG, 0, MPI_COMM_WORLD);
}



// Read the demographic file (containing the labels of vertices) in parallel
void EdgeSwitch::ReadDemographicFile(const char *demFileName) { // 1st param: demographic (labeled) file name

	int uLabel, max = 0;
	LI  i, u;
	string  line;
	
	ifstream dfp(demFileName);
	if(!dfp.is_open()) {
		cout << "Can not open file: " << demFileName << endl;
		exit(1);
	}
	
	Timer tread(1), t(1);
	
	label = new int[noNodes];
	int *tempLabel = new int[noNodes];
	for(i=0; i<noNodes; i++)
		label[i] = tempLabel[i] = 0;
	
	in.value += t.getsec();
	
	if(myRank == 0) { cout << "Reading demographic file: " << demFileName << " ... "; cout.flush(); }
	
	dfp.seekg(0, dfp.end);
    fileSize  = dfp.tellg();
	byteQuota = LI(ceil(double(fileSize)/noProcs));
	byteStart = myRank * byteQuota;
	byteEnd   = byteStart + byteQuota;	// for last processor, it may exceed fileSize, hence need to use ifp.good()
	dfp.clear();
	dfp.seekg(byteStart, dfp.beg);
	
	if(myRank > 0 && dfp.good())	// may fall in broken line
		getline(dfp, line);
	
	byteStart = dfp.tellg();
	
	// sending msg containing byte address up to which the previous processor should read
	if(myRank > 0)	// except first processor
		MPI_Send(&byteStart, 1, MPI_LONG, myRank-1, 99, MPI_COMM_WORLD);	// sending the start address of graph read of this processor to previous processor
	
	if(myRank != noProcs-1)	// except last processor
		MPI_Recv(&byteEnd, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	else
		byteEnd = Min(byteEnd, fileSize);
	
	dfp.clear();
	dfp.seekg(byteStart, dfp.beg);
	nextBytePos = byteStart;
	
	while(nextBytePos < byteEnd && dfp.good()) {
		dfp >> u >> uLabel;
		tempLabel[u] = uLabel;
		max = Max(max, uLabel);
		nextBytePos = 1 + dfp.tellg();	// adding 1 required for the new-line/carriage-return
	}
	
	dfp.close();
	MPI_Allreduce(tempLabel, label, noNodes, MPI_INT, MPI_SUM, MPI_COMM_WORLD);	// make the labels of all nodes available at each processor
	MPI_Allreduce(&max, &maxLabel, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);		// maxLabel = maximum label
	
	FreeMem(tempLabel);
	if(myRank == 0) { cout << "done.\nReadDemographicFile: ";}
	PrintAvgMaxMinTime(tread.getsec());
}


// allocate and initialize necessary data structures
void EdgeSwitch::AllocateDataStructure() {

	LI i, j, k;
	Timer t(1);
	bidx = new DT[maxLabel+1];	// additional 1 for label zero
	bidx[0] = 0;
	for (i=1; i<= maxLabel; i++)// bidx is used to save time for computing index 
		bidx[i] = bidx[i-1] + i;

	noBins	  = bidx[maxLabel] + maxLabel + 1;	// number of bins
	bsize	  = new DT[noBins];
	tempbsize = new DT[noBins];
	
	for(i=0; i<noBins; i++)
		tempbsize[i]  = 0;
	in.value += t.getsec();
	if(myRank == 0) printf("AllocateDataStructure: ");
	PrintAvgMaxMinTime(t.getsec());
}


// READ graph file in parallel to count the number and size of the bins ...
void EdgeSwitch::ReadGraphToComputeBinSize(const char *graphFileName) { // 1st param: original graph name
	
	DT  u, lu, v, lv;
	LI  i, deg, wordCount;
	LD  	number;
	string  line;
	istringstream istr;
	
	ifstream ifp(graphFileName);
	if(!ifp.is_open()) {
		cout << "Can not open file: " << graphFileName << endl;
		exit(1);
	}
	
	Timer tread(1);
	if(myRank == 0) { cout << "Reading graph file (1st phase--count binsize): " << graphFileName << " ... "; cout.flush(); }
	
	ifp.seekg(0, ifp.end);
    fileSize  = ifp.tellg();
	byteQuota = LI(ceil(double(fileSize)/noProcs));
	byteStart = myRank * byteQuota;
	byteEnd   = byteStart + byteQuota;	// for last processor, it may exceed fileSize, hence need to use ifp.good()
	ifp.clear();
	ifp.seekg(byteStart, ifp.beg);
	nextBytePos = byteStart;
	maxLocDeg = 0;
	
	if(ifp.good())	// may fall in broken line
		getline(ifp, line);
	
	while(ifp.good()) {
		byteStart = ifp.tellg();
		getline(ifp, line);
		istr.clear();
		istr.str(line);
		
		wordCount = 0;
		while(istr >> number)
			wordCount++;
		if(wordCount == 2)	// Found the beginning of an adjacency list of a vertex in Galib format: <vertexID> <Degree>
			break;
	}
	
	// sending msg containing byte address up to which the previous processor should read
	if(myRank > 0)	// except first processor
		MPI_Send(&byteStart, 1, MPI_LONG, myRank-1, 89, MPI_COMM_WORLD);	// sending the start address of graph read of this processor to previous processor
	
	if(myRank != noProcs-1)	// except last processor
		MPI_Recv(&byteEnd, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	else
		byteEnd = Min(byteEnd, fileSize);
	
	ifp.clear();
	ifp.seekg(byteStart, ifp.beg);
	nextBytePos = byteStart;
	
	while(nextBytePos < byteEnd) {
		ifp >> u >> deg;
		if(!ifp.good() || ifp.tellg() >= byteEnd)	// may reach EOF
			break;
		lu  = label[u];
		maxLocDeg = Max(maxLocDeg, deg);
		for(i=0; i < deg; i++) {
			ifp >> v;
			lv = label[v];
			if(lu > lv || (lu == lv && u > v))
				tempbsize[bidx[lu]+lv]++;
		}
		nextBytePos = 1 + ifp.tellg();	// adding 1 for the new-line/carriage-return
	}
	// READ graph file complete ...
	
	// accumulate the bin size across all processors
	MPI_Allreduce(tempbsize, bsize, noBins, MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&maxLocDeg, &maxDeg, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);			// maxDeg = maximum degree of a node
	FreeMem(tempbsize);
	
	if(myRank == 0) { cout << "done.\nReadGraphToCountBinSize: "; }
	PrintAvgMaxMinTime(tread.getsec());
}


void EdgeSwitch:: MultinomialAndSort() {
	
	LI i;
	Timer t(1);
	noEdges = 0;
	for(i=0; i<noBins; i++)
		noEdges += bsize[i];
	
	LD *prob 		= new LD[noBins];
	LI *tempNSwitch = new LI[noBins];
	
	// for parallel sorting purpose, array size should be a multiple of <noProcs>
	LI   ceilNBin	= LI(ceil(1.0L*noBins/noProcs)*noProcs);
	LI  *noSwitch  	= new LI[ceilNBin];
	int *index 		= new int[ceilNBin];
	
	// calculate the probability of an edge switch operation at each edge bin
	for (i=0; i<noBins; i++) {
		prob[i]  = LD(bsize[i])/noEdges;
		index[i] = i;
	}
	
	// determine the number of samples
	if (visitRate < 0.9999)
		noShuffle = LI( - double(noEdges)/2.0 * log(1.0 - visitRate) );
	else
		noShuffle = LI( double(noEdges)/2.0 * log(double(noEdges)) );
	in.value += t.getsec();
	
	if(myRank == 0) {
		cout << "|V| = " << noNodes << ", |E| = " << noEdges << ", fraction = " << visitRate << ", #shuffle = " << noShuffle << ", #procs = " << noProcs << endl;
		cout << "Max label = " << maxLabel << ", No. of bins = " << noBins << ", maxLocDeg (at proc. 0) = " << maxLocDeg << ", maxDeg = " << maxDeg << endl << endl; cout.flush();
	}
	
	// MULTINOMIAL distribution ...
	// determine which bin performs how many edge switch operations by generating multinomial random variables
	MPI_Barrier(MPI_COMM_WORLD);	// not necessary though ... using it to measure time spent properly ...
	if(myRank == 0) { cout << "Generating multinomial random variables ... "; cout.flush(); }
	t.start();
	LI multAtThisProc = (myRank < noShuffle%noProcs) ? noShuffle/noProcs + 1 : noShuffle/noProcs;
	Multinomial(multAtThisProc, prob, noBins, tempNSwitch);
	MPI_Allreduce(tempNSwitch, noSwitch, noBins, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	in.value += t.getsec();
	in1.value += t.getsec();	// multinomial time
	timeMult = t.getsec();
	if(myRank == 0) {
		// OPTIONAL: check whether multinomial random variables were generated successfully
		LI dummys = 0;
		for(i=0; i<noBins; i++)
			dummys += noSwitch[i];
		if(dummys != noShuffle) {
			cout << "Error in multinomial distribution ...\n";
			exit(1);
		}
		
		cout << "done.\nMultinomial: ";
	}
	
	PrintAvgMaxMinTime(t.getsec());
	
	t.start();
	FreeMem(prob);	// prob is no longer needed ... so free it
	FreeMem(tempNSwitch);
	
	// for parallel sorting purpose, fillup the rest of the array
	for(i=noBins; i<ceilNBin; i++) {
		noSwitch[i] = 0;
		index[i]   	= i;
	}
	in.value  += t.getsec();
	
	// sort in descending order ... based on number of bin switches (noSwitch)
	if(myRank == 0) { cout << "Sorting ... "; cout.flush(); }
	MPI_Barrier(MPI_COMM_WORLD);	// not necessary though ... using it to measure time spent properly ...
	Timer tsort(1);
	ParMergeSortD(noSwitch, index, int(ceilNBin), myRank, noProcs);	// parallel merge sort in descending order
	
	in.value  += tsort.getsec();
	in1.value += tsort.getsec();	// sorting time
	timeSort   = tsort.getsec();
	if(myRank == 0)	{ cout<< "done.\nParallel-merge-sort: "; }
	PrintAvgMaxMinTime(tsort.getsec());
	
	Timer tbcast(1);
	MPI_Bcast(noSwitch, noBins, MPI_LONG, 0, MPI_COMM_WORLD);	// root - proc 0, since the merged sorted list is at proc 0
	MPI_Bcast(index, noBins, MPI_INT, 0, MPI_COMM_WORLD);		// root - proc 0, since the merged sorted list is at proc 0
	in.value += tbcast.getsec();
	//in1.value += tbcast.getsec();
	
	if(myRank == 0)	cout<< "Message broadcasting: ";
	PrintAvgMaxMinTime(tbcast.getsec());
	if(myRank == 0)	cout<< "Sorting in total (par-sort + message broadcast): ";
	PrintAvgMaxMinTime(tsort.getsec() + tbcast.getsec());
	
	t.start();
	binMetaData = new BinMetaData[noBins]; 	// keeps metadata of bins ...
	bin2Proc.reserve(noBins);		 		// keeps the mapping info of which bin belongs to which processor range [sProc, eProc]
	binContainer.reserve(LI(ceil(1.25*noBins/noProcs)));
	
	for(i=0; i<noBins; i++) {
		binMetaData[i].noSwitch = noSwitch[i];
		binMetaData[i].index 	= index[i];
		binMetaData[i].noEdge 	= bsize[index[i]];
	}
	
	FreeMem(index);
	FreeMem(noSwitch);
	FreeMem(bsize);
	in.value += t.getsec();
	if(myRank == 0) cout << "AllocateDataStructure2: ";
	PrintAvgMaxMinTime(t.getsec());
}



void EdgeSwitch::Partition() {
	
	LI i;
	
	if(myRank == 0) cout << "partitioning the bins ...\n\n";
	
	noGlobSwitch  = 0.0;				// no. of global switch
	#ifndef MANUAL_PARTITION_INPUT
		maxLocBinSize = (LI)ceil(partFactor*(double(noShuffle)/noBins));	// maximum local bin size in terms of edge switch
		//maxLocBinSize = (LI)ceil(partFactor*(double(noShuffle)/noProcs));	// maximum local bin size in terms of edge switch
	#endif
	
	// compute required no. of processors for edge switches in the bins that has #switches(per bin) > max no. of switch per processor
	// each of these bins will be partitioned
	if(maxLocBinSize >= binMetaData[0].noSwitch && noBins < noProcs) {
		if(myRank == 0) printf("Since maxLocBinSize = %ld is greater than binMetaData[0].noSwitch = %ld, partFactor is set to 1\n", maxLocBinSize, binMetaData[0].noSwitch);
		//maxLocBinSize = LI(ceil(double(noShuffle)/noProcs));
		maxLocBinSize = (LI)ceil((double(noShuffle)/noBins));
		if(myRank == 0) printf("New maxLocBinSize = %ld\n", maxLocBinSize);
	}
	if(myRank == 0) printf("maxLocBinSize = %ld, binMetaData[0].noSwitch = %ld\n", maxLocBinSize, binMetaData[0].noSwitch);
	
	for(i=0; i<noBins; i++)	// bins [0,i-1] will be partitioned and considered as global bins requiring communication. bins [i,noBins-1] will be local bins requiring no communication
		if(binMetaData[i].noSwitch <= maxLocBinSize)
			break;
		else
			noGlobSwitch += binMetaData[i].noSwitch;
	
	if(myRank == 0) cout << "\n\ni = " << i << ", noGlobSwitch = " << noGlobSwitch << " (" << noGlobSwitch/noShuffle*100 << "%), noLocSwitch = " << noShuffle - noGlobSwitch << ", maxLocBinSize = " << maxLocBinSize << "\n\n";
	
	//// WILL NEED POLICY TO DETERMINE <noLocProc> and <noGlobProc> automatically
	if(i == 0) {
		noLocProc  = noProcs;
		noGlobProc = 0;
	}
	else {
		#ifndef MANUAL_PARTITION_INPUT
			noLocProc = (noShuffle - noGlobSwitch) > 0 ? Max((LI)ceil((noShuffle - noGlobSwitch)/maxLocBinSize), 1L) : 0;
			//noLocProc = int(noProcs/2);
			//noLocProc  = 40;
			noGlobProc = noProcs - noLocProc;
		#endif
	}
	
	// Will it be required in new strategy???
	if(i > 0 && noGlobProc < 2) {
		i = 0;
		noLocProc  = noProcs;
		noGlobProc = 0;
		//printf("There are not enough number of processors to perform global edge switch\ni = %ld, noGlobProc = %ld\n", i, noGlobProc);
		//exit(1);
	}
	noTotGlobBins 	= i;
	
	if(myRank == 0) printf("noTotGlobBins = %ld, noLocProc = %d, splitAtProc = noGlobProc = %d, noGlobSwitch = %.0lf\n", noTotGlobBins, noLocProc, noGlobProc, noGlobSwitch);
	
	MPI_Barrier(MPI_COMM_WORLD);	// used for debugging purpose;  not necessary though
	Timer tPartition(1);
	
	AssignLocalBins();
	//cout.flush();					// used for debugging purpose;  not necessary though
	//MPI_Barrier(MPI_COMM_WORLD);	// used for debugging purpose;  not necessary though
	AssignGlobalBins();
	
	in.value  += tPartition.getsec();
	in1.value += tPartition.getsec();
	timePart   = tPartition.getsec();
	if(myRank == 0) cout << "done.\nPartition: ";
	PrintAvgMaxMinTime(tPartition.getsec());
}



void EdgeSwitch::AssignLocalBins() {
	
	LI i, hpsz = noLocProc;			// heap-size
	MinHeap<LI, LI> minHeap(hpsz);	// create a min-heap for the remaining number of processors
	
	for(i=0; i<hpsz; i++) {
		minHeap.heap[i]  = i;	// processor's RELATIVE rank
		minHeap.value[i] = 0;	// initial #switches at each processor
		minHeap.index[i] = i;
	}
		
	minHeap.heapsize = minHeap.length-1;
	//minHeap.MakeHeap();	// not necessary since all of the initial values are 0
	
	for(i=noTotGlobBins; i<noBins; i++) {
		if(binMetaData[i].noEdge == 0)
			continue;
		LI relativeProcRank = minHeap.FindMinimum();	// find the processor with minimum no of edge switches in the minheap
		int realProcRank	= relativeProcRank + noGlobProc;
		bin2Proc[binMetaData[i].index] = (Bin2Proc){realProcRank, realProcRank};
		
		//if(myRank == 0) cout << "Assigning local bin " << i << ", index = " << binMetaData[i].index << ", noSwitch = " << binMetaData[i].noSwitch << ", noEdge = " << binMetaData[i].noEdge << " to processor " << realProcRank << " with value before assignment = " << minHeap.value[relativeProcRank] << endl;
		
		if(myRank == realProcRank) {
			binContainer.push_back(int(i));
			noLocBins++;
		}
		minHeap.IncreaseKey(relativeProcRank, minHeap.value[relativeProcRank]+binMetaData[i].noSwitch);	
	}
}



void EdgeSwitch::AssignGlobalBins() {
	
	int curProc = 0, availNoGlobProc = noGlobProc;
	LI remnNoGlobSwitch = LI(noGlobSwitch);
	
	if(availNoGlobProc <= noTotGlobBins) {
		if(myRank == 0) printf("Error inside AssignGlobalBins function ... \navailNoGlobProc <= noTotGlobBins, availNoGlobProc = %d, noTotGlobBins = %ld\n", availNoGlobProc, noTotGlobBins);
		exit(1);
	}
	
	//if(myRank == 0) printf("noTotGlobBins = %ld, noGlobSwitch = %.0lf, availNoGlobProc = %d\n", noTotGlobBins, noGlobSwitch, availNoGlobProc);
	
	for(int i = 0; i < noTotGlobBins; i++) {
		// availNoGlobProc can become 0!
		// enough number of processors left for the later bins!!!
		// reqNoProc --> take floor/ceil/round?
		if(availNoGlobProc <= 0) {
			if(myRank == 0) printf("Error inside AssignGlobalBins function ... \navailNoGlobProc = %d while there are unassigned global bins, #unassignedNoGlobBins = %ld whereas #assignedNoGlobBins = %d\n", availNoGlobProc, noTotGlobBins-i, i);
			exit(1);
		}
		
		LI idealNoGlobSwitch = (LI)ceil(double(remnNoGlobSwitch)/availNoGlobProc);
		int reqNoProc 	  = int(round(double(binMetaData[i].noSwitch)/idealNoGlobSwitch));	// required no. of processors for this bin
		int sProc 		  = curProc;		// this bin is assigned from processor <sProc> to <eProc>
		int eProc 		  = curProc + reqNoProc - 1;
		//if(myRank == sProc) printf("\n\nremnNoGlobSwitch = %ld, remnNoGlobBin = %ld, availNoGlobProc = %d, idealNoGlobSwitch = %ld\nBinMetaData %ld (index = %d, noEdge = %ld, noSwitch = %ld) is assigned from proc %d to %d with perProcSwitch = %ld, perProcEdge = %ld\n", remnNoGlobSwitch, noTotGlobBins-i, availNoGlobProc, idealNoGlobSwitch, i, binMetaData[i].index, binMetaData[i].noEdge, binMetaData[i].noSwitch, sProc, eProc, binMetaData[i].noSwitch/reqNoProc, binMetaData[i].noEdge/reqNoProc);
		
		bin2Proc[binMetaData[i].index] = (Bin2Proc){sProc,eProc};
		if(myRank >= sProc && myRank <= eProc)
			globalBinMetaData = (GlobalBinMetaData){binMetaData[i], i, sProc, eProc};
		
		curProc 			 = eProc + 1;
		remnNoGlobSwitch 	-= binMetaData[i].noSwitch;
		availNoGlobProc		-= reqNoProc;
	}
}



void EdgeSwitch::ReadGraphToPopulateBins(const char *graphFileName) {
	
	DT		u, v, proc;
	LI 		i, j, k, deg, lu, lv, bi, li, lj ;
	
	ifstream ifp(graphFileName);
	if(!ifp.is_open()) {
		cout << "Can not open file: " << graphFileName << endl;
		exit(1);
	}
	
	AllocateDataStructure3();
	
	// Reading the graph file (second phase) ...
	if(myRank == 0) { cout << "maxRecvBufSize = " << maxRecvBufSize*2 << "\nReading graph file (2nd phase--populate the bins): " << graphFileName << " ... "; cout.flush(); }
	Timer tread(1);
	
	ifp.seekg(byteStart, ifp.beg);
	nextBytePos = byteStart;
	
	while(nextBytePos < byteEnd) {
		ifp >> u >> deg;
		if(!ifp.good() || ifp.tellg() >= byteEnd)	// may reach EOF
			break;

		lu  = label[u];
		
		for(i=0; i<noProcs; i++)
			nodeInArr[i] = false;
		
		for(i=0; i < deg; i++) {
			ifp >> v;
			lv = label[v];
			
			if(lu > lv || (lu == lv && u > v)) {
				li 		  	 = bidx[lu] + lv;
				int sProc 	 = bin2Proc[li].sProc;
				int eProc 	 = bin2Proc[li].eProc;
				int reqNProc = eProc - sProc + 1;
				proc 		 = sProc + u % reqNProc;
				
				//printf("Proc %d is assigning edge (%d,%d) with label (%d,%d) and li = %d and to processor %d\n", myRank, u, v, lu, lv, li, proc);
				
				if(proc == myRank) {	// edge (u,v) belongs to this processor
					if(noLocBins > 0) {
						bi = binMap.at(li);
						bin[bi][binSize[bi]++] = (Edge){u,v};
					}
					else {
						globBin.push_back((Edge){u,v});
						isBusy.push_back(false);
					}
				} 
				else {	// edge (u,v) belongs to other processor
					if(!nodeInArr[proc]) {		// putting the node <u>'s adj. edges in the send array for the first time!
						nodeInArr[proc] = true;
						sendBuffer[proc][0]++;							// one more node <u>
						sendBuffer[proc][sendBufIndex[proc]++] = u;		// put node <u>
						nodeCount[proc] = sendBufIndex[proc];			// save the index which holds the number of edges adjacenct to <u>
						sendBuffer[proc][sendBufIndex[proc]++] = 0;		// number of edges adjacenct to <u>
					}
					sendBuffer[proc][sendBufIndex[proc]++] = v;
					sendBuffer[proc][nodeCount[proc]]++;
				}
			}
		}
		
		for(proc=0; proc<noProcs; proc++) {
			if(maxSendBufSize - sendBufIndex[proc] < 2)	{ // may not have enough space
				MPI_Isend(sendBuffer[proc], sendBufIndex[proc], MPI_TYPE, proc, 81, MPI_COMM_WORLD, &request);
				ServeGraphReadMsgsTillSend();
				sendBuffer[proc][0] = 0;
				sendBufIndex[proc]  = 1;
			}
		}
		
		int flag = 1;
		while(flag) {
			MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
			if(flag)
				ServeGraphReadMsgs();
		}
		
		nextBytePos = 1 + ifp.tellg();	// adding 1 for the new-line/carriage-return
	}
	
	ifp.close();
	
	for(i=0; i<noProcs; i++) {	// send the remaining edges
		if(sendBuffer[i][0] && i != myRank) {
			MPI_Isend(sendBuffer[i], sendBufIndex[i], MPI_TYPE, i, 81, MPI_COMM_WORLD, &request);
			ServeGraphReadMsgsTillSend();
		}
	}
	
	for(i=0; i<noProcs; i++) {	// send acknowledgement
		if(i != myRank)	{
			MPI_Isend(sendBuffer[i], 1, MPI_TYPE, i, 85, MPI_COMM_WORLD, &request);		// acknowledgement
			ServeGraphReadMsgsTillSend();
		}
	}
	
	while(ackCounter != noProcs-1)	// until receiving ack from every processor, receive messages from other processors
		ServeGraphReadMsgs();
	
	if(myRank == 0) { cout << "done.\nReadGraphToPopulateBins: "; }
	PrintAvgMaxMinTime(tread.getsec());
	// Reading the graph file (second phase) complete ...
	
	//if(!noLocBins) printf("Proc %d, GlobalBinsize(or #edges) = %ld\n", myRank, LI(globBin.size()));
	
	/*if(noLocBins) {
		for(i=0; i<noLocBins; i++) {	// print edges in the local bins
			printf("Proc %d, i = %ld, binsize = %ld, edges are: ", myRank, i, binSize[i]);
			for(j=0; j<binSize[i]; j++)
				printf("(%d,%d), ", bin[i][j].u, bin[i][j].v);
			printf("\n");
		}
	}
	else {	// print edges in the global bins
		printf("Proc %d, GlobalBinsize = %ld, edges are: ", myRank, LI(globBin.size()));
		for(j=0; j<globBin.size(); j++)
			printf("(%d,%d), ", globBin[j].u, globBin[j].v);
		printf("\n");
	}*/
	
	Timer t(1);
	FreeMem(label);
	FreeMem(bidx);
	FreeMem(binSize);
	FreeMem(nodeCount);
	FreeMem(nodeInArr);
	FreeMem(recvBuffer);
	FreeMem(sendBufIndex);
	FreeAll(sendBuffer, noProcs);
	vector<Bin2Proc>().swap(bin2Proc);
	binMap.clear();
	in.value += t.getsec();
}


void EdgeSwitch::AllocateDataStructure3() {
	
	Timer t(1);
	BinMetaData *tmpbin = binMetaData;
	
	if(myRank >= noGlobProc) {	// processors performing local edge switch operations ...
		binSize 	= new LI[noLocBins];
		bin 		= new Edge*[noLocBins];	// allocate the edge bins
		binMetaData = new BinMetaData[noLocBins];
		for(LI i=0; i<noLocBins; i++) {
			binMetaData[i] = tmpbin[binContainer[i]];
			binMap[binMetaData[i].index] = i;
			bin[i]		= new Edge[binMetaData[i].noEdge];
			binSize[i] 	= 0;
		}
		vector<int>().swap(binContainer);
	}
	else {
		binMetaData   = NULL;
		int reqNoProc = globalBinMetaData.endProc - globalBinMetaData.startProc + 1;
		LI bSize 	  = LI(1.25*ceil(LD(globalBinMetaData.binMD.noEdge)/reqNoProc));
		globBin.reserve(bSize);
		isBusy.reserve(bSize);
	}
	FreeMem(tmpbin);
	
	ackCounter 	 	= 0;
	maxSendBufSize 	= 2 * Max(maxLocDeg + 3, MSG_LEN);	// maximum send 	buffer size, extra 3 because of #node, node-id, #edges-adjacenct-to-that-node, multiplied by 2 so that the next node's adjacency list (in worst case, the entire list) fits into the array
	maxRecvBufSize 	= 2 * Max(maxDeg + 3, MSG_LEN);		// maximum receive buffer size, extra 3 because of #node, node-id, #edges-adjacenct-to-that-node
	recvBuffer	 	= new DT[maxRecvBufSize];
	sendBufIndex 	= new DT[noProcs];
	sendBuffer  	= new DT*[noProcs];
	nodeInArr  		= new bool[noProcs];
	nodeCount   	= new DT[noProcs];
	
	for(int i=0; i<noProcs; i++) {
		sendBuffer[i] 	 = new DT[maxSendBufSize];	// buffering messages that will be sent to processor <i>
		sendBuffer[i][0] = 0;	// how many edges will be sent to processor <i>
		sendBufIndex[i]  = 1;	// next available index/slot in the buffer for processor <i>
	}
	maxSendBufSize /= 2;
	maxRecvBufSize /= 2;
	
	in.value += t.getsec();
	if(myRank == 0) printf("AllocateDataStructure3: ");
	PrintAvgMaxMinTime(t.getsec());
}



void EdgeSwitch::ServeGraphReadMsgsTillSend() {
	
	int sendFlag = 0;
	while(!sendFlag) {	// until true meaning safe to use the send buffer again
		MPI_Test(&request, &sendFlag, &status);
		int recvFlag = 1;
		while(recvFlag) {
			MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &recvFlag, &status);
			if(recvFlag)
				ServeGraphReadMsgs();
		}
	}
}


void EdgeSwitch::ServeGraphReadMsgs() {
	
	MPI_Recv(recvBuffer, maxRecvBufSize*2, MPI_TYPE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	if(status.MPI_TAG == 85)
		ackCounter++;
	else
		ServeAReadMsg();
}


void EdgeSwitch::ServeAReadMsg() {
	
	LI nextNodeIndex = 1;
	for(LI i=0; i<recvBuffer[0]; i++) {
		DT u    	 = recvBuffer[nextNodeIndex++];
		DT noAdjEdge = recvBuffer[nextNodeIndex++];
		DT lu   	 = label[u];
		
		for(LI j=0; j<noAdjEdge; j++) {
			DT v  = recvBuffer[nextNodeIndex++];		
			if(noLocBins > 0) {
				DT lv = label[v];
				LI li = bidx[lu] + lv;
				LI bi = binMap.at(li);
				bin[bi][binSize[bi]++] = (Edge){u,v};
			}
			else {
				globBin.push_back((Edge){u,v});
				isBusy.push_back(false);
			}
		}
	}
}


void EdgeSwitch::Shuffle() {
	
	ShufflePreProcessing();
	if(myRank == 0) cout << "Switching edges ... ";
	MPI_Barrier(MPI_COMM_WORLD);	// used to determine precise shuffle time;  not necessary though
	
	if(noLocBins != 0)
		ShuffleAtLocalBin();
	else
		ShuffleAtPartialBin();
	
	if(myRank == 0) cout << "done.\n";
	
	LI tempNoSwitch, noSwitchPerformed = noLocSwitchP + myNoGlobShuffle;
	MPI_Reduce(&noSwitchPerformed, &tempNoSwitch, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	if(myRank == 0) { cout << "No. edge switch performed = " << tempNoSwitch << " (" << tempNoSwitch*1.0/noShuffle*100 << "%) out of total noShuffle = " << noShuffle << endl << endl; cout.flush(); }	
	
	if(myRank == 0) cout << "Total computation time: ";
	PrintAvgMaxMinTime(in.value);
	
	if(myRank == 0) cout << "Computation time for <multinomial + sort + partition + shuffle> (excluding bin pre-process, memory allocation, deallocation): ";
	int p = PrintAvgMaxMinTime(in1.value);
	
	fflush(stdout);
	cout.flush();
	
	if(myRank == p) {
		cout << "Summary of time (sec.) taken by Proc " << myRank << ":\n\n\tMultinomial:\t" << timeMult << endl;
		cout << "\tSort:\t\t" << timeSort << endl;
		cout << "\tPartition:\t" << timePart << endl;
		cout << "\tShuffle:\t" << timeShuffle << endl << endl;
	}
	
	fflush(stdout);
	cout.flush();
	
	if(myRank == 0) cout << "Multinomial: ";
	PrintAvgMaxMinTime(timeMult);
	if(myRank == 0) cout << "Sort: ";
	PrintAvgMaxMinTime(timeSort);
	if(myRank == 0) cout << "Partition: ";
	PrintAvgMaxMinTime(timePart);
	if(myRank == 0) cout << "Shuffle: ";
	PrintAvgMaxMinTime(timeShuffle);
}


void EdgeSwitch::ShufflePreProcessing() {

	int i, reqNoProc, *ranks, tempRank = myRank, sProc, *recvCounts;
	double p, *prob;
	
	Timer t(1);
	MPI_Group origGroup, newGroup, dummyGroup;
	MPI_Comm  newComm, dummyComm;
	MPI_Comm_group(MPI_COMM_WORLD, &origGroup);
	
	if(noLocBins == 0) {	// global processors containing global bin ...
		compSwitch 	= 0;
		sProc		= globalBinMetaData.startProc;
		reqNoProc 	= globalBinMetaData.endProc - sProc + 1;
		prob		= new double[reqNoProc];
		ranks 		= new int[reqNoProc];
		recvCounts 	= new int[reqNoProc];
		globBinSize = globBin.size();
		
		for(i=0; i<PEDGERSIZE; i++)
			pendEdge[i].reserve(PEDGECSIZE);
		
		for(i=0; i<reqNoProc; i++) {
			ranks[i] 	  = i + sProc;
			recvCounts[i] = 1;
		}
	}
	
	for(int i=0; i<noTotGlobBins; i++) {
		if(noLocBins == 0 && globalBinMetaData.binNo == i) {
			MPI_Group_incl(origGroup, reqNoProc, ranks, &newGroup);
			MPI_Comm_create(MPI_COMM_WORLD, newGroup, &newComm);
		}
		else {
			MPI_Group_incl(origGroup, 1, &tempRank, &dummyGroup);
			MPI_Comm_create(MPI_COMM_WORLD, dummyGroup, &dummyComm);
		}
	}
	
	if(noLocBins == 0) {	// global processors containing global bin ...
		p = LD(globBin.size())/globalBinMetaData.binMD.noEdge;
		MPI_Allgather(&p, 1, MPI_DOUBLE, prob, 1, MPI_DOUBLE, newComm);
		
		LI  *temp 	 = new LI[reqNoProc];
		LI no_switch = globalBinMetaData.binMD.noSwitch;
		no_switch 	 = (myRank - sProc < no_switch%reqNoProc) ? no_switch/reqNoProc + 1 : no_switch/reqNoProc;
		Multinomial(no_switch, prob, reqNoProc, temp);
		MPI_Reduce_scatter(temp, &myNoGlobShuffle, recvCounts, MPI_LONG, MPI_SUM, newComm);
		
		cumProb = prob;
		for(i=1; i<reqNoProc; i++)
			cumProb[i] += cumProb[i-1];
		
		FreeMem(temp);
		FreeMem(ranks);
		FreeMem(recvCounts);
	}
	in.value += t.getsec();
}



void EdgeSwitch::ShuffleAtLocalBin() {
	
	LI i, k, bi, idx, x, y, randnum1, randnum2;
	Vertex u, v, u2, v2;
	Timer t(1);
	
	for(bi=0; bi < noLocBins; bi++) {
		
		noLocEdges 	+= binMetaData[bi].noEdge;
		noLocSwitch += binMetaData[bi].noSwitch;
		LI bSize 	 = binMetaData[bi].noEdge;
		
		for (i=0; i < binMetaData[bi].noSwitch; ) {
			for(k=0; k < MAX_LOCL_ATTEMPT; k++) {
				
				// pick the first edge uniformly at random
				randnum1 = (LI)rand() % bSize;
				u = bin[bi][randnum1].u;
				v = bin[bi][randnum1].v;
			
				// pick the second edge uniformly at random
				randnum2 = (LI)rand() % bSize;
				u2 = bin[bi][randnum2].u;
				v2 = bin[bi][randnum2].v;
				
				// check for loop and edge switch that does not change the edges
				if (u==u2 || u==v2 || v==u2 || v==v2)
					continue;
				
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
				if(flag) continue;
				
				for(idx = randnum1+1; idx < bSize && bin[bi][idx].u == x; idx++) {
					if(bin[bi][idx].v == y) {
						flag = true;
						break;
					}
				}
				if(flag) continue;
				
				// check for parallel edge (u2,v)
				x = u2;
				y = v;				
				for(idx = randnum2-1; idx >= 0 && bin[bi][idx].u == x; idx--) {
					if(bin[bi][idx].v == y) {
						flag = true;
						break;
					}
				}
				if(flag) continue;
				
				for(idx = randnum2+1; idx < bSize && bin[bi][idx].u == x; idx++) {
					if(bin[bi][idx].v == y) {
						flag = true;
						break;
					}
				}
				if(flag) continue;
				
				Swap(bin[bi][randnum1].v, bin[bi][randnum2].v);	// swap v and v2
				
				i++;	// increase the number of edge switch
				break;
			}
			
			if(k == MAX_LOCL_ATTEMPT) {
				//printf("Maximum attempt taken for one edge switch at bin %ld: #edge-at-bin = %ld, #switch-should-be-performed = %ld, #switch-performed = %ld ... Continuing to next bin ...\n", binMetaData[bi].index, binMetaData[bi].noEdge, binMetaData[bi].noSwitch, i);
				break;
			}
		}
		
		noLocSwitchP += i;
	}
	
	in.value  	+= t.getsec();
	in1.value 	+= t.getsec();
	timeShuffle  = t.getsec();
}


void EdgeSwitch::ShuffleAtPartialBin() {
	
	Timer t(1);
	int sProc 		= globalBinMetaData.startProc;
	int reqNoProc 	= globalBinMetaData.endProc - globalBinMetaData.startProc + 1;
	
	for(LI k, ni = 0; ni < myNoGlobShuffle; ) {
		for(k=0; k < MAX_GLOB_ATTEMPT; k++) {
			Pj  = sProc + BRSearch(cumProb, 0, reqNoProc-1, double(rand())/RAND_MAX);
			tag = 0;
			
			if(Pj == myRank)
				LocalShuffleAtPartialBin();
			else
				GlobalShuffleAtPartialBin();
			
			if(tag%10 != 3) {	// successful edge switch
				ni++;
				break;
			}
		}
		
		if(k == MAX_GLOB_ATTEMPT) {
			cout << "Maximum attempt taken at processor " << myRank << " for a global edge switch operation at bin " << globalBinMetaData.binMD.index << endl;
			cout << "No. of edge switch operation performed = " << ni << " whereas expected #switch operation was = " << myNoGlobShuffle << endl;
			exit(1);
		}
	}
	
	SendTerminateSignal();
	
	while(terminateCount != reqNoProc-1) {	// receive terminate signal from every proc in the group except myself
		MPI_Recv(recv, 2, MPI_TYPE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		ServeRecvdMessages();
	}
	
	in.value 	+= t.getsec();
	in1.value 	+= t.getsec();
	timeShuffle  = t.getsec();
}


void EdgeSwitch::LocalShuffleAtPartialBin() {
	
	LI randnum1, randnum2;
	Vertex  u, v, u2, v2;
	
	// pick the first edge uniformly at random
	randnum1 = (LI)rand() % globBinSize;
	u = globBin[randnum1].u;
	v = globBin[randnum1].v;

	// pick the second edge uniformly at random
	randnum2 = (LI)rand() % globBinSize;
	u2 = globBin[randnum2].u;
	v2 = globBin[randnum2].v;

	// check for busy edge, loop, parallel edge and edge switch that does not change the edges
	if (isBusy[randnum1] || isBusy[randnum2] || u==u2 || u==v2 || v==u2 || v==v2 || CreateParallelEdgeInGBin(randnum1,v2) || CreateParallelEdgeInGBin(randnum2,v) || EdgeIsPending(u,v2) || EdgeIsPending(u2,v))
		tag = 3;
	else
		Swap(globBin[randnum1].v, globBin[randnum2].v);	// swap v and v2	
}


void EdgeSwitch::GlobalShuffleAtPartialBin() {

	rn = (LI)rand() % globBinSize;	// pick the first edge uniformly at random
	send[0] = p1 = globBin[rn].u;
	send[1] = q1 = globBin[rn].v;

	if(isBusy[rn])	// (p1,q1) is busy, discard this edge switch
		tag = 3;
	else {	// send this edge to the processor <Pj> with tag = 1
		isBusy[rn] = true;	// make (p1,q1) busy
		MPI_Send(send, 2, MPI_TYPE, Pj, 1, MPI_COMM_WORLD);
		
		bool  temp  = false;
		while(temp == false) {
			MPI_Recv(recv, 2, MPI_TYPE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			temp = ServeRecvdMessages();
		}
	}
}


void EdgeSwitch::SendTerminateSignal() {
	for (int proc = globalBinMetaData.startProc; proc <= globalBinMetaData.endProc; proc++)
		if(proc != myRank)
			MPI_Send(send, 1, MPI_TYPE, proc, 128, MPI_COMM_WORLD);
}


bool EdgeSwitch::ServeRecvdMessages() {
	
	DT  x1, y1, x2, y2, k, randnum, index, source;
	
	if(status.MPI_TAG == 1) {
		x1 = recv[0];
		y1 = recv[1];
		source = status.MPI_SOURCE;
		
		for(k = 0; k < MAX_SERV_ATTEMPT; k++) {		
			randnum = rand() % globBinSize;	// select 2nd edge
			x2  = globBin[randnum].u;
			y2  = globBin[randnum].v;
			
			if(isBusy[randnum] || x1 == x2 || y1 == y2 || x1 == y2 || x2 == y1 || CreateParallelEdgeInGBin(randnum,y1) || EdgeIsPending(x2,y1))
				continue;
			
			isBusy[randnum] = true;
			pendEdge[x2%PEDGERSIZE].push_back((pEdge){x2,y1,randnum});
			break;
		}
		
		if(k == MAX_SERV_ATTEMPT) {
			cout << "No of edge at proc " << myRank << ": " << globBinSize << endl;
			cout << "Maximum attempt taken at processor " << myRank << " for one edge switch inside ServeRecvdMessages function" << endl << "Exiting ...";
			exit(1);	
		}
		
		send[0] = x2; send[1] = y2;
		MPI_Send(send, 2, MPI_TYPE, source, 2, MPI_COMM_WORLD);
	}
	
	else if(status.MPI_TAG == 2) {
		
		send[0] = p2 = recv[0];
		q2 = recv[1];
		send[1] = q1;
		
		tag = (CreateParallelEdgeInGBin(rn,q2) || EdgeIsPending(p1,q2))? 13 : 5;	// 13 -> discard, 5 -> success and acknowledgement
		
		if(tag == 5)
			globBin[rn].v = q2;	// replace (p1,q1) by (p1,q2)
		isBusy[rn] = false;
		MPI_Send(send, 2, MPI_TYPE, Pj, tag, MPI_COMM_WORLD);	// send confirmation/acknowledgement/discard message
		return true;
	}
	
	else if(status.MPI_TAG == 5) {	// acknowledgement
		x2 = recv[0];
		y1 = recv[1];
		index = RemoveFromPending(x2,y1);
		globBin[index].v = y1;	// replace (x2,y2) by (x2,y1)
		isBusy[index] = false;
	}
	
	else if(status.MPI_TAG == 13) {		// discard
		x2 = recv[0];
		y1 = recv[1];
		index = RemoveFromPending(x2,y1);
		isBusy[index] = false;
	}
	
	else if(status.MPI_TAG == 128)
		terminateCount++;
	
	return false;
}


// check for parallel edge (x,y) where <x> is the <u> value of the edge at <index> position of <globBin>, <y> is the <v> value of the potential new edge
bool EdgeSwitch::CreateParallelEdgeInGBin(DT index, DT y) {

	DT x = globBin[index].u, idx;
	
	for(idx = index-1; idx >= 0 && globBin[idx].u == x; idx--) {
		if(globBin[idx].v == y)
			return true;
	}
	
	for(idx = index+1; idx < globBinSize && globBin[idx].u == x; idx++) {
		if(globBin[idx].v == y)
			return true;
	}
	
	return false;
}


// check whether the edge (x2,y1) is a pending edge
bool EdgeSwitch::EdgeIsPending(DT x2, DT y1) {
	
	int i = x2%PEDGERSIZE;
	for(int j=0; j<pendEdge[i].size(); j++)
		if(x2 == pendEdge[i][j].x && y1 == pendEdge[i][j].y) 
			return true;
	return false;
}



DT EdgeSwitch::RemoveFromPending(DT x2, DT y1) {
	
	DT i = x2%PEDGERSIZE, ind = -1;
	for(DT j=0; j<pendEdge[i].size(); j++) {
		if(x2 == pendEdge[i][j].x && y1 == pendEdge[i][j].y) {
			ind = pendEdge[i][j].indx;
			pendEdge[i][j] = pendEdge[i].back();
			pendEdge[i].pop_back();
		}
	}
	return ind;
}



void EdgeSwitch::WriteGraph(const char *oname) { // 1st param: output file/graph name

	LI u, sizeu, i, j;
	char *str = new char[strlen(oname) + 10];
	stringstream ss;
	ss << myRank;
	const char *ch = ss.str().c_str();
	strcpy(str, ch);
	strcat(str,".");
	strcat(str, oname);
	
	ofstream ofp(str);
	if(!ofp.is_open()) {
		cout << "Can not open file: " << str << endl;
		exit(1);
	}
	
	if(noLocBins > 0) {
		for(i=0; i<noLocBins; i++)
			for(j=0; j<binMetaData[i].noEdge; j++)
				ofp << bin[i][j].u << "\t" << bin[i][j].v << endl;
	}
	else {
		for(i=0; i<globBinSize; i++)
			ofp << globBin[i].u << "\t" << globBin[i].v << endl;
	}
	
	ofp.close();
	FreeMem(str);
}
// =================== END: WriteGraph function ===============================



// ==================== PrintAvgMaxMinTime function ===================================
int EdgeSwitch::PrintAvgMaxMinTime(double t) {
	
	double time_taken;
	recordTime dIn, dOut;
	dIn.index = myRank;
	dIn.value = t;
						// MAX-TIME
	MPI_Allreduce(&dIn, &dOut, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
	if(myRank == 0)	printf("Max. = %.4lf sec. (%d), ", dOut.value, dOut.index);
	int p = dOut.index;
						// MIN-TIME
	MPI_Reduce(&dIn, &dOut, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
	if(myRank == 0)	printf("Min. = %.4lf sec. (%d), ", dOut.value, dOut.index);
	
						// AVG-TIME
	MPI_Reduce(&dIn.value, &time_taken, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(myRank == 0)	printf("Avg. = %.4lf sec.\n\n", time_taken/noProcs);
	
	return p;
}
// =================== END: PrintAvgMaxMinTime function ===============================


// ======================= FreeGraph function =================================
void EdgeSwitch::FreeGraph() {
	
	FreeMem(cumProb);
	FreeMem(binMetaData);
	FreeAll(bin, noLocBins);
	vector<Edge>().swap(globBin);
	vector<bool>().swap(isBusy);
}
// ====================== END: FreeGraph function =============================


#endif /* #ifndef LABEL_SHUFFLE_H_ */