/* 
 * File	  : Par-Lab-Shuffle.cpp 
 * Purpose: Shuffle labeled edges parallely (MPI program)
 * Author : Hasanuzzaman Bhuiyan
 * Email  : mhb@vbi.vt.edu
 * Date	  : November 14, 2014
 */

#include "Par-Lab-Shuffle-Com.hpp"
#include <sys/time.h>

//#define WRITE_AUX_OUTPUT_FILES

struct timerec {
    double 	value;
    int   	index;
};

using namespace std;

/*
	This program takes following 6 arguments/parameters as input
			1st arg : input graph file name
			2nd arg : demographic file name
			3rd arg : shuffle ratio
			4th arg : partitioning-factor
			5th arg : output graph file name
			6th arg : write the output file?
				 0 - No
				 1 - Yes, in edge list format
*/

int main(int argc, char **argv) {	
	int rank, size;
	double time_taken;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	EdgeSwitch g;
	
	#ifdef MANUAL_PARTITION_INPUT
		if (argc != 9) {
			if(rank == 0) {
				cout << "\nWrong command line argument input\n\nUsage: " << argv[0] << " <input-graph> <demographic-file> <shuffle-ratio> <partitioning-factor> <output-graph> <write-flag> <max-loc-bin-size> <no-global-procs-fraction>" << endl;
				cout << "Total parameters passed " << argc << endl;
			}
			MPI_Finalize();
		}
		g.Initialize(rank, size, atof(argv[3]), atof(argv[4]), atol(argv[7]), atof(argv[8]));	// INITIALIZE
	#else
		if (argc != 7) {
			if(rank == 0) {
				cout << "\nWrong command line argument input\n\nUsage: " << argv[0] << " <input-graph> <demographic-file> <shuffle-ratio> <partitioning-factor> <output-graph> <write-flag>" << endl;
				cout << "Total parameters passed " << argc << endl;
			}
			MPI_Finalize();
		}
		g.Initialize(rank, size, atof(argv[3]), atof(argv[4]));	// INITIALIZE
	#endif
	
	g.GraphPreProcessing(argv[1], argv[2]);	// READ DEMOGRAPHIC FILE IN PARALLEL

							// SHUFFLE
	srand((unsigned)time(NULL));
	MPI_Barrier(MPI_COMM_WORLD);	// used to determine precise shuffle time;  not necessary though
	g.Shuffle();
	
	
	#ifdef WRITE_AUX_OUTPUT_FILES
		char str[64];
		stringstream ss;
		ss << rank;
		const char *ch = ss.str().c_str();
		
		// PRINT - No. of local bin at this processor
		strcpy(str, ch);
		strcat(str, ".nobin");
		ofstream nbinFile(str);
		nbinFile << rank << "\t" << Max(g.noLocBins,1) << endl;
		nbinFile.close();
		
		// PRINT - No. of edge at this processor
		const char *ch1 = ss.str().c_str();
		strcpy(str, ch1);
		strcat(str, ".noedge");
		ofstream nedgeFile(str);
		nedgeFile << rank << "\t" << g.noLocEdges + g.globBinSize << endl;
		nedgeFile.close();
		
		// PRINT - No. of edge-switch performed by this processor
		const char *ch2 = ss.str().c_str();
		strcpy(str, ch2);
		strcat(str, ".noswitch");
		ofstream nswitchFile(str);
		nswitchFile << rank << "\t" << g.noLocSwitchP + g.myNoGlobShuffle << endl;
		nswitchFile.close();
		
		// PRINT - time required by this processor
		const char *ch3 = ss.str().c_str();
		strcpy(str, ch3);
		strcat(str, ".exectime");
		ofstream timeFile(str);
		timeFile << rank << "\t" << g.timeShuffle << endl;
		timeFile.close();
								// WRITE OUTPUT GRAPH
		if(atoi(argv[6]) == 1) {
			if(rank == 0)	cout << "\nWriting output graph ... ";
			g.WriteGraph(argv[5]);	// for each processor, output in edge list format ... append ".<rank>" at the end of file name for each processor
			if(rank == 0)	cout << "done.\n";
		}
	#endif
	
	
	if(rank == 0) cout << "\nExiting after successful job completion\n";
	
	MPI_Finalize();
	
	return 0;
}
