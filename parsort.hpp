/* merge sort */
#include <cstdio>
#include <iostream>
#include <mpi.h>
#include <ctime>
#include <cstdlib>

using namespace std;
typedef  long int PLI;

void showElapsed(int id, char *m);
void showVector(int *v, int n, int id);
void merge(PLI *A, int *AI, int asize, PLI *B, int *BI, int bsize, PLI*, int*);
void m_sort(PLI *A, int *I, int min, int max);

double startT, stopT;
double startTime;

void showElapsed(int id, char *m)
{
	printf("%d: %s %f secs\n",id,m,(clock()-startTime)/CLOCKS_PER_SEC);
}

void showVector(int *v, int n, int id)
{
	int i;
	printf("%d: ",id);
	for(i=0;i<n;i++)
		printf("%d ",v[i]);
	putchar('\n');
}

void merge(PLI *A, int *IA, int asize, PLI *B, int *IB, int bsize, PLI *C, int *IC) {
	int ai, bi, ci, i;
	int csize = asize + bsize;
	
	ai = 0;
	bi = 0;
	ci = 0;

	/* printf("asize=%d bsize=%d\n", asize, bsize); */

	while ((ai < asize) && (bi < bsize)) {
		if (A[ai] >= B[bi]) {
			C[ci]  = A[ai];
			IC[ci] = IA[ai];
			ci++; ai++;
		} else {
			C[ci]  = B[bi];
			IC[ci] = IB[bi];
			ci++; bi++;
		}
	}

	if (ai >= asize)
		for (i = ci; i < csize; i++, bi++) {
			C[i]  = B[bi];
			IC[i] = IB[bi];
		}
	else if (bi >= bsize)
		for (i = ci; i < csize; i++, ai++) {
			C[i]  = A[ai];
			IC[i] = IA[ai];
		}

	for (i = 0; i < asize; i++) {
		A[i]  = C[i];
		IA[i] = IC[i];
	}
	for (i = 0; i < bsize; i++) {
		B[i]  = C[asize+i];
		IB[i] = IC[asize+i];
	}
}

void m_sort(PLI *A, int *IA, int min, int max)
{
	int mid = (min+max)/2;
	int lowerCount = mid - min + 1;
	int upperCount = max - mid;
	PLI *C; int *IC;	// mergedData pointer
	/* If the range consists of a single element, it's already sorted */
	if (max == min) {
		return;
	} else {
		/* Otherwise, sort the first half */
		m_sort(A, IA, min, mid);
		/* Now sort the second half */
		m_sort(A, IA, mid+1, max);
		/* Now merge the two halves */
		C  = new  PLI[lowerCount+upperCount];
		IC = new int[lowerCount+upperCount];
		merge(A + min, IA + min, lowerCount, A + mid + 1, IA + mid + 1, upperCount, C, IC);
	}
}

// parallel merge sort in descending order ...
void ParMergeSortD(PLI *&data, int *&index, int dataSize, int rank, int noProcs)
{
	PLI  *dataChunk, *mergedData;
	PLI  *dataOther;
	int *indexChunk, *mergedIndex;
	int *indexOther;
	int m, n;
	int id, p;
	int s;
	int i, j, k;
	int step;
	MPI_Status status;

	//startT = clock();
	
	n  = dataSize;
	id = rank;
	p  = noProcs;
	s  = n/p;
	
	if(n%p) {
		cout << "\n\nError! Datasize must be a multiple of number of processors\n\n";
		exit(1);
	}
	
	dataChunk  = new PLI[s];
	indexChunk = new int[s];
	
	for(i=0, j=id*s; i<s; i++, j++) {
		dataChunk[i]  = data[j];
		indexChunk[i] = index[j];
	}
	
	m_sort(dataChunk, indexChunk, 0, s-1);

	step = 1;
	while(step < p)
	{
		if(id%(2*step)==0)
		{
			if(id+step<p)
			{
				MPI_Recv(&m, 1, MPI_INT, id+step, 0, MPI_COMM_WORLD, &status);
				dataOther  = new PLI[m];
				indexOther = new int[m];
				MPI_Recv(dataOther, m, MPI_LONG, id+step, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(indexOther, m, MPI_INT, id+step, 0, MPI_COMM_WORLD, &status);
				mergedData  = new  PLI[s+m];
				mergedIndex = new int[s+m];
				merge(dataChunk, indexChunk, s, dataOther, indexOther, m, mergedData, mergedIndex);
				dataChunk  = mergedData;
				indexChunk = mergedIndex;
				s += m;
			}
		}
		else
		{
			int near = id-step;
			MPI_Send(&s, 1, MPI_INT, near, 0, MPI_COMM_WORLD);
			MPI_Send(dataChunk, s, MPI_LONG, near, 0, MPI_COMM_WORLD);
			MPI_Send(indexChunk, s, MPI_INT, near, 0, MPI_COMM_WORLD);
			break;
		}
		step = step*2;
	}
	
	if(rank == 0) {
		data  = dataChunk;
		index = indexChunk;
	}

	/*
	stopT = clock();
	if(id==0)
	{
		//printf("%d; %d processors; %f secs\n",s,p,(stopT-startT)/CLOCKS_PER_SEC);

		FILE *fout = fopen("result","w");
		for(i=0;i<s;i++) {
			fprintf(fout, "i = %d, data = %ld, index = %d\n", i, dataChunk[i], indexChunk[i]);
			if(i < s-1 && dataChunk[i] < dataChunk[i+1])
				fprintf(fout, "\n\nERROR\n\n");
		}
		fclose(fout);
	}
	*/
}
