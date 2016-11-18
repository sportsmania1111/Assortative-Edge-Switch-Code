/* 
 * Modified for parallel labeled edge switches
 * File  : MinHeap.hpp
 * Author: Md Hasanuzzaman Bhuiyan
 * Email : mhb@vbi.vt.edu
 * Created on February 13, 2012
 * 
 * This is the implementation of "Minimum Heap" specifically for the Graph type algorithms
 * Here, array heap's indices are from <0 ... n-1> (0 to n-1) whereas in the Cormen book, the indices are assumed from <1 ... n> (1 to n)
 * It is the user's responsibility to copy the data initially into the "heap" array before calling any function
 *
 */

#ifndef _MinHeap_HPP
#define	_MinHeap_HPP

#include "utility.hpp"
#include "utility2.hpp"

#ifndef NULL
	#define  NULL 0
#endif

// Parent(i) = floor((i-1)/2) ... assuming heap starts from index 0 and ends at index (n-1)
#define Parent(i) 	(i-1) >> 1
// Left(i) = 2*i + 1
#define Left(i) 	(i << 1) + 1
// Right(i) = 2*i + 2
#define Right(i) 	(i << 1) + 2

//typedef  unsigned long int  ULI;
typedef  long int  LI;



// ================================ Class Declaration =========================


// NType -> node type (i.e. vertex number) and DType -> data type (i.e. value or distance data type)

template <typename NType, typename DType>
class MinHeap
{
	public:		
		NType 	*heap;			// constructs the main tree ... holds the indices of vertices
		DType 	*value;			// value (distance) of the vertices ... data type can be float/double/integer ... all the comparisons are made on this value of vertices
		LI 		*index;			// index holds the mapping of the vertices to the heap ... that means which vertex is at which position in the heap (or binary tree)
		LI 		length;			// total number of vertices in the heap
		LI 		heapsize;		// indicates that the heap property is maintained from "0" to position "heapsize"
		
		MinHeap()
		{
			heap 	 = NULL;
			value 	 = NULL;
			index 	 = NULL;
			length 	 = 0;                            
			heapsize = -1; 		// IMPORTANT: initially it is set to -1
		}
             
		MinHeap(LI n)
		{
			Initialize(n);
		}
             
		~MinHeap()
		{
			FreeMem(heap);
			FreeMem(value);
			FreeMem(index);
		}
		
		NType &operator[](LI i) { // "[]" operator overloaded returning reference to the heap value at position "i" of the heap
			return heap[i];
		}
             
		void Initialize(LI n);
		void MinHeapify(LI i);
		NType FindMinimum();
		void MakeHeap();             				
		NType ExtractMin();
		void IncreaseKey(NType i, DType key);
		void DecreaseKey(NType i, DType key);
		void Insert(NType element, DType key);		
		void Delete(NType element);
};



// ================================ Initialize ================================



// "" initializes the min heap data structure 
// IMPORTANT: 
//		1. heap, value and index array are just created here ... 
//		2. they are not initialized with the data ...
//		3. it is user's responsibility to initialize them with data ...

template <typename NType, typename DType>
void MinHeap<NType, DType>::Initialize(LI n)
{
    length 	 = n;
    heapsize = 0;
    heap 	 = new NType[n];  
	value 	 = new DType[n];
	index 	 = new LI[n];
}



// ================================ MinHeapify ================================



// maintains the minimum heap property at the position "i" of the heap
// here, "i" is the index/position of the heap array where this function will be applied ... it is "NOT" the vertex number 
// it assumes that the subtree rooted at the left and right child are already maintaining the min heap property

template <typename NType, typename DType>
void MinHeap<NType, DType>::MinHeapify (LI i)
{
	LI l, r, smallest;
	while(i < heapsize)
	{
		l = Left(i);
		r = Right(i);

		if(l <= heapsize && value[heap[l]] < value[heap[i]])
			smallest = l;
		else
			smallest = i;

		if(r <= heapsize && value[heap[r]] < value[heap[smallest]])
			smallest = r;

		if(smallest == i)	// the min heap property is already maintained ... nothing to check more ...
			break;
					
		Swap(index[heap[i]], index[heap[smallest]]); // updating the mapping info
		Swap(heap[i], heap[smallest]);	// swapping the elements/vertices
		i = smallest;		
	}
}



// ================================ MakeHeap ==================================


// "MakeHeap" ensures that the heap is maintaining the MIN_HEAP property
// that is the parent's value is less than or equal to the both children's value through out the whole heap
// [equivalent to "BuildMinHeap" of the Cormen's book]

template <typename NType, typename DType>
void MinHeap<NType, DType>::MakeHeap()
{
	heapsize = length-1;

	for(LI i = heapsize >> 1; i >= 0; i--)
		MinHeapify(i);
}



// ================================ FindMinimum ===============================


// returns the "index" of the vertex (NOT the value/distance) containing the least value/distance in the heap ... does "not" extract it from heap

template <typename NType, typename DType>
NType MinHeap<NType, DType>::FindMinimum()
{
    return heap[0];    
}



// ================================ ExtractMin ================================


// Extracts the "vertex with minimum value" from the heap and returns the "index" of that vertex

template <typename NType, typename DType>
NType MinHeap<NType, DType>::ExtractMin()
{        
    if(heapsize < 0)
    {
		printf("Error : Heap Underflow\n");
		exit(1);
	}    
	NType min = heap[0];
	heap[0] = heap[heapsize];
	index[heap[0]] = 0;
	heapsize--;
	if(heapsize > 0)
               MinHeapify(0);	
	return min;
}



// ================================ IncreaseKey ===============================


// "IncreaseKey" increases the value of the <i> -th vertex to <key> ... <i> is the "index" or "vertex no." of that vertex
// IMPORTANT: Here, <i> is "NOT" the position of that vertex in the heap
// index[i] = position of the <i> -th vertex in the heap

template <typename NType, typename DType>
void MinHeap<NType, DType>::IncreaseKey(NType i, DType key)
{
	if(value[i] > key)
	{
		printf("Error : New key is smaller than current key\n");
		return;
	}	
	value[i] = key;		
	MinHeapify(index[i]); // MinHeapify takes the "position of heap" as parameter, where it will maintain the MinHeap-property
}



// ================================ DecreaseKey ===============================


// "DecreaseKey" decreases the value of the <i>-th vertex to <key> ... <i> is the "index" or "vertex no." of that vertex
// IMPORTANT: Here, <i> is "NOT" the position of that vertex in the heap
// index[i] = position of the "i"-th vertex in the heap

template <typename NType, typename DType>
void MinHeap<NType, DType>::DecreaseKey(NType i, DType key)
{
	if(value[i] < key)
	{
		printf("Error : New key is larger than current key\n");
		return;
	}
	
	value[i] = key;	
	
	while(index[i] > 0 && value[heap[Parent(index[i])]] > value[i])
	{
		LI p = heap[Parent(index[i])];
		LI pindx = index[i];				
		Swap(heap[index[i]], heap[Parent(index[i])]);	// swapping	
		index[i] = Parent(index[i]);	// updating the mapping info
		index[p] = pindx;		
	}
}



// ================================ Insert ====================================


// inserts a new vertex with "vertex number or index" <element> and "value" <key>
// IMPORTANT: 
//		1. it assumes that the <element> is between [0 .. (n-1)]
// 		2. if index[element] == maxvalue(element), then this <element> is already NOT in the heap ... so insert it
//		   else it is already in the heap ... Error !!!

template <typename NType, typename DType>
void MinHeap<NType, DType>::Insert(NType element, DType key)
{
	if(heapsize + 1 >= length) // checking whether the heap is already full, and is there any space to insert?
	{
		printf("Error: attempting to insert into an already fulfilled heap\n"); // heap is full ... no insertion possible
		return;
	}	
	if(index[element] != maxvalue(element))
	{
		printf("Error: attempting to insert an element that already exist in the heap\n");
		return;
	}
	heapsize++;
	heap[heapsize] = element; 
	index[element] = heapsize;
	value[element] = maxValue(key); // equivalent to putting infinity for this new element
	HeapDecreaseKey(element, key);	
}



// ================================ Delete ====================================


// deletes the vertex with "vertex number or index" <element>
// IMPORTANT: 
//		1. it assumes that the <element> is between [0 .. (n-1)]
// 		2. if index[element] < maxValue(element), it means that the element exist in the heap and mapped to position index[element] of the heap,
//		   else it does not exist

template <typename NType, typename DType>
void MinHeap<NType, DType>::Delete(NType element)
{
	if(heapsize < 0) // checking whether the heap is already empty?
	{
		printf("Error: attempting to delete from an empty heap\n"); // heap is empty ... no deletion possible
		return;
	}	
	if(index[element] == maxValue(element)) // if the element's mapping is set to infinity, it means that this element does not exist in the heap
	{
		printf("Error: attempting to delete an element that does not exist\n");
		return;
	}	
	heap[index[element]] = heap[heapsize];
	index[heapsize] = index[element];	
	heapsize--;
	MinHeapify(index[element]);
	index[element] = maxValue(element); // setting the element's mapping to infinity indicating that this element does not exist in the heap	
}


#endif /* _MinHeap_HPP */ 
