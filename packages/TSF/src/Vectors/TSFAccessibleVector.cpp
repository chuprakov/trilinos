#include "TSFAccessibleVector.h"
#include "TSFError.h"
#include "TSFBlas.h"
#include "TSFRandomNumberGenerator.h"
#include "TSFVector.h"


#ifndef DBL_MAX
#define DBL_MAX 1.0e100
#endif

using namespace TSF;

TSFAccessibleVector::TSFAccessibleVector(const TSFVectorSpace& space)
	: TSFVectorBase(space) 
{;}


void TSFAccessibleVector::axpy(const TSFReal& a, const TSFVector& x)
{
	const TSFAccessibleVector* xPtr 
		= dynamic_cast<const TSFAccessibleVector*>(x.ptr());

	if (xPtr==0) TSFError::raise("cast error in axpy()");
	
	rewindChunkIterator();
	xPtr->rewindChunkIterator();
	int myChunkSize;
	int otherChunkSize;
	
	while(hasMoreChunks())
		{
			TSFReal* myChunk = getNextChunk(myChunkSize);
			const TSFReal* otherChunk = xPtr->getNextChunk(otherChunkSize);
			if (myChunkSize != otherChunkSize)
				{
					TSFError::raise("mismatched chunk sizes in doBinaryOperation");
				}
			TSFReal* yy = myChunk;
			TSFReal* xx = const_cast<TSFReal*>(otherChunk);
			int one=1;
			TSFBlas<TSFReal>::axpy(&myChunkSize, &a, xx, &one, yy, &one);
			restoreChunk(myChunk, myChunkSize);
		}
}

void TSFAccessibleVector::acceptCopyOf(const TSFVector& x)
{
	const TSFAccessibleVector* xPtr 
		= dynamic_cast<const TSFAccessibleVector*>(x.ptr());

	if (xPtr==0) TSFError::raise("cast error in acceptCopyOf()");
	
	rewindChunkIterator();
	xPtr->rewindChunkIterator();
	int myChunkSize;
	int otherChunkSize;
	
	while(hasMoreChunks())
		{
			TSFReal* myChunk = getNextChunk(myChunkSize);
			const TSFReal* otherChunk = xPtr->getNextChunk(otherChunkSize);
			if (myChunkSize != otherChunkSize)
				{
					TSFError::raise("mismatched chunk sizes in doBinaryOperation");
				}
			TSFReal* yy = myChunk;
			TSFReal* xx = const_cast<TSFReal*>(otherChunk);
			int one = 1;
			TSFBlas<TSFReal>::copy(&myChunkSize, xx, &one, yy, &one);
			restoreChunk(myChunk, myChunkSize);
		}
}

void TSFAccessibleVector::scalarMult(const TSFReal& a)
{
	rewindChunkIterator();
	int myChunkSize;
	
	while(hasMoreChunks())
		{
			TSFReal* myChunk = getNextChunk(myChunkSize);
			TSFReal* x = myChunk;
			int one = 1;
			TSFBlas<TSFReal>::scal(&myChunkSize, &a, x, &one);
			restoreChunk(myChunk, myChunkSize);
		}
}


TSFReal TSFAccessibleVector::dot(const TSFVector& other) const
{
	const TSFAccessibleVector* yPtr 
		= dynamic_cast<const TSFAccessibleVector*>(other.ptr());

	if (yPtr==0) TSFError::raise("cast error in axpy()");
	
	rewindChunkIterator();
	yPtr->rewindChunkIterator();
	int myChunkSize;
	int otherChunkSize;
	
	TSFReal sum = 0.0;
	while(hasMoreChunks())
		{
			const TSFReal* myChunk = getNextLocalChunk(myChunkSize);
			const TSFReal* otherChunk = yPtr->getNextLocalChunk(otherChunkSize);
			if (myChunkSize != otherChunkSize)
				{
					TSFError::raise("mismatched chunk sizes in doBinaryOperation");
				}
			TSFReal* x = const_cast<TSFReal*>(myChunk);
			TSFReal* y = const_cast<TSFReal*>(otherChunk);
			int one = 1;
			sum += TSFBlas<TSFReal>::dot(&myChunkSize, x, &one, y, &one);
		}
	return reduce(sum, TSFSumReduce);
}

TSFReal TSFAccessibleVector::norm1() const
{
	rewindChunkIterator();
	int myChunkSize;
	
	TSFReal sum = 0.0;
	while(hasMoreChunks())
		{
			const TSFReal* myChunk = getNextLocalChunk(myChunkSize);
			TSFReal* x = const_cast<TSFReal*>(myChunk);
			for (int i=0; i<myChunkSize; i++, x++) 
				{
					sum += fabs(*x);
				}
		}
	return reduce(sum, TSFSumReduce);
}

TSFReal TSFAccessibleVector::normInf() const
{
	rewindChunkIterator();
	int myChunkSize;
	
	TSFReal mx = DBL_MAX;
	while(hasMoreChunks())
		{
			const TSFReal* myChunk = getNextLocalChunk(myChunkSize);
			TSFReal* x = const_cast<TSFReal*>(myChunk);
			TSFReal f;
			for (int i=0; i<myChunkSize; i++, x++) 
				{
					if ((f=fabs(*x))>mx) mx = f;
				}
		}
	return reduce(mx, TSFMaxReduce);
}

TSFReal TSFAccessibleVector::sumElements() const
{
	rewindChunkIterator();
	int myChunkSize;
	
	TSFReal sum = 0.0;
	while(hasMoreChunks())
		{
			const TSFReal* myChunk = getNextLocalChunk(myChunkSize);
			TSFReal* x = const_cast<TSFReal*>(myChunk);
			for (int i=0; i<myChunkSize; i++, x++) 
				{
					sum += *x;
				}
		}
	return reduce(sum, TSFSumReduce);
}

void TSFAccessibleVector::setScalar(const TSFReal& a)
{
	rewindChunkIterator();
	int myChunkSize;
	
	while(hasMoreChunks())
		{
			TSFReal* myChunk = getNextChunk(myChunkSize);
			TSFReal* x = myChunk;
			for (int i=0; i<myChunkSize; i++, x++) 
				{
					*x = a;
				}
			restoreChunk(myChunk, myChunkSize);
		}
}

void TSFAccessibleVector::randomize(const TSFRandomNumberGenerator& r)
{
	rewindChunkIterator();
	int myChunkSize;
	
	while(hasMoreChunks())
		{
			TSFReal* myChunk = getNextChunk(myChunkSize);
			TSFReal* x = myChunk;
			r.generateRandomNumbers(myChunkSize, x);
			restoreChunk(myChunk, myChunkSize);
		}
	invalidateGhostValues();
	synchronizeGhostValues();
}

