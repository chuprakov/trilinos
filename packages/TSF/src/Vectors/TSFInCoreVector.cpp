#include "TSFInCoreVector.h"

using namespace TSF;

TSFInCoreVector::TSFInCoreVector(const TSFVectorSpace& space)
	: TSFAccessibleVector(space), hasMoreChunks_(true)
{}

const TSFReal* TSFInCoreVector::getNextChunk(int& chunkSize) const
{
	chunkSize = nLocal() + nGhost();
	hasMoreChunks_ = false;
	return dataPointer();
}

TSFReal* TSFInCoreVector::getNextChunk(int& chunkSize) 
{
	chunkSize = nLocal() + nGhost();
	hasMoreChunks_ = false;
	return dataPointer();
}

const TSFReal* TSFInCoreVector::getNextLocalChunk(int& chunkSize) const
{
	chunkSize = nLocal() ;
	hasMoreChunks_ = false;
	return dataPointer();
}

