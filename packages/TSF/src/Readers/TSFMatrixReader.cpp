#include "TSFMatrixReader.h"


using namespace TSF;

TSFMatrixReader::TSFMatrixReader(TSFMatrixReaderBase* ptr)
	: ptr_(ptr) 
{}

TSFLinearOperator TSFMatrixReader::read(const TSFVectorType& vectorType) const
{
	return ptr_->read(vectorType);
}
