#include "TSFMatrixView.h"

#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"

using namespace TSF;

TSFMatrixView::TSFMatrixView(const TSFLinearOperator& op)
	: matrix_(0), op_(op)
{
	matrix_ = dynamic_cast<TSFMatrixOperator*>(&(*(op_.getPtr())));
	if (matrix_==0) 
		{
			TSFError::raise("attempted to create TSFMatrixView from a non-matrix operator");
		}
	
}

