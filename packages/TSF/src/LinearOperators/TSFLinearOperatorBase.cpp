#include "TSFError.h"
#include "TSFLinearOperator.h"
#include "TSFLinearOperatorBase.h"
#include "TSFMatrixOperator.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"

#include "TSFInverseOperator.h"
#include "TSFInverseAdjointOperator.h"
#include "TSFAdjointOperator.h"
#include "TSFLinearSolver.h"

using namespace TSF;
using std::ostream;

TSFLinearOperatorBase::TSFLinearOperatorBase()
	: domain_(), range_()
{
	;
}

TSFLinearOperatorBase::TSFLinearOperatorBase(const TSFVectorSpace& domain,
																						 const TSFVectorSpace& range)
	: domain_(domain), range_(range)
{
	;
}

void TSFLinearOperatorBase::getBlock(int /* i */, int /* j */, 
																		 TSFLinearOperator& /* sub */) const 
{
	TSFError::raise("TSFLinearOperatorBase::getBlock called for non-block operator");
}

void TSFLinearOperatorBase::setBlock(int /* i */, int /* j */, 
																		 const TSFLinearOperator& /* sub */)
{
	TSFError::raise("TSFLinearOperatorBase::setBlock called for non-block operator");
}



void TSFLinearOperatorBase::applyInverse(const TSFVector& /* arg */,
																				 TSFVector& /* out */) const
{
	TSFError::raise("TSFLinearOperatorBase::applyInverse not implememented in base class");
}

void TSFLinearOperatorBase::applyInverse(const TSFLinearSolver& /* solver */,
																				 const TSFVector& arg,
																				 TSFVector& out) const
{
	applyInverse(arg, out);
}

void TSFLinearOperatorBase::applyAdjoint(const TSFVector& /* arg */,
																				 TSFVector& /* out */) const
{
	TSFError::raise("TSFLinearOperatorBase::applyAdjoint not implememented in base class");

}

void TSFLinearOperatorBase::applyInverseAdjoint(const TSFVector& /* arg */,
																								TSFVector& /* out */) const
{
	TSFError::raise("TSFLinearOperatorBase::applyInverseAdjoint "
									" not implememented in base class");
}

void TSFLinearOperatorBase::applyInverseAdjoint(const TSFLinearSolver& /*solver*/,
																								const TSFVector& arg,
																								TSFVector& out) const
{
	applyInverseAdjoint(arg, out);
}

void TSFLinearOperatorBase::getInverse(const TSFLinearOperator& self,
																			 TSFLinearOperator& inv) const 
{
	try
		{
			inv = new TSFInverseOperator(self);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFLinearOperatorBase::getInverse");
		}
}

void TSFLinearOperatorBase::getInverse(const TSFLinearSolver& solver,
																			 const TSFLinearOperator& self,
																			 TSFLinearOperator& inv) const 
{
	try
		{
			inv = new TSFInverseOperator(self, solver);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFLinearOperatorBase::getInverse");
		}
}

void TSFLinearOperatorBase::getAdjoint(const TSFLinearOperator& self,
																			 TSFLinearOperator& adj) const 
{
	try
		{
			adj = new TSFAdjointOperator(self);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFLinearOperatorBase::getAdjoint");
		}
}

void TSFLinearOperatorBase::getInverseAdjoint(const TSFLinearOperator& self,
																							TSFLinearOperator& invAdj) const 
{
	try
		{
			invAdj = new TSFInverseAdjointOperator(self);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFLinearOperatorBase::getInverseAdjoint");
		}
}

void TSFLinearOperatorBase::getInverseAdjoint(const TSFLinearSolver& solver,
																							const TSFLinearOperator& self,
																							TSFLinearOperator& invAdj) const 
{
	try
		{
			invAdj = new TSFInverseAdjointOperator(self, solver);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFLinearOperatorBase::getInverseAdjoint");
		}
}


const TSFSmartPtr<const TSFMatrixOperator>
TSFLinearOperatorBase::getMatrix() const
{
	return TSFSmartPtr<const TSFMatrixOperator>();
}

void TSFLinearOperatorBase::print(ostream& os) const 
{
	os << "TSFLinearOperator[]";
}
