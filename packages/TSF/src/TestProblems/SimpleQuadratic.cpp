#include "SimpleQuadratic.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpace.h"
#include "TSFVector.h"

#include "TSFError.h"
#include "DenseSerialVectorSpace.h"
#include "LAPACKGeneralMatrix.h"
#include <math.h>

using namespace TSF;

SimpleQuadratic::SimpleQuadratic(double a, double b, double c)
	: TSFNonlinearOperatorBase(new DenseSerialVectorSpace(1),
														 new DenseSerialVectorSpace(1)),
	a_(a), b_(b), c_(c)
{}

void SimpleQuadratic::apply(const TSFVector& in, 
									TSFVector& out) const
{
	TSFReal x = in[0];
	
	out = range().createMember();

	out[0] = (a_*x + b_)*x + c_;
}

TSFLinearOperator SimpleQuadratic::derivative(const TSFVector& evalPt) const
{
	TSFReal x = evalPt[0];

	TSFReal dfdx = 2.0*a_*x + b_;

	LAPACKGeneralMatrix* J = new LAPACKGeneralMatrix(domain(), range());

	J->setElement(0, 0, dfdx);

	return J;
}
