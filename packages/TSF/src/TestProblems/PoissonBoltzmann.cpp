#include "PoissonBoltzmann.h"
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

PoissonBoltzmann::PoissonBoltzmann(int n)
	: TSFNonlinearOperatorBase(new DenseSerialVectorSpace(n),
														 new DenseSerialVectorSpace(n)),
	n_(n), uRight_(log(cosh(1.0))), h_(1.0/((TSFReal) n_-1))
{;}

void PoissonBoltzmann::apply(const TSFVector& in, 
														 TSFVector& out) const
{
	out = range().createMember();
	out[0] = in[0];
	out[n_-1] = in[n_-1] - uRight_;

	for (int i=1; i<n_-1; i++)
		{
			TSFReal uMinus = in[i-1];
			TSFReal u = in[i];
			TSFReal uPlus = in[i+1];
			TSFReal boltzmannTerm = exp(-0.5*u);
			out[i] = (uMinus + uPlus - 2.0*u)/h_/h_ - boltzmannTerm;
		}
}

TSFLinearOperator PoissonBoltzmann::derivative(const TSFVector& evalPt) const  
{
	LAPACKGeneralMatrix* J = new LAPACKGeneralMatrix(domain(), range());

	J->zero();

	J->setElement(0, 0, 1.0);
	J->setElement(n_-1, n_-1, 1.0);

	for (int i=1; i<n_-1; i++)
		{
			J->setElement(i, i-1, 1.0/h_/h_);
			J->setElement(i, i+1, 1.0/h_/h_);
			TSFReal boltzmannTerm = exp(-0.5*evalPt[i]);
			J->setElement(i, i, -2.0/h_/h_ - boltzmannTerm);
		}
	return J;
}
