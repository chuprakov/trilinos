#include "StokesRightPreconditionerFactory.h"
#include "ILUKPreconditionerFactory.h"
#include "TSFPreconditionerFactory.h"
#include "GenericRightPreconditioner.h"

#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFMatrixOperator.h"
#include "TSFIdentityOperator.h"
#include "TSFBlockLinearOperator.h"
#include "TSFUtils.h"

using namespace TSF;
using std::string;

StokesRightPreconditionerFactory
::StokesRightPreconditionerFactory(const TSFLinearSolver& innerSolver)
	: TSFPreconditionerFactoryBase(), innerSolver_(innerSolver)
{}

TSFPreconditioner StokesRightPreconditionerFactory
::createPreconditioner(const TSFLinearOperator& A) const
{
	TSFLinearOperator F = A.getBlock(0,0);
	//	TSFLinearOperator Finv = F.inverse(innerSolver_);

	TSFPreconditionerFactory pf = new ILUKPreconditionerFactory(1);
	TSFPreconditioner fp = pf.createPreconditioner(F);
	TSFLinearOperator Finv = fp.left();

	TSFLinearOperator Bt = A.getBlock(0,1);
	TSFLinearOperator C = A.getBlock(1,0);

	TSFLinearOperator I00 = new TSFIdentityOperator(F.domain());
	TSFLinearOperator I11= new TSFIdentityOperator(Bt.domain());

	TSFLinearOperator P1 = new TSFBlockLinearOperator(A.domain(), A.range());
	TSFLinearOperator P2 = new TSFBlockLinearOperator(A.domain(), A.range());
	TSFLinearOperator P3 = new TSFBlockLinearOperator(A.domain(), A.range());

	P1.setBlock(0, 0, Finv);
	P1.setBlock(1, 1, I11);

	P2.setBlock(0, 0, I00);
	P2.setBlock(0, 1, Bt);
	P2.setBlock(1, 1, -I11);

	P3.setBlock(0, 0, I00);
	P3.setBlock(1, 1, I11);

	return new GenericRightPreconditioner(P1*P2*P3);
	
}

string StokesRightPreconditionerFactory::toString() const 
{
	return "StokesRight factory";
}
