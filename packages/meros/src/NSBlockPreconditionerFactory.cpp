#include "NSBlockPreconditionerFactory.h"
#include "TSFPreconditionerFactory.h"
#include "GenericRightPreconditioner.h"

#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFMatrixOperator.h"
#include "TSFIdentityOperator.h"
#include "TSFBlockLinearOperator.h"
#include "TSFUtils.h"

#include "RightBlockNSOperatorSource.h"
#include "KayLoghinRightOperatorSource.h"

using namespace SPP;
using namespace TSF;
using std::string;

NSBlockPreconditionerFactory
::NSBlockPreconditionerFactory(const TSFLinearSolver& Fsolver, const SchurFactory& sfac)
 	: TSFPreconditionerFactoryBase(), Fsolver_(Fsolver), sfac_(sfac)
{}

TSFPreconditioner NSBlockPreconditionerFactory
::createPreconditioner(const TSFLinearOperator& op) const
{
  TSFError::raise("NSBlockPreconditionerFactory::createPreconditioner called "
                  "with TSFLinearOperator instead of TSFOperatorSource");
  return 0;
}

TSFPreconditioner NSBlockPreconditionerFactory
::createPreconditioner(const TSFOperatorSource& opSource) const
{
  // Kay & Loghin style block preconditioner
  // 
  //      |  inv(F)   0   | |   I    -Bt  | |   I      0     |
  //      |    0      I   | |   0     I   | |   0   -inv(X)  |
  //

  // Get the saddle point operator S
  TSFLinearOperator S = opSource.getOp();
  
  // Get the blocks from S
  TSFLinearOperator F = S.getBlock(0,0);
  TSFLinearOperator Finv = F.inverse(Fsolver_);

  TSFLinearOperator Bt = S.getBlock(0,1);
  //	TSFLinearOperator C = S.getBlock(1,0);
  	
  // Build the Schur complement approx using the SchurFactory
  TSFLinearOperator Xinv = sfac_.getSchurInvApprox(opSource);

  // Make identity matrices on the right spaces
  TSFLinearOperator I00 = new TSFIdentityOperator(F.domain());
	TSFLinearOperator I11= new TSFIdentityOperator(Bt.domain());

  // Make block operators for the three preconditioner factors
	TSFLinearOperator P1 = new TSFBlockLinearOperator(S.domain(), S.range());
	TSFLinearOperator P2 = new TSFBlockLinearOperator(S.domain(), S.range());
	TSFLinearOperator P3 = new TSFBlockLinearOperator(S.domain(), S.range());

  // Put the pieces together into the composed block preconditioner

  P1.setBlock(0, 0, Finv);
	P1.setBlock(1, 1, I11);

	P2.setBlock(0, 0, I00);
	P2.setBlock(0, 1, -Bt);
	P2.setBlock(1, 1, I11);

	P3.setBlock(0, 0, I00);
	P3.setBlock(1, 1, -Xinv);

  return new GenericRightPreconditioner(P1*P2*P3);

}


string NSBlockPreconditionerFactory::toString() const 
{
	return "NS block preco factory";
}





