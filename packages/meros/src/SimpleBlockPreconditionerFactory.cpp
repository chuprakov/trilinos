#include "SimpleBlockPreconditionerFactory.h"
#include "TSFPreconditionerFactory.h"
#include "GenericRightPreconditioner.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFMatrixOperator.h"
#include "TSFIdentityOperator.h"
#include "TSFBlockLinearOperator.h"
#include "TSFUtils.h"
#include "RightBlockNSOperatorSource.h"
#include "SimpleOperatorSource.h"
#include "BlockForwardsolver.h"

using namespace SPP;
using namespace TSF;
using std::string;
//Need to update to reflect the Schur solver (pick a default solver)
SimpleBlockPreconditionerFactory
::SimpleBlockPreconditionerFactory(const TSFLinearSolver& Fsolver, const SchurFactory& sfac)
 	: TSFPreconditionerFactoryBase(), Fsolver_(Fsolver), sfac_(sfac)
{}

SimpleBlockPreconditionerFactory
::SimpleBlockPreconditionerFactory(const TSFLinearSolver& Fsolver, const SchurFactory& sfac, const TSFLinearSolver& Schursolver)
  : TSFPreconditionerFactoryBase(), Fsolver_(Fsolver), sfac_(sfac), Schursolver_(Schursolver)
{}

TSFPreconditioner SimpleBlockPreconditionerFactory
::createPreconditioner(const TSFLinearOperator& op) const
{
  TSFError::raise("SimpleBlockPreconditionerFactory::createPreconditioner called "
                  "with TSFLinearOperator instead of TSFOperatorSource");
  return 0;
}

TSFPreconditioner SimpleBlockPreconditionerFactory
::createPreconditioner(const TSFOperatorSource& opSource) const
{
  // SIMPLE style block preconditioner
  // 
  //      | I   -inv(diag(F))Bt | |   F        0               |-1
  //      | 0        I          | |   B     -B inv(diag(F)) Bt |
  //

  // Get the saddle point operator S
  TSFLinearOperator S = opSource.getOp();

   const SimpleOperatorSource* diagSrcPtr = dynamic_cast<const SimpleOperatorSource*>(opSource.ptr());

  // Get the blocks from S
  TSFLinearOperator F = S.getBlock(0,0);
  TSFLinearOperator Finv = diagSrcPtr->getDinv(); 
  
  TSFLinearOperator Bt = S.getBlock(0,1);
  TSFLinearOperator B = S.getBlock(1,0);
  //	TSFLinearOperator C = S.getBlock(1,1);
  
  TSFLinearOperator FinvBt = Finv*Bt;
  	
  // Build the Schur complement approx using the SchurFactory
  TSFLinearOperator X = B*Finv*Bt;

  // Make identity matrices on the right spaces
  TSFLinearOperator Iv = new TSFIdentityOperator(F.domain());
  TSFLinearOperator Ip = new TSFIdentityOperator(Bt.domain());

  // Make block operators for the three preconditioner factors
  TSFLinearOperator upper = new TSFBlockLinearOperator(S.domain(), S.range());
  TSFLinearOperator lower = new TSFBlockLinearOperator(S.domain(), S.range());

  // Put the pieces together into the composed block preconditioner
   lower.setBlock(0,0,F);
   lower.setBlock(1,0,B);
   lower.setBlock(1,1,X);

   upper.setBlock(0,0,Iv);
   upper.setBlock(0,1,FinvBt);
   upper.setBlock(1,1,Ip);	

   //Set up the TSF lower triangular solver with the F and Schur complement solvers
   TSFLinearSolver sol = new BlockForwardsolver(Fsolver_,Schursolver_); 
   TSFLinearOperator Linv = lower.inverse(sol);
   TSFLinearOperator prec = upper*Linv;

  return new GenericRightPreconditioner(prec);
}

string SimpleBlockPreconditionerFactory::toString() const 
{
	return "NS block preco factory";
}





