#include "BJBSchurFactory.h"

#include "SchurFactory.h"

#include "SchurFactoryBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFLinearOperatorBase.h"
#include "TSFLinearSolverBase.h"
#include "TSFIdentityOperator.h"
#include "TSFOperatorSourceBase.h"
#include "BJBRightOperatorSource.h"

using namespace SPP;
using namespace TSF;
using std::string;

BJBSchurFactory::BJBSchurFactory(TSFLinearSolver& ApSolver)
  : ApSolver_(ApSolver)
{}



TSFLinearOperator BJBSchurFactory
::getSchurInvApprox(const TSFOperatorSource& OpSrc) const
{
  TSFLinearOperator Xinv; // get the right domain and range space for Xinv?
 
  // Cast OpSrc to BJBRightOperatorSource type.
  const BJBRightOperatorSource* bjbSrcPtr = 
    dynamic_cast<const BJBRightOperatorSource*>(OpSrc.ptr());

  // Add check that operator source is really the right type.

  // Get saddle operator if we need it:
  TSFLinearOperator S = bjbSrcPtr->getOp();

  TSFLinearOperator F = S.getBlock(0,0);
  TSFLinearOperator Bt = S.getBlock(0,1);
  TSFLinearOperator B = S.getBlock(1,0);

  //  cout << "Describe S:" << endl;
  //  S.describe();

  // Get Ap operator and set up Apinv using ApSolver
  TSFLinearOperator Ap = bjbSrcPtr->getAp();
  TSFLinearOperator Apinv = Ap.inverse(ApSolver_);

  Xinv = -Apinv * B * F * Bt * Apinv;
  return Xinv;
}

string BJBSchurFactory::toString() const
{
	return "BJBSchurFactory toString";
}
