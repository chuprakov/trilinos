#include "KayLoghinSchurFactory.h"

#include "SchurFactory.h"

#include "SchurFactoryBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFLinearOperatorBase.h"
#include "TSFLinearSolverBase.h"
#include "TSFIdentityOperator.h"
#include "TSFOperatorSourceBase.h"
#include "KayLoghinRightOperatorSource.h"

using namespace SPP;
using namespace TSF;
using std::string;

KayLoghinSchurFactory::KayLoghinSchurFactory(TSFLinearSolver& ApSolver)
  : ApSolver_(ApSolver)
{}



TSFLinearOperator KayLoghinSchurFactory
::getSchurInvApprox(const TSFOperatorSource& OpSrc) const
{
  TSFLinearOperator Xinv; // get the right domain and range space for Xinv?
 
  // Cast OpSrc to KayLoghinRightOperatorSource type.
  const KayLoghinRightOperatorSource* klSrcPtr = 
    dynamic_cast<const KayLoghinRightOperatorSource*>(OpSrc.ptr());

  // Add check that operator source is really the right type.

  // Get saddle operator if we need it:
  //  TSFLinearOperator S = klSrcPtr->getOp();
  //  cout << "Describe S:" << endl;
  //  S.describe();

  // Get Ap operator and set up Apinv using ApSolver
  TSFLinearOperator Ap = klSrcPtr->getAp();
  TSFLinearOperator Apinv = Ap.inverse(ApSolver_);

  // Get Fp operator
  TSFLinearOperator Fp = klSrcPtr->getFp();

  // Using Mp = I for now
  TSFLinearOperator Mpinv = new TSFIdentityOperator(Fp.domain());

  Xinv = -Mpinv * Fp * Apinv;
  return Xinv;
}

string KayLoghinSchurFactory::toString() const
{
	return "KayLoghinSchurFactory toString";
}
