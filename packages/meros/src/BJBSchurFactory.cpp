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
  TSFLinearOperator C = S.getBlock(1,1);

  //  cout << "Describe S:" << endl;
  //  S.describe();

 if (C.isZeroOperator())
   {
  // Get Ap operator and set up Apinv using ApSolver
     cerr << "\n C is zero!";
  TSFLinearOperator Ap = bjbSrcPtr->getAp();
  TSFLinearOperator Apinv = Ap.inverse(ApSolver_);
  Xinv = -Apinv * B * F * Bt * Apinv;
   }
 else
   {
     TSFReal alpha, beta;
     bool noscal = bjbSrcPtr->getScal();
     if(noscal)
       {
   alpha = bjbSrcPtr->getAlpha();
   beta = bjbSrcPtr->getBeta();
   cerr << "\nC is not zero! creating preconditioner with beta= " << beta << " and alpha= " << alpha << endl;
       }
     else
       {
       alpha = 1.0;
       beta = 1.0;  
       cerr << "\n Alpha and beta not given, calculating ... alpha is " << alpha << " and beta is " << beta << endl;
       }
  TSFLinearOperator Dinv = bjbSrcPtr->getDinv();
  TSFLinearOperator Cinv = C.inverse(ApSolver_);
  Xinv = -beta*Dinv + alpha*alpha*Cinv*(B*F*Bt)*Cinv;
   }
 Xinv.describe();
  return Xinv;
}

string BJBSchurFactory::toString() const
{
	return "BJBSchurFactory toString";
}
