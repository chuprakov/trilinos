#include "BJBRightOperatorSource.h"
#include "RightBlockNSOperatorSource.h"
#include "TSFOperatorSourceBase.h"
#include "TSFBlockLinearOperator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "PetraMatrix.h"

using namespace TSF;
using namespace SPP;
using std::string;


BJBRightOperatorSource::
BJBRightOperatorSource(TSFLinearOperator& S)
 	: S_(S), Ap_(), hasAp_(false)
{;}

BJBRightOperatorSource::
BJBRightOperatorSource(TSFLinearOperator& S,
                       TSFLinearOperator& Ap)
 	: S_(S), Ap_(Ap), hasAp_(true)
{;}


TSFLinearOperator BJBRightOperatorSource
::getOp() const
{
  return S_;
}

TSFLinearOperator BJBRightOperatorSource
::getAp() const
{

  // if Ap_ already exists (passed in or already built), return it
  if (hasAp_)
    return Ap_;
  
  // if Ap_ doesn't exist yet, need to get it from S
  else
    {
      cerr << "BJBRightOperatorSource: in right part of if" << endl;
      
      // Using Ap = C
      Ap_ = S_.getBlock(1,1);
      cerr << "BJBRightOperatorSource: got here" << endl;
      
      hasAp_ = true;
      return(Ap_);
      
    }
}


string BJBRightOperatorSource::toString() const 
{
	return "BJB operator source";
}
