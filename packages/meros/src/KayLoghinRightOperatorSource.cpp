#include "KayLoghinRightOperatorSource.h"
#include "RightBlockNSOperatorSource.h"
#include "TSFOperatorSourceBase.h"

using namespace TSF;
using namespace SPP;
using std::string;

// KayLoghinRightOperatorSource::
// KayLoghinRightOperatorSource(TSFLinearOperator& S)
//	: S_(S), Fp_(), Ap_(), Mp_()
//{;}

KayLoghinRightOperatorSource::
KayLoghinRightOperatorSource(TSFLinearOperator& S,
                             TSFLinearOperator& Fp,
                             TSFLinearOperator& Ap)
 	: S_(S), Fp_(Fp), Ap_(Ap), Mp_()
{;}

// KayLoghinRightOperatorSource::
// KayLoghinRightOperatorSource(TSFLinearOperator& S,
//                              TSFLinearOperator& Fp,
//                              TSFLinearOperator& Ap,
//                              TSFLinearOperator& Mp)
// 	: S_(S), Fp_(Fp), Ap_(Ap), Mp_(Mp)
// {;}

// KayLoghinRightOperatorSource::
// KayLoghinRightOperatorSource(TSFLinearOperator& S,
//                              TSFLinearOperator& Ap)
// 	: S_(S), Fp_(), Ap_(Ap), Mp_()
// {;}
// // ...

TSFLinearOperator KayLoghinRightOperatorSource
::getOp() const
{
  return S_;
}

TSFLinearOperator KayLoghinRightOperatorSource
::getAp() const
{
  return Ap_;
}

TSFLinearOperator KayLoghinRightOperatorSource
::getFp() const
{
  return Fp_;
}
TSFLinearOperator KayLoghinRightOperatorSource
::getMp() const
{
  return Mp_;
}


string KayLoghinRightOperatorSource::toString() const 
{
	return "Kay Loghin operator source";
}
