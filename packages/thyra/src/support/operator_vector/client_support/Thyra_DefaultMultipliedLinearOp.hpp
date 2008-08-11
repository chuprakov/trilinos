// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_MULTIPLICATIVE_LINEAR_OP_HPP
#define THYRA_DEFAULT_MULTIPLICATIVE_LINEAR_OP_HPP

#include "Thyra_DefaultMultipliedLinearOpDecl.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_Utils.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template<class Scalar>
DefaultMultipliedLinearOp<Scalar>::DefaultMultipliedLinearOp()
{}


template<class Scalar>
DefaultMultipliedLinearOp<Scalar>::DefaultMultipliedLinearOp(
  const int                                                   numOps
  ,const Teuchos::RCP<LinearOpBase<Scalar> >          Ops[]
  )
{
  initialize(numOps,Ops);
}


template<class Scalar>
DefaultMultipliedLinearOp<Scalar>::DefaultMultipliedLinearOp(
  const int                                                   numOps
  ,const Teuchos::RCP<const LinearOpBase<Scalar> >    Ops[]
  )
{
  initialize(numOps,Ops);
}


template<class Scalar>
void DefaultMultipliedLinearOp<Scalar>::initialize(
  const int                                                   numOps
  ,const Teuchos::RCP<LinearOpBase<Scalar> >          Ops[]
  )
{
  Ops_.resize(numOps);
  for( int k = 0; k < numOps; ++k )
    Ops_[k].initialize(Ops[k]);
  validateOps();
  setupDefaultObjectLabel();
}


template<class Scalar>
void DefaultMultipliedLinearOp<Scalar>::initialize(
  const int                                                   numOps
  ,const Teuchos::RCP<const LinearOpBase<Scalar> >    Ops[]
  )
{
  Ops_.resize(numOps);
  for( int k = 0; k < numOps; ++k )
    Ops_[k].initialize(Ops[k]);
  validateOps();
  setupDefaultObjectLabel();
}


template<class Scalar>
void DefaultMultipliedLinearOp<Scalar>::uninitialize()
{
  Ops_.resize(0);
  setupDefaultObjectLabel();
}


// Overridden form MultipliedLinearOpBase


template<class Scalar>
int DefaultMultipliedLinearOp<Scalar>::numOps() const
{
  return Ops_.size();
}


template<class Scalar>
bool DefaultMultipliedLinearOp<Scalar>::opIsConst(const int k) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( 0 <= k && k < numOps() ) );
#endif
  return Ops_[k].isConst();
}


template<class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
DefaultMultipliedLinearOp<Scalar>::getNonconstOp(const int k)
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( 0 <= k && k < numOps() ) );
#endif
  return Ops_[k].getNonconstObj();
}


template<class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultMultipliedLinearOp<Scalar>::getOp(const int k) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( 0 <= k && k < numOps() ) );
#endif
  return Ops_[k].getConstObj();
}


// Overridden from LinearOpBase


template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultMultipliedLinearOp<Scalar>::range() const
{
  assertInitialized();
  return getOp(0)->range();
}


template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultMultipliedLinearOp<Scalar>::domain() const
{
  assertInitialized();
  return getOp(numOps()-1)->domain();
}


template<class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultMultipliedLinearOp<Scalar>::clone() const
{
  return Teuchos::null; // Not supported yet but could be!
}


// Overridden from Teuchos::Describable

                                                
template<class Scalar>
std::string DefaultMultipliedLinearOp<Scalar>::description() const
{
  assertInitialized();
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  std::ostringstream oss;
  oss << Teuchos::Describable::description() << "{numOps = "<<numOps()<<"}";
  return oss.str();
}

template<class Scalar>
void DefaultMultipliedLinearOp<Scalar>::describe(
  Teuchos::FancyOStream                &out_arg
  ,const Teuchos::EVerbosityLevel      verbLevel
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  using Teuchos::RCP;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  assertInitialized();
  RCP<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  const int numOps = Ops_.size();
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      *out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      *out
        << Teuchos::Describable::description() << "{"
        << "rangeDim=" << this->range()->dim()
        << ",domainDim="<< this->domain()->dim() << "}\n";
      OSTab tab(out);
      *out
        <<  "numOps = "<< numOps << std::endl
        <<  "Constituent LinearOpBase objects for M = Op[0]*...*Op[numOps-1]:\n";
      OSTab tab2(out);
      for( int k = 0; k < numOps; ++k ) {
        *out << "Op["<<k<<"] = " << Teuchos::describe(*getOp(k),verbLevel);
      }
      break;
    }
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
}


// protected


// Overridden from SingleScalarLinearOpBase


template<class Scalar>
bool DefaultMultipliedLinearOp<Scalar>::opSupported(EOpTransp M_trans) const
{
  bool opSupported = true;
  for( int k = 0; k < static_cast<int>(Ops_.size()); ++k )
    if(!Thyra::opSupported(*getOp(k),M_trans)) opSupported = false;
  return opSupported;
  // ToDo: Cache these?
}


template<class Scalar>
void DefaultMultipliedLinearOp<Scalar>::apply(
  const EOpTransp                     M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "DefaultMultipliedLinearOp<Scalar>::apply(...)",*this,M_trans,X,Y
    );
#endif // TEUCHOS_DEBUG  
  const int numOps = Ops_.size();
  const Index m = X.domain()->dim();
  if( real_trans(M_trans)==NOTRANS ) {
    //
    // Y = alpha * M * X + beta*Y
    // =>
    // Y = alpha * op(Op[0]) * op(Op[1]) * ... * op(Op[numOps-1]) * X + beta*Y
    //
    RCP<MultiVectorBase<Scalar> > T_kp1, T_k; // Temporary propagated between loops 
    for( int k = numOps-1; k >= 0; --k ) {
      RCP<MultiVectorBase<Scalar> >         Y_k;
      RCP<const MultiVectorBase<Scalar> >   X_k;
      if(k==0)        Y_k = rcp(Y,false);  else Y_k = T_k = createMembers(getOp(k)->range(),m);
      if(k==numOps-1) X_k = rcp(&X,false); else X_k = T_kp1;
      if( k > 0 )
        Thyra::apply(*getOp(k),M_trans,*X_k,&*Y_k);
      else
        Thyra::apply(*getOp(k),M_trans,*X_k,&*Y_k,alpha,beta);
      T_kp1 = T_k;
    }
  }
  else {
    //
    // Y = alpha * M' * X + beta*Y
    // =>
    // Y = alpha * Op[numOps-1]' * Op[1]' * ... * Op[0]' * X + beta * Y
    //
    RCP<MultiVectorBase<Scalar> > T_km1, T_k; // Temporary propagated between loops 
    for( int k = 0; k <= numOps-1; ++k ) {
      RCP<MultiVectorBase<Scalar> >         Y_k;
      RCP<const MultiVectorBase<Scalar> >   X_k;
      if(k==numOps-1)   Y_k = rcp(Y,false);  else Y_k = T_k = createMembers(getOp(k)->domain(),m);
      if(k==0)          X_k = rcp(&X,false); else X_k = T_km1;
      if( k < numOps-1 )
        Thyra::apply(*getOp(k),M_trans,*X_k,&*Y_k);
      else
        Thyra::apply(*getOp(k),M_trans,*X_k,&*Y_k,alpha,beta);
      T_km1 = T_k;
    }
  }
}


// private


template<class Scalar>
void DefaultMultipliedLinearOp<Scalar>::validateOps()
{
  using Teuchos::toString;
#ifdef TEUCHOS_DEBUG
  try {
    const int numOps = Ops_.size();
    for( int k = 0; k < numOps; ++k ) {
      TEST_FOR_EXCEPT( Ops_[k]().get() == NULL );
      if( k < numOps-1 ) {
        THYRA_ASSERT_LINEAR_OP_TIMES_LINEAR_OP_SPACES_NAMES(
          "DefaultMultipliedLinearOp<Scalar>::initialize(...)"
          ,*Ops_[k].getConstObj(),NOTRANS,("Ops["+toString(k)+"]")
          ,*Ops_[k+1].getConstObj(),NOTRANS,("Ops["+toString(k+1)+"]")
          );
      }
    }
  }
  catch(...) {
    uninitialize();
    throw;
  }
#endif
}


template<class Scalar>
void DefaultMultipliedLinearOp<Scalar>::setupDefaultObjectLabel()
{
  std::ostringstream label;
  const int numOps = Ops_.size();
  for( int k = 0; k < numOps; ++k ) {
    std::string Op_k_label = Ops_[k].getConstObj()->getObjectLabel();
    if (Op_k_label.length() == 0)
      Op_k_label = "ANYM";
    if (k > 0)
      label << "*";
    label << "("<<Op_k_label<<")";
  }
  this->setObjectLabel(label.str());
  validateOps();
}


}	// end namespace Thyra


template<class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstMultiply(
  const Teuchos::RCP<LinearOpBase<Scalar> > &A,
  const Teuchos::RCP<LinearOpBase<Scalar> > &B,
  const std::string &M_label
  )
{
  using Teuchos::arrayArg;
  using Teuchos::RCP;
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
    multOp = Teuchos::rcp(
    new DefaultMultipliedLinearOp<Scalar>(
      2
      ,arrayArg<RCP<LinearOpBase<Scalar> > >(A,B)()
      )
    );
  if(M_label.length())
    multOp->setObjectLabel(M_label);
  return multOp;
}


template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::multiply(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &A,
  const Teuchos::RCP<const LinearOpBase<Scalar> > &B,
  const std::string &M_label
  )
{
  using Teuchos::arrayArg;
  using Teuchos::RCP;
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
    multOp = Teuchos::rcp(
    new DefaultMultipliedLinearOp<Scalar>(
      2
      ,arrayArg<RCP<const LinearOpBase<Scalar> > >(A,B)()
      )
    );
  if(M_label.length())
    multOp->setObjectLabel(M_label);
  return multOp;
}


template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::multiply(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &A,
  const Teuchos::RCP<const LinearOpBase<Scalar> > &B,
  const Teuchos::RCP<const LinearOpBase<Scalar> > &C,
  const std::string &M_label
  )
{
  using Teuchos::arrayArg;
  using Teuchos::RCP;
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
    multOp = Teuchos::rcp(
    new DefaultMultipliedLinearOp<Scalar>(
      3
      ,arrayArg<RCP<const LinearOpBase<Scalar> > >(A,B,C)()
      )
    );
  if(M_label.length())
    multOp->setObjectLabel(M_label);
  return multOp;
}


#endif	// THYRA_DEFAULT_MULTIPLICATIVE_LINEAR_OP_HPP
