/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
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
// **********************************************************************/
/* @HEADER@ */

#ifndef TSFNONLINEAROPERATORBASE_HPP
#define TSFNONLINEAROPERATORBASE_HPP

#include "TSFConfigDefs.hpp"
#include "TSFHandleable.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "TSFVectorDecl.hpp"

namespace TSFExtended
{
  using TSFCore::Index;
  using namespace Teuchos;

  /** 
   * Base class for nonlinear operators
   */
  template <class Scalar>
  class NonlinearOperatorBase 
    : public Handleable<NonlinearOperatorBase<Scalar> >,
    public ObjectWithVerbosity<NonlinearOperatorBase<Scalar> >
    {
    public:
      /** Empty ctor, for contexts in which we don't know the
       * domain and range spaces at the beginning of construction time */
      NonlinearOperatorBase() 
        : domain_(), range_(), 
          jacobianIsValid_(false),
          residualIsValid_(false),
          currentEvalPt_(),
          currentFunctionValue_(),
          currentJ_()
      {;}

      /** Construct a nonlinear operator with a domain and range */
      NonlinearOperatorBase(const VectorSpace<Scalar>& domain,
                            const VectorSpace<Scalar>& range) 
        : domain_(domain.ptr()), range_(range.ptr()), 
          jacobianIsValid_(false),
          residualIsValid_(false),
          currentEvalPt_(),
          currentFunctionValue_(),
          currentJ_()
      {;}
                            
      /** Return the domain space */
      const RefCountPtr<const TSFCore::VectorSpace<Scalar> >& domain() const 
      {return domain_;}

      /** Return the range space */
      const RefCountPtr<const TSFCore::VectorSpace<Scalar> >& range() const 
      {return range_;}

      /** Set the evaluation point */
      void setEvalPt(const Vector<Scalar>& x)
      {
        if (verbosity() > VerbLow)
          {
            cerr << "NonlinearOperatorBase Setting new eval pt";
            if (verbosity() > VerbHigh)
              {
                cerr << " to " << endl ;
                x.print(cerr);
              }
            cerr << endl;
          }
        jacobianIsValid_ = false;
        residualIsValid_ = false;

//         TEST_FOR_EXCEPTION(!x.space().isCompatible(*domain()),
//                            runtime_error,
//                            "evaluation point " << x
//                            << " for nonlinear operator is not in the "
//                            "operator's domain space ");
        
        currentEvalPt_ = x.copy();
      }

      /** Get the current point at which the function is to be 
       * evaluated */
      const Vector<double>& currentEvalPt() const {return currentEvalPt_;}

      /** Return the Jacobian at the current evaluation point */
       LinearOperator<double> getJacobian() const 
      {
        if (verbosity() > VerbLow)
          {
            cerr << "NonlinearOperatorBase getting Jacobian" << endl;
          }
        if (!jacobianIsValid_)
          {
            if (verbosity() > VerbLow)
              {
                cerr << "...computing new J and F" << endl;
              }
            currentJ_ 
              = computeJacobianAndFunction(currentFunctionValue_);
            jacobianIsValid_ = true;
            residualIsValid_ = true;
          }
        else
          {
            if (verbosity() > VerbLow)
              {
                cerr << "...reusing valid J" << endl;
              }
          }
        if (verbosity() > VerbHigh)
          {
            cerr << "J is " << endl;
            currentJ_.print(cerr);
            cerr << endl;
          }
        return currentJ_;
      }

      

      /** Return the function value at the current evaluation point */
      Vector<double> getFunctionValue() const 
      {
        if (verbosity() > VerbLow)
          {
            cerr << "NonlinearOperatorBase getting function value" << endl;
          }
        if (!residualIsValid_)
          {
            if (verbosity() > VerbLow)
              {
                cerr << "...computing new F" << endl;
              }
            currentFunctionValue_ = computeFunctionValue();
            residualIsValid_ = true;
          }
        else
          {
            if (verbosity() > VerbLow)
              {
                cerr << "...reusing valid F" << endl;
              }
          }

        if (verbosity() > VerbHigh)
          {
            cerr << "F is " << endl;
            currentFunctionValue_.print(cerr);
            cerr << endl;
          }
        return currentFunctionValue_;
      }


      /** Return an initial guess appropriate to this problem */
      virtual Vector<double> getInitialGuess() const = 0 ;


    protected:

      /** Compute the Jacobian at the current eval point */
      virtual LinearOperator<Scalar> computeJacobianAndFunction(Vector<double>& functionValue) const = 0 ;

      /** Compute the function value at the current eval point */
      virtual Vector<Scalar> computeFunctionValue() const 
      {
        computeJacobianAndFunction(currentFunctionValue_);
        return currentFunctionValue_;
      }

      
      /** Set the domain and range. This is protected so that solver
       * developers don't try to change the spaces on the fly */
      void setDomainAndRange(const VectorSpace<Scalar>& domain,
                             const VectorSpace<Scalar>& range)
      {
        domain_ = domain.ptr();
        range_ = range.ptr();
      }


    private:
      /** */
      RefCountPtr<const TSFCore::VectorSpace<Scalar> > domain_;

      /** */
      RefCountPtr<const TSFCore::VectorSpace<Scalar> > range_;

      /** */
      mutable bool jacobianIsValid_;

      /** */
      mutable bool residualIsValid_;

      /** */
      Vector<double> currentEvalPt_;

      /** */
      mutable Vector<double> currentFunctionValue_;

      /** */
      mutable LinearOperator<double> currentJ_;
    };



 
}


#endif
