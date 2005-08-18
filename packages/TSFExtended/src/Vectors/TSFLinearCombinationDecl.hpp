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

#ifndef TSFLINEARCOMBINATIONDECL_HPP
#define TSFLINEARCOMBINATIONDECL_HPP

#include "TSFConfigDefs.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFLinearOperatorDecl.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace TSFExtendedOps
{
  /** 
   * 
   */
  enum LCSign {LCAdd = 1, LCSubtract = -1};

  template <class Scalar>
  class ConvertibleToVector
  {
  public:
    virtual ~ConvertibleToVector();
    /** */
    virtual TSFExtended::Vector<Scalar> eval() const = 0 ;

    /** */
    Scalar operator*(const ConvertibleToVector<Scalar>& other) const ;

    /** */
    Scalar norm1() const ;

    /** */
    Scalar norm2() const ;

    /** */
    Scalar normInf() const ;
  };

  /** 
   * Class LC1 holds the information required to perform
   * a 1-term linear combination, i.e., a scalar times a vector
   * or an operator times a vectir.
   */
  template <class Scalar>
  class LC1 : public ConvertibleToVector<Scalar>
  {
  public:
    /** */
    LC1(const TSFExtended::Vector<Scalar>& x) ;

    /** */
    LC1(const Scalar& alpha, const TSFExtended::Vector<Scalar>& x);

    /** */
    LC1(const Scalar& alpha,
        const TSFExtended::LinearOperator<Scalar>& op, 
        const TSFExtended::Vector<Scalar>& x);

    /** 
     * Evaluate the term into the argument vector, overwriting 
     * the previous value of the argument. */
    void evalInto(TSFExtended::Vector<Scalar>& result) const ;

    /** Add the term into the argument vector */
    void addInto(TSFExtended::Vector<Scalar>& result, 
                 LCSign sign = LCAdd) const ;

    /** Evaluate the term and return its value */
    virtual TSFExtended::Vector<Scalar> eval() const ;

    /** Determine whether this term contains the given vector */
    bool containsVector(const Thyra::VectorBase<Scalar>* vec) const ;

    /** */
    LC1<Scalar> operator*(const Scalar& beta) const ;
    
  private:
    Scalar alpha_;

    TSFExtended::LinearOperator<Scalar> op_;

    TSFExtended::Vector<Scalar> x_;
  };


  /**
   * Class LCN contains an n-term linear combination. 
   */
  template <class Scalar, class Node1, class Node2>
  class LCN : public ConvertibleToVector<Scalar>
  {
    public:
    /** */
    LCN(const Node1& x1, const Node2& x2, LCSign sign = LCAdd);

    /** */
    operator TSFExtended::Vector<Scalar>() const ;
      
    /** */
    void evalInto(TSFExtended::Vector<Scalar>& result) const ;

    /** */
    void addInto(TSFExtended::Vector<Scalar>& result, 
                 LCSign sign = LCAdd) const ;

    /** */
    virtual TSFExtended::Vector<Scalar> eval() const ;

    /** */
    bool containsVector(const Thyra::VectorBase<Scalar>* vec) const ;
    
    private:
    Node1 x1_;

    Node2 x2_;

    LCSign sign_;
  };

}

#endif  /* DOXYGEN_DEVELOPER_ONLY */


#endif
