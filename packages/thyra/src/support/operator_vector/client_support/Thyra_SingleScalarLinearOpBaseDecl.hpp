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

#ifndef THYRA_SINGLE_SCALAR_LINEAR_OP_BASE_DECL_HPP
#define THYRA_SINGLE_SCALAR_LINEAR_OP_BASE_DECL_HPP

#include "Thyra_LinearOpDefaultBaseDecl.hpp"

namespace Thyra {

/** \brief Base class for linear operators that can only handle a single
 * scalar type.
 *
 * This class is meant to provide an easier way for subclasses to provide
 * implementations for <tt>LinearOpBase::apply()</tt> and
 * <tt>LinearOpBase::applyTranspose()</tt> and is not meant to be used as an
 * client interface.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class SingleScalarLinearOpBase : virtual public LinearOpDefaultBase<Scalar> {
public:

  using LinearOpDefaultBase<Scalar>::apply;

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief . */
  bool applySupports( const EConj conj ) const;
  /** \brief . */
  bool applyTransposeSupports( const EConj conj ) const;
  /** \brief . */
  void apply(
    const EConj                       conj
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;
  /** \brief . */
  void applyTranspose(
    const EConj                       conj
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;

  //@}

protected:

  /** @name Pure virtual functions (must be overridden by subclass) */
  //@{

  /** \brief Return if the <tt>M_trans</tt> operation of <tt>apply()</tt> is
   * supported or not.
   *
   * Preconditions:<ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt>
   * </ul>
   *
   * Note that an operator must support at least one of the values
   * of <tt>ETrans</tt> (i.e. the transposed or the non-transposed
   * operations must be supported, both can not be unsupported)
   */
  virtual bool opSupported(EOpTransp M_trans) const = 0;

  /** \brief Apply the linear operator (or its transpose).
   *
   * ToDo: Finish documentation!
   *
   * Preconditions:<ul>
   * <li><tt>this->opSupported(M_trans)==true</tt>
   * </ul>
   */
  virtual void apply(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    MultiVectorBase<Scalar> *Y,
    const Scalar alpha,
    const Scalar beta
    ) const = 0;

  //@}

};	// end class SingleScalarLinearOpBase

}	// end namespace Thyra

#endif	// THYRA_SINGLE_SCALAR_LINEAR_OP_BASE_DECL_HPP
