// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
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

// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreLinOp.hpp

#ifndef TSFCORE_LIN_OP_HPP
#define TSFCORE_LIN_OP_HPP

#include "TSFCoreLinearOp.hpp"

namespace TSFCore {

///
/** Base class for simple aggregation of a <tt>LinearOp</tt> and its
 * natural (logical) transpose argument (see <tt>LinOpPersisting</tt>,
 * <tt>LinOpNonPersisting</tt> for concrete types).
 *
 * This class is nothing more than a base clas for a silly aggregation
 * of a <tt>LinearOp</tt> object and its <tt>ETransp</tt> value that
 * defines the mathematical definition of the non-transposed operator.
 *
 * This simple type has all of the operations of a type with value
 * semantics in that can be default constructed, copy constructed, and
 * assigned.  These operations are supplied by the compiler
 * automatically since the compiler-supplied versions have the correct
 * behavior.
 */
template<class Scalar>
class LinOpBase {
public:

	/** @name Public initializers/accessors */
	//@{

	///
	/** Set the transpose argument
	 *
	 * @param  defaultTrans  [in] Default definition of transpose (default: <tt>NOTRANS</tt>).
	 *
	 * Postconditions:<ul>
	 * <li>this->defaultTrans() == defaultTrans</tt>
	 * </ul>
	 */
	void defaultTrans( const ETransp defaultTrans );

	///
	/** Return the default defintion of mathematical transpose.
	 */
	ETransp defaultTrans() const;

	//@}

	/** @name LinearOp wrappers */
	//@{

	///
	/** Return the range space of the logical linear operator.
	 *
	 * Simply returns: \code

   return ( this->defaultTrans()==NOTRANS ? this->op()->range() : this->op()->domain() );
   \endcode
	 */
	Teuchos::RefCountPtr<const VectorSpace<Scalar> > range() const;

	///
	/** Return the domain space of the logical linear operator.
	 *
	 * Simply returns: \code

   return ( this->defaultTrans()==NOTRANS ? this->op()->domain() : this->op()->range() );
   \endcode
	 */
	Teuchos::RefCountPtr<const VectorSpace<Scalar> > domain() const;

	///
	/** Return if the operation is supported on the logical linear operator.
	 *
	 * Simply returns: \code

   return this->op()->opSupported(trans_trans(this->defaultTrans(),M_trans));
   \endcode
	 */
	bool opSupported(ETransp M_trans) const;

	///
	/** Apply the logical linear operator (or its transpose) to a vector:
	 * <tt>y = alpha*op(M)*x + beta*y</tt>.
	 *
	 * @param  M_trans
	 *                [in] Determines whether the logical transposed or non-trnasposed
	 *                operator is applied as:
	 *                <ul>
	 *                <li> <tt>op(M) = M</tt>, for <tt>M_trans==NOTRANS</tt>
	 *                <li> <tt>op(M) = M'</tt>, for <tt>M_trans==TRANS</tt>
	 *                </ul>
	 * @param  x      [in] See <tt>LinearOp::apply()</tt>.
	 * @param  y      [in/out] See <tt>LinearOp::apply()</tt>.
	 * @param  alpha  [in] See <tt>LinearOp::apply()</tt>.
	 * @param  beta   [in] See <tt>LinearOp::apply()</tt>.
	 * 
	 * Preconditions:<ul>
	 * <li> See <tt>LinearOp::apply()</tt>.
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> See <tt>LinearOp::apply()</tt>.
	 * </ul>
	 *
	 * Simply calls: \code

   this->op()->apply(trans_trans(M_trans,this->defaultTrans()),x,y,alpha,beta)
   \endcode
	 */
	void apply(
		const ETransp            M_trans
		,const Vector<Scalar>    &x
		,Vector<Scalar>          *y
		,const Scalar            alpha = 1.0
		,const Scalar            beta  = 0.0
		) const;

	///
	/** Apply the linear operator (or its transpose) to a multi-vector :
	 * <tt>Y = alpha*op(M)*X + beta*Y</tt>.
	 *
	 * @param  M_trans
	 *                [in] Determines whether the logical transposed or non-trnasposed
	 *                operator is applied as:
	 *                <ul>
	 *                <li> <tt>op(M) = M</tt>, for <tt>M_trans==NOTRANS</tt>
	 *                <li> <tt>op(M) = M'</tt>, for <tt>M_trans==TRANS</tt>
	 *                </ul>
	 * @param  X      [in] See <tt>LinearOp::apply()</tt>.
	 * @param  Y      [in/out] See <tt>LinearOp::apply()</tt>.
	 * @param  alpha  [in] See <tt>LinearOp::apply()</tt>.
	 * @param  beta   [in] See <tt>LinearOp::apply()</tt>.
	 * 
	 * Preconditions:<ul>
	 * <li> See <tt>LinearOp::apply()</tt>.
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> See <tt>LinearOp::apply()</tt>.
	 * </ul>
	 *
	 * Simply calls: \code

   this->op()->apply(trans_trans(M_trans,this->defaultTrans()),X,Y,alpha,beta)
	 \endcode
	 */
	void apply(
		const ETransp                 M_trans
		,const MultiVector<Scalar>    &X
		,MultiVector<Scalar>          *Y
		,const Scalar                 alpha = 1.0
		,const Scalar                 beta  = 0.0
		) const;

	//@}

protected:

	/** @name Protected Constructors/initializers/accessors */
	//@{

	///
	LinOpBase();
	///
	LinOpBase(
		const Teuchos::RefCountPtr<const LinearOp<Scalar> >   &op
		,const ETransp                                        defaultTrans
		);

	//@}

	/** @name Protected accessors */
	//@{

	///
	Teuchos::RefCountPtr<const LinearOp<Scalar> > op() const;

	//@}

private:

	Teuchos::RefCountPtr<const LinearOp<Scalar> >  op_;
	ETransp                                        defaultTrans_;

};

///
/** Simple aggregation class of a <tt>LinearOp</tt> and its natural
 * (logical) transpose argument using a persisting relationship.
 *
 * Objects of this type should be passed to functions where a
 * non-persisting relationship (as defined in the RefCountPtr
 * beginner's guide) is used with the underlying <tt>LinearOp</tt>
 * object.
 *
 * Objects of this type should never, never be used as data members in
 * a class since this is by definition a persisting relationship.  For
 * persisting relationships, use a <tt>LinOpPersisting</tt> object
 * instead.
 */
template<class Scalar>
class LinOpPersisting : public LinOpBase<Scalar> {
public:

	/** @name Public constructors/accessors */
	//@{

	///
	/** Default construct.
	 *
	 * Postconditions:<ul>
	 * <li>this->op().get() == NULL</tt>
	 * <li>this->defaultTrans() == NOTRANS</tt>
	 * </ul>
	 */
	LinOpPersisting();

	///
	/** Construct with an operator and optionally its default transpose argument.
	 *
	 * @param  op            [in] Smart pointer to linear operator (persisting relationship).
	 * @param  defaultTrans  [in] Default definition of transpose (default: <tt>NOTRANS</tt>).
	 *
	 *
	 * Preconditions:<ul>
	 * <li>op.get() != NULL</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li>this->op().get() == op.get()</tt>
	 * <li>this->defaultTrans() == defaultTrans</tt>
	 * </ul>
	 */
	LinOpPersisting(
		const Teuchos::RefCountPtr<const LinearOp<Scalar> >   &op
		,const ETransp                                        defaultTrans = NOTRANS
		);

	///
	/** Return smart pointer to underlying linear operator.
	 */
	Teuchos::RefCountPtr<const LinearOp<Scalar> > op() const;

	//@}

};

///
/** Simple aggregation class of a <tt>LinearOp</tt> and its natural
 * (logical) transpose argument using a non-persisting relationship.
 *
 * Objects of this type should be passed to functions where a
 * persisting relationship (as defined in the RefCountPtr beginner's
 * guide) is being created with the underlying <tt>LinearOp</tt>
 * object.  Objects of this type also should be used a private data
 * members to aggregate a <tt>LinearOp</tt> object and its
 * mathematical defintion of the non-transposed operator.
 *
 * This class has value semantics!
 */
template<class Scalar>
class LinOpNonPersisting : public LinOpBase<Scalar> {
public:

	/** @name Public constructors/accessors */
	//@{

	///
	/** Default construct.
	 *
	 * Postconditions:<ul>
	 * <li>this->op() == &op</tt>
	 * <li>this->defaultTrans() == NOTRANS</tt>
	 * </ul>
	 */
	LinOpNonPersisting();

	///
	/** Construct with an operator and optionally its default transpose argument.
	 *
	 * @param  op            [in] Raw reference to linear operator (non-persisting relationship).
	 * @param  defaultTrans  [in] Default definition of transpose (default: <tt>NOTRANS</tt>).
	 *
	 * Postconditions:<ul>
	 * <li>this->op() == &op</tt>
	 * <li>this->defaultTrans() == defaultTrans</tt>
	 * </ul>
	 */
	LinOpNonPersisting(
		const LinearOp<Scalar>                                &op
		,const ETransp                                        defaultTrans = NOTRANS
		);

	///
	/** Construct for a <tt>LinOpPersisting</tt> object.
	 *
	 * @param  op  [in]
	 *
	 * Postconditions:<ul>
	 * <li>this->op() == &op.op().get()</tt>
	 * <li>this->defaultTrans() == op.defaultTrans()</tt>
	 * </ul>
	 *
	 * This constructor essentially allows an implicit conversion from a
	 * <tt>LinOpPersisting</tt> object to a <tt>LinOpNonPersisting</tt>.
	 * This seems reasonable since one piece of code that maintains a
	 * persisting relationship to a <tt>LinearOp</tt> object
	 * (i.e. through a <tt>LinOpPersisting</tt> object) should be able
	 * to give a non-persisting relationship to the <tt>LinearOp</tt>
	 * object (i.e. through a <tt>LinOpNonPersisting</tt> object)
	 */
	LinOpNonPersisting( const LinOpPersisting<Scalar> &op );

	///
	/** Return raw pointer to underlying linear operator.
	 */
	const LinearOp<Scalar>* op() const;

	//@}

};

// ///////////////////////////////////
// Implementation

//
// LinOpBase
//

template<class Scalar>
inline
LinOpBase<Scalar>::LinOpBase()
	:defaultTrans_(NOTRANS)
{}

template<class Scalar>
inline
LinOpBase<Scalar>::LinOpBase(
	const Teuchos::RefCountPtr<const LinearOp<Scalar> >   &op
	,const ETransp                                        defaultTrans
	)
	:op_(op)
	,defaultTrans_(defaultTrans)
{}

template<class Scalar>
inline
void LinOpBase<Scalar>::LinOpBase<Scalar>::defaultTrans( const ETransp defaultTrans )
{
	defaultTrans_ = defaultTrans;
}

template<class Scalar>
inline
Teuchos::RefCountPtr<const LinearOp<Scalar> > LinOpBase<Scalar>::op() const
{
	return op_;
}

template<class Scalar>
inline
ETransp LinOpBase<Scalar>::defaultTrans() const
{
	return defaultTrans_;
}

template<class Scalar>
inline
Teuchos::RefCountPtr<const VectorSpace<Scalar> >
LinOpBase<Scalar>::range() const
{
	return ( defaultTrans_==NOTRANS ? op_->range() : op_->domain() );
}

template<class Scalar>
inline
Teuchos::RefCountPtr<const VectorSpace<Scalar> >
LinOpBase<Scalar>::domain() const
{
	return ( defaultTrans_==NOTRANS ? op_->domain() : op_->range() );
}

template<class Scalar>
inline
bool LinOpBase<Scalar>::opSupported(ETransp M_trans) const
{
	return op_->opSupported(trans_trans(M_trans,defaultTrans_));
}

template<class Scalar>
inline
void LinOpBase<Scalar>::apply(
	const ETransp            M_trans
	,const Vector<Scalar>    &x
	,Vector<Scalar>          *y
	,const Scalar            alpha
	,const Scalar            beta
	) const
{
	op_->apply(trans_trans(M_trans,defaultTrans_),x,y,alpha,beta);
}

template<class Scalar>
inline
void LinOpBase<Scalar>::apply(
	const ETransp                 M_trans
	,const MultiVector<Scalar>    &X
	,MultiVector<Scalar>          *Y
	,const Scalar                 alpha
	,const Scalar                 beta
	) const
{
	op_->apply(trans_trans(M_trans,defaultTrans_),X,Y,alpha,beta);
}

//
// LinOpPersisting
//

template<class Scalar>
inline
LinOpPersisting<Scalar>::LinOpPersisting()
	:LinOpBase<Scalar>()
{}

template<class Scalar>
inline
LinOpPersisting<Scalar>::LinOpPersisting(
	const Teuchos::RefCountPtr<const LinearOp<Scalar> >   &op
	,const ETransp                                        defaultTrans
	)
	:LinOpBase<Scalar>(op,defaultTrans)
{
	TEST_FOR_EXCEPT( op.get() == NULL );
}

template<class Scalar>
inline
Teuchos::RefCountPtr<const LinearOp<Scalar> >
LinOpPersisting<Scalar>::op() const
{
	return LinOpBase<Scalar>::op();
}

//
// LinOpNonPersisting
//

template<class Scalar>
inline
LinOpNonPersisting<Scalar>::LinOpNonPersisting()
	:LinOpBase<Scalar>()
{}

template<class Scalar>
inline
LinOpNonPersisting<Scalar>::LinOpNonPersisting(
	const LinearOp<Scalar>                                &op
	,const ETransp                                        defaultTrans
	)
	:LinOpBase<Scalar>(Teuchos::rcp(&op,false),defaultTrans)
{}

template<class Scalar>
inline
LinOpNonPersisting<Scalar>::LinOpNonPersisting( const LinOpPersisting<Scalar> &op )
	:LinOpBase<Scalar>(op.op(),op.defaultTrans())
{}

template<class Scalar>
inline
const LinearOp<Scalar>*
LinOpNonPersisting<Scalar>::op() const
{
	return LinOpBase<Scalar>::op().get();
}

}	// end namespace TSFCore

#endif	// TSFCORE_LIN_OP_HPP
