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
// TSFCoreLinearOpHandle.hpp

#ifndef TSFCORE_LINEAR_OP_HANDLE_HPP
#define TSFCORE_LINEAR_OP_HANDLE_HPP

#include "TSFCoreLinearOp.hpp"

namespace TSFCore {

///
/** \brief Simple aggregation class of a <tt>LinearOp</tt>, its
 * "mathematical" transpose argument and its "mathematical" scalar
 * multiplier.
 *
 * This simple concrete handle-like class should be used by clients
 * that do not want to worrry about whether what the mathematical
 * definition of the un-transposed scaled operator is.  Using this
 * handle class removes the need to create "decorator" subclasses to
 * represent different transpose arguments or to scale the operator
 * and therefore preserves the ability to dynamically cast the
 * underlying <tt>LinearOp</tt> object to query it for extended
 * interfaces.
 *
 * \ingroup TSFCore_ANA_Development_grp
 */
template<class Scalar>
class LinearOpHandle {
public:

	/** @name Public constructors/accessors */
	//@{

	///
	/** \brief Default construct to NULL handle.
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->op().get() == NULL</tt>
	 * <li><tt>this->defaultTrans() == NOTRANS</tt>
	 * <li><tt>this->defaultAlpha() == Teuchos::ScalarTraits<Scalar>::one()</tt>
	 * </ul>
	 */
	LinearOpHandle();

	///
	/** \brief Construct with an operator and optionally its default transpose
	 * and scaling arguments.
	 *
	 * @param  op            [in] Smart pointer to linear operator (persisting relationship).
	 * @param  defaultTrans  [in] Default definition of transpose (default: <tt>NOTRANS</tt>).
	 * @param  defaultAlpha  [in] Default value for scalar <tt>alpha</tt>
	 *                            (default: <tt>Teuchos::ScalarTraits<Scalar>::one()</tt>).
	 *
	 *
	 * Preconditions:<ul>
	 * <li><tt>op.get() != NULL</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->op().get() == op.get()</tt>
	 * <li><tt>this->defaultTrans() == defaultTrans</tt>
	 * <li><tt>this->defaultAlpha() == defaultAlpha</tt>
	 * </ul>
	 *
	 * Note that this constructor defines an implicit conversion from a
	 * <tt>Teuchos::RefCountPtr<const LinearOp<Scalar> ></tt> object to
	 * a <tt>LinearOpHandle<Scalar></tt> object.  This is reasonable
	 * behavior and is therefore allowed.
	 */
	LinearOpHandle(
		const Teuchos::RefCountPtr<const LinearOp<Scalar> >   &op
		,const ETransp                                        defaultTrans  = NOTRANS
		,const Scalar                                         &defaultAlpha = Teuchos::ScalarTraits<Scalar>::one()
		);

	///
	/** Return smart pointer to underlying linear operator.
	 */
	Teuchos::RefCountPtr<const LinearOp<Scalar> > op() const;

	///
	/** Set the mathematical transpose argument
	 *
	 * @param  defaultTrans  [in] Default definition of transpose (default: <tt>NOTRANS</tt>).
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->defaultTrans() == defaultTrans</tt>
	 * </ul>
	 */
	void defaultTrans( const ETransp defaultTrans );

	///
	/** Return the default defintion of mathematical transpose.
	 */
	ETransp defaultTrans() const;

	///
	/** Set the default scalar multiplier.
	 *
	 * @param  defaultAlpha  [in] The default scalar multiplier <tt>alpha</tt>.
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->defaultAlpha() == defaultAlpha</tt>
	 * </ul>
	 */
	void defaultAlpha( const Scalar &defaultAlpha );

	///
	/** Return the default defintion of mathematical transpose.
	 */
	Scalar defaultAlpha() const;

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
	 * @param  M_trans [in] See <tt>LinearOp::apply()</tt>.
	 * @param  x       [in] See <tt>LinearOp::apply()</tt>.
	 * @param  y       [in/out] See <tt>LinearOp::apply()</tt>.
	 * @param  alpha   [in] See <tt>LinearOp::apply()</tt>.
	 * @param  beta    [in] See <tt>LinearOp::apply()</tt>.
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
   this->op()->apply(trans_trans(M_trans,this->defaultTrans()),x,y,(this->defaultAlpha()*alpha),beta)
   \endcode
	 */
	void apply(
		const ETransp            M_trans
		,const Vector<Scalar>    &x
		,Vector<Scalar>          *y
		,const Scalar            alpha = Teuchos::ScalarTraits<Scalar>::one()
		,const Scalar            beta  = Teuchos::ScalarTraits<Scalar>::zero()
		) const;

	///
	/** Apply the linear operator (or its transpose) to a multi-vector :
	 * <tt>Y = alpha*op(M)*X + beta*Y</tt>.
	 *
	 * @param  M_trans [in] See <tt>LinearOp::apply()</tt>.
	 * @param  X       [in] See <tt>LinearOp::apply()</tt>.
	 * @param  Y       [in/out] See <tt>LinearOp::apply()</tt>.
	 * @param  alpha   [in] See <tt>LinearOp::apply()</tt>.
	 * @param  beta    [in] See <tt>LinearOp::apply()</tt>.
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
   this->op()->apply(trans_trans(M_trans,this->defaultTrans()),X,Y,(this->defaultAlpha()*alpha),beta)
	 \endcode
	 */
	void apply(
		const ETransp                 M_trans
		,const MultiVector<Scalar>    &X
		,MultiVector<Scalar>          *Y
		,const Scalar                 alpha = Teuchos::ScalarTraits<Scalar>::one()
		,const Scalar                 beta  = Teuchos::ScalarTraits<Scalar>::zero()
		) const;

	//@}

private:

	Teuchos::RefCountPtr<const LinearOp<Scalar> >  op_;
	ETransp                                        defaultTrans_;
	Scalar                                         defaultAlpha_;

};

// ///////////////////////////////////
// Implementation

template<class Scalar>
inline
LinearOpHandle<Scalar>::LinearOpHandle()
	:defaultTrans_(NOTRANS),defaultAlpha_(Teuchos::ScalarTraits<Scalar>::one())
{}

template<class Scalar>
inline
LinearOpHandle<Scalar>::LinearOpHandle(
	const Teuchos::RefCountPtr<const LinearOp<Scalar> >   &op
	,const ETransp                                        defaultTrans
	,const Scalar                                         &defaultAlpha
	)
	:op_(op)
	,defaultTrans_(defaultTrans)
	,defaultAlpha_(defaultAlpha)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT(op.get()==NULL);
#endif
}

template<class Scalar>
inline
void LinearOpHandle<Scalar>::defaultTrans( const ETransp defaultTrans )
{
	defaultTrans_ = defaultTrans;
}

template<class Scalar>
inline
ETransp LinearOpHandle<Scalar>::defaultTrans() const
{
	return defaultTrans_;
}

template<class Scalar>
inline
void LinearOpHandle<Scalar>::defaultAlpha( const Scalar &defaultAlpha )
{
	defaultAlpha_ = defaultAlpha_;
}

template<class Scalar>
inline
Scalar LinearOpHandle<Scalar>::defaultAlpha() const
{
	return defaultAlpha_;
}

template<class Scalar>
inline
Teuchos::RefCountPtr<const LinearOp<Scalar> > LinearOpHandle<Scalar>::op() const
{
	return op_;
}

template<class Scalar>
inline
Teuchos::RefCountPtr<const VectorSpace<Scalar> >
LinearOpHandle<Scalar>::range() const
{
	return ( defaultTrans_==NOTRANS ? op_->range() : op_->domain() );
}

template<class Scalar>
inline
Teuchos::RefCountPtr<const VectorSpace<Scalar> >
LinearOpHandle<Scalar>::domain() const
{
	return ( defaultTrans_==NOTRANS ? op_->domain() : op_->range() );
}

template<class Scalar>
inline
bool LinearOpHandle<Scalar>::opSupported(ETransp M_trans) const
{
	return op_->opSupported(trans_trans(M_trans,defaultTrans_));
}

template<class Scalar>
inline
void LinearOpHandle<Scalar>::apply(
	const ETransp            M_trans
	,const Vector<Scalar>    &x
	,Vector<Scalar>          *y
	,const Scalar            alpha
	,const Scalar            beta
	) const
{
	op_->apply(trans_trans(M_trans,defaultTrans_),x,y,(defaultAlpha_*alpha),beta);
}

template<class Scalar>
inline
void LinearOpHandle<Scalar>::apply(
	const ETransp                 M_trans
	,const MultiVector<Scalar>    &X
	,MultiVector<Scalar>          *Y
	,const Scalar                 alpha
	,const Scalar                 beta
	) const
{
	op_->apply(trans_trans(M_trans,defaultTrans_),X,Y,(defaultAlpha_*alpha),beta);
}

}	// end namespace TSFCore

#endif	// TSFCORE_LINEAR_OP_HANDLE_HPP
