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
// TSFCoreLinearOpDecl.hpp

#ifndef TSFCORE_LINEAR_OP_DECL_HPP
#define TSFCORE_LINEAR_OP_DECL_HPP

#include "TSFCoreOpBase.hpp"

namespace TSFCore {

///
/** Base class for all linear operators.
 *
 * A linear operator can perform the following operations:
 *
 * <ul>
 * <li><tt>y = alpha*op(M)*x + beta*y</tt>  (vector version)
 * <li><tt>Y = alpha*op(M)*X + beta*Y</tt>  (multi-vector version)
 * </ul>
 *
 * through the <tt>apply()</tt> methods where <tt>y</tt> and
 * <tt>x</tt> are <tt>Vector</tt> objects while <tt>Y</tt> and
 * <tt>X</tt> are <tt>MultiVector</tt> objects.  The reason for the
 * exact form of the above operations is that there are direct BLAS
 * and equivalent versions of these operations.
 *
 * A linear operator has <tt>domain()</tt> and <tt>range()</tt> vector
 * spaces associated with it for the vectors <tt>x</tt> and <tt>y</tt>
 * that lie in the domain and the range spaces of the non-transposed
 * linear operator <tt> y = M*x </tt>.
 *
 * Note that it is strictly forbidden to alias the input/output
 * objects <tt>y</tt> and <tt>Y</tt> with the input objects <tt>x</tt>
 * and <tt>X</tt>.
 *
 * If a <tt>%LinearOp</tt> subclass can not support a particular
 * value of <tt>M_tans</tt> in the <tt>apply()</tt> methods, then the
 * method <tt>opSupported()</tt> returns <tt>false</tt> for that
 * particular value of <tt>M_trans</tt>.
 *
 * <b>Notes for subclass developers</b>
 *
 * There are only three methods that a subclass is required to
 * override: <tt>domain()</tt>, <tt>range()</tt> and <tt>apply()</tt>.
 * Note that the methods <tt>domain()</tt> and <tt>range()</tt> should
 * simply return <tt>VectorSpace</tt> objects for subclasses that are
 * already defined for the vectors that the linear operator interacts
 * with through the method <tt>apply()</tt>.  Therefore, given that
 * approprate <tt>VectorSpace</tt> and <tt>Vector</tt> subclasses exist,
 * The only real work involved in implementing a <tt>LinearOp</tt> 
 * is defining a single method <tt>apply()</tt>.
 *
 * If a <tt>LinearOp</tt> subclass can not support a particular value
 * of the transpose argument <tt>M_trans</tt> in the <tt>apply()</tt>
 * methods, then the method <tt>opSuported(M_trans)</tt> must be
 * overriden to return <tt>false</tt> for this value of
 * <tt>M_trans</tt>.
 *
 * If possible, the subclass should also override the <tt>clone()</tt>
 * method with allows clients to create copies of a <tt>LinearOp</tt>
 * object.  This functionality is very important in some
 * circumstances.  However, this functionality is not required and
 * <tt>clone()</tt> returns a null smart pointer object.
 *
 * If multi-vectors are supported in general by the application and
 * linear algebra library then, if possible, the subclass should also
 * override the multi-vector version of <tt>apply()</tt>.  In many
 * cases, a specialized multi-vector version will outperform the
 * default implementation (which is based on the single vector
 * version) in this class.
 */
template<class Scalar>
class LinearOp : virtual public OpBase<Scalar> {
public:

	/** @name Pure virtual methods (must be overridden by subclass) */

	///
	/** Apply the linear operator (or its transpose) to a vector:
	 * <tt>y = alpha*op(M)*x + beta*y</tt>.
	 *
	 * @param  M_trans
	 *                [in] Determines whether the transposed or non-trnasposed
	 *                operator is applied as:
	 *                <ul>
	 *                <li> <tt>op(M) = M</tt>, for <tt>M_trans==NOTRANS</tt>
	 *                <li> <tt>op(M) = M'</tt>, for <tt>M_trans==TRANS</tt>
	 *                </ul>
	 *                where <tt>M == *this</tt>
	 * @param  x      [in] The right hand side vector 
	 * @param  y      [in/out] The target vector being transformed
	 * @param  alpha  [in] Scalar multiplying <tt>M</tt>, where <tt>M==*this</tt>.
     *                The default value of <tt>alpha</tt> is </tt>1.0</tt>
	 * @param  beta   [in] The multiplier for the target vector <tt>y</tt>.
	 *                The default value of <tt>beta</tt> is <tt>0.0</tt>.
	 * 
	 * Preconditions:<ul>
	 * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> <tt>this->opSupported(M_trans)==true</tt> (throw <tt>Exceptions::OpNotSupported</tt>)
	 * <li> <tt>y->space()->isCompatible(M_trans==NOTRANS ? *this->range() : *this->domain()) == true</tt>
	 *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>x.space()->isCompatible(M_trans==NOTRANS ? *this->domain() : *this->range()) == true</tt>
	 *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>y</tt> can not alias <tt>x</tt>.  It is up to the client to ensure that <tt>y</tt>
	 *      and <tt>x</tt> are distinct since in general this can not be verified by the implementation until,
	 *      perhaps, it is too late.  If possible, an exception will be thrown if aliasing is detected.
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> Is it not obvious?  After the method returns the vector <tt>y</tt>
	 *      is transformed as indicated above.
	 * </ul>
	 */
	virtual void apply(
		const ETransp            M_trans
		,const Vector<Scalar>    &x
		,Vector<Scalar>          *y
		,const Scalar            alpha = 1.0
		,const Scalar            beta  = 0.0
		) const = 0;

	//@}

	/** @name Virtual functions with default implemenations */
	//@{

	///
	/** Clone the linear operator object (if supported).
	 *
	 * The primary purpose for this method is to allow a client to
	 * capture the current state of a linear operator object and be
	 * guaranteed that some other client will not alter its behavior.
	 * A smart implementation will use reference counting and lazy
	 * evaluation internally and will not actually copy any large
	 * amount of data unless it has to.
	 *
	 * The default implementation returns <tt>return.get()==NULL</tt>
	 * which is allowable by this specification.  A linear operator
	 * object is not required to return a non-NULL value but almost
	 * every good linear operator implementation should and will.
	 */
	virtual Teuchos::RefCountPtr<const LinearOp<Scalar> > clone() const;

	///
	/** Apply the linear operator (or its transpose) to a multi-vector :
	 * <tt>Y = alpha*op(M)*X + beta*Y</tt>.
	 *
	 * @param  M_trans
	 *                [in] Determines whether the transposed or non-trnasposed
	 *                operator is applied as:
	 *                <ul>
	 *                <li> <tt>op(M) = M</tt>, for <tt>M_trans==NOTRANS</tt>
	 *                <li> <tt>op(M) = M'</tt>, for <tt>M_trans==TRANS</tt>
	 *                </ul>
	 *                where <tt>M == *this</tt>
	 * @param  X      [in] The right hand side multi-vector 
	 * @param  Y      [in/out] The target multi-vector being transformed
	 * @param  alpha  [in] Scalar multiplying <tt>M</tt>, where <tt>M==*this</tt>.
     *                The default value of <tt>alpha</tt> is </tt>1.0</tt>
	 * @param  beta   [in] The multiplier for the target multi-vector <tt>y</tt>.
	 *                The default value of <tt>beta</tt> is <tt>0.0</tt>.
	 * 
	 * Preconditions:<ul>
	 * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> <tt>this->opSupported(M_trans)==true</tt> (throw <tt>Exceptions::OpNotSupported</tt>)
	 * <li> <tt>Y->range()->isCompatible(M_trans==NOTRANS ? *this->range() : *this->domain()) == true</tt>
	 *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>Y->domain()->isCompatible(*X.domain()) == true</tt>
	 *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>X.range()->isCompatible(M_trans==NOTRANS ? *this->domain() : *this->range()) == true</tt>
	 *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>Y</tt> can not alias <tt>X</tt>.  It is up to the client to ensure that <tt>Y</tt>
	 *      and <tt>X</tt> are distinct since in general this can not be verified by the implementation until,
	 *      perhaps, it is too late.  If possible, an exception will be thrown if aliasing is detected.
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> Is it not obvious?  After the method returns the multi-vector <tt>Y</tt>
	 *      is transformed as indicated above.
	 * </ul>
	 *
	 * This method has a default implementation in terms of the
	 * <tt>apply()</tt> method for vectors.
	 */
	virtual void apply(
		const ETransp                 M_trans
		,const MultiVector<Scalar>    &X
		,MultiVector<Scalar>          *Y
		,const Scalar                 alpha = 1.0
		,const Scalar                 beta  = 0.0
		) const;

	//@}

};	// end class LinearOp

}	// end namespace TSFCore

#endif	// TSFCORE_LINEAR_OP_DECL_HPP
