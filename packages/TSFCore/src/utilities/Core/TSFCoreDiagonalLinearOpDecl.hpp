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
// TSFCoreDiagonalLinearOpDecl.hpp

#ifndef TSFCORE_DIAGONAL_LINEAR_OP_DECL_HPP
#define TSFCORE_DIAGONAL_LINEAR_OP_DECL_HPP

#include "TSFCoreLinOp.hpp"

namespace TSFCore {

///
/** Concrete <tt>LinearOp</tt> subclass for diagonal linear operators.
 *
 * This class represents a diagonal linear operator <tt>M</tt> of the form:
 \verbatim
 
 M = diag(diag)

 \endverbatim
 *
 * where <tt>diag</tt> is a <tt>Vector</tt> object.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class DiagonalLinearOp : virtual public LinearOp<Scalar> {
public:

	/** @name Constructors/initializers/accessors */
	//@{

	///
	/** Constructs to uninitialized.
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getDiag().get()==NULL</tt>
	 * </ul>
	 */
	DiagonalLinearOp();

	/// Calls <tt>initialize()</tt>
	DiagonalLinearOp(
		const Teuchos::RefCountPtr<const Vector<Scalar> >   &diag
		);

	///
	/** Initialize given the diagonal.
	 *
	 * @param  diag  [in] Smart pointer to diagonal vector. 
	 *
	 * Precconditions:<ul>
	 * <li><tt>diag.get()!=NULL</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getDiag().get()==diag.get()</tt>
	 * <li><tt>this->this->domain().get() == diag->space().get()</tt>
	 * <li><tt>this->this->range().get() == diag->space().get()</tt>
	 * </ul>
	 */
	void initialize(
		const Teuchos::RefCountPtr<const Vector<Scalar> >   &diag
		);

	///
	/** Returns the current number of constutient operators.
	 *
	 * A return value of <tt>0</tt> indicates that <tt>this</tt> is not
	 * fully initialized.
	 */
	Teuchos::RefCountPtr<const Vector<Scalar> >  getDiag() const;

	///
	/** Set to uninitialized.
	 *
	 * @param  diag  [out] Optional pointer to smart pointer for diagonal.
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getDiag().get()==NULL</tt>
	 * </ul>
	 */
	void uninitialize(
		Teuchos::RefCountPtr<const Vector<Scalar> >  *diag = NULL
		);

	//@}

	/** @name Overridden from OpBase */
	//@{

	///
	/** Returns <tt>this->getDiag()->space()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getDiag().get()!=NULL</tt>
	 * </ul>
	 */
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > domain() const;
	///
	/** Returns <tt>this->getDiag()->space()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getDiag().get()!=NULL</tt>
	 * </ul>
	 */
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > range() const;
	///
	/** Returns <tt>true</tt>.
	 */
	bool opSupported(ETransp M_trans) const;

	//@}

	/** @name Overridden from LinearOp */
	//@{

	///
	void apply(
		const ETransp            M_trans
		,const Vector<Scalar>    &x
		,Vector<Scalar>          *y
		,const Scalar            alpha
		,const Scalar            beta
		) const;
	///
	void apply(
		const ETransp                 M_trans
		,const MultiVector<Scalar>    &X
		,MultiVector<Scalar>          *Y
		,const Scalar                 alpha
		,const Scalar                 beta
		) const;
	///
	Teuchos::RefCountPtr<const LinearOp<Scalar> > clone() const;

	//@}

private:

	Teuchos::RefCountPtr< const Vector<Scalar> > diag_;

	void assertInitialized() const;

};

// /////////////////////////////////
// Inline members

template<class Scalar>
inline
Teuchos::RefCountPtr<const Vector<Scalar> >  
DiagonalLinearOp<Scalar>::getDiag() const
{
	return diag_;
}

template<class Scalar>
inline
void DiagonalLinearOp<Scalar>::assertInitialized() const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( diag_.get() == NULL );
#endif
}

}	// end namespace TSFCore

#endif	// TSFCORE_DIAGONAL_LINEAR_OP_DECL_HPP
