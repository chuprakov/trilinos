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
// TSFCoreMultiplicativeLinearOpDecl.hpp

#ifndef TSFCORE_MULTIPLICATIVE_LINEAR_OP_DECL_HPP
#define TSFCORE_MULTIPLICATIVE_LINEAR_OP_DECL_HPP

#include "TSFCoreLinOp.hpp"

namespace TSFCore {

///
/** Concrete composite <tt>LinearOp</tt> subclass that creates a
 * multiplicative linear operator out of one or more constituent
 * <tt>#LinearOp</tt> objects.
 *
 * This class represents a multiplicative linear operator <tt>M</tt> of the form:
 \verbatim
 
 M = Op[0] * Op[1] * ... * Op[numOps-1]

 \endverbatim
 *
 * where <tt>Op[]</tt> is an array of <tt>numOps</tt>
 * <tt>LinOpPersisting</tt> objects.  Of course the operator
 * <tt>M</tt> is not constructed explicitly but instead just applies
 * the constituent linear operators accordingly using temporaries.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class MultiplicativeLinearOp : virtual public LinearOp<Scalar> {
public:

	/** @name Constructors/initializers/accessors */
	//@{

	///
	/** Constructs to uninitialized.
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->numOps()==0</tt>
	 * </ul>
	 */
	MultiplicativeLinearOp();

	/// Calls <tt>initialize()</tt>
	MultiplicativeLinearOp(
		const int                        numOps
		,const LinOpPersisting<Scalar>   Ops[]
		);

	///
	/** Initialize given a list of linear operators.
	 *
	 * @param  numOps  [in] Number of constituent opeators.
	 * @param  Ops     [in] Array (length <tt>numOps</tt>) of
	 *                 constituent linear operators and their
	 *                 aggregated default definitions of the
	 *                 non-transposed operator.
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->numOps()==numOps</tt>
	 * <li><tt>this->getOp(k).op().get()==Ops[k].op().get()</tt>
	 * </ul>
	 */
	void initialize(
		const int                        numOps
		,const LinOpPersisting<Scalar>   Ops[]
		);

	///
	/** Returns the current number of constutient operators.
	 *
	 * A return value of <tt>0</tt> indicates that <tt>this</tt> is not
	 * fully initialized.
	 */
	int numOps() const;

	///
	/** Return the <tt>k</tt>th constituent operator.
	 *
	 * @param  k  [in] The zero-based index of the constituent operator to return.
	 *
	 * Preconditions:<ul>
	 * <li><tt> 0 <= k < this->numOps()</tt>
	 * </ul>
	 */
	const LinOpPersisting<Scalar>& getOp(const int k) const;

	///
	/** Set to uninitialized.
	 *
	 * @param  numOps  [in] Number of operators (must be <tt>numOps==this->numOps()</tt>).
	 * @param  Ops     [out] Array (length <tt>numOps</tt>) that if <tt>Ops!=NULL</tt>
	 *                 then <tt>Ops[k]</tt> will be set to <tt>this->getOp(k)</tt>, for
	 *                 <tt>k=0...numOps-1</tt>.
	 *
	 * Precconditions:<ul>
	 * <li>[<tt>Ops!=NULL</tt>] <tt>numOps==this->numOps()</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->numOps()==0</tt>
	 * </ul>
	 */
	void uninitialize(
		const int                  numOps   = 0
		,LinOpPersisting<Scalar>   Ops       = NULL
		);

	//@}

	/** @name Overridden from OpBase */
	//@{

	///
	/** Returns <tt>this->getOp(this->numOps()-1).domain()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->numOps()==0</tt>
	 * </ul>
	 */
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > domain() const;
	///
	/** Returns <tt>this->getOp(0).range()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->numOps()==0</tt>
	 * </ul>
	 */
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > range() const;
	///
	/** Returns <tt>true</tt> only if all constituent opeators support
	 * <tt>M_trans</tt>.
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

	std::vector<LinOpPersisting<Scalar> >  Ops_;

	void assertInitialized() const;

};

// /////////////////////////////////
// Inline members

template<class Scalar>
inline
int MultiplicativeLinearOp<Scalar>::numOps() const
{
	return Ops_.size();
}

template<class Scalar>
inline
const LinOpPersisting<Scalar>& MultiplicativeLinearOp<Scalar>::getOp(const int k) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( !( 0 <= k && k < numOps() ) );
#endif
	return Ops_[k];
}

template<class Scalar>
inline
void MultiplicativeLinearOp<Scalar>::assertInitialized() const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( !( numOps() > 0 ) );
#endif
}

}	// end namespace TSFCore

#endif	// TSFCORE_MULTIPLICATIVE_LINEAR_OP_DECL_HPP
