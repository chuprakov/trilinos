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

#include "TSFCoreLinearOpHandle.hpp"

namespace TSFCore {

///
/** Concrete composite <tt>LinearOp</tt> subclass that creates a
 * multiplicative linear operator out of one or more constituent
 * <tt>LinearOp</tt> objects.
 *
 * This class represents a multiplicative linear operator <tt>M</tt> of the form:
 \verbatim
 
 M = gamma * Op[0] * Op[1] * ... * Op[numOps-1]
 \endverbatim
 *
 * where <tt>Op[]</tt> is an array of <tt>numOps</tt>
 * <tt>LinearOpHandle</tt> objects and <tt>gamma</tt> is a scalar.
 * Of course the operator <tt>M</tt> is not constructed explicitly but
 * instead just applies the constituent linear operators accordingly
 * using temporaries.
 *
 * In other words, this class defines <tt>apply()</tt> as:
 *
 \verbatim

 y = alpha*M*x + beta*y
   = (alpha*gamma) * ( Op[0] * ( Op[1] * ( .... ( Op[numOps-1] * x ) ... ) ) ) + beta * y
 \endverbatim
 *
 * for the case where <tt>M_trans==NOTRANS</tt> and as:
 *
 \verbatim

 y = alpha*M'*x + beta*y
   = (alpha*gamma) * ( Op[numOps-1]' * ( Op[numOps-2]' * ( .... ( Op[0]' * x ) ... ) ) ) + beta * y
 \endverbatim
 *
 * for the case where <tt>M_trans!=NOTRANS</tt> (where the transpose
 * <tt>'</tt> either defines <tt>TRANS</tt> or <tt>CONJTRANS</tt>).
 *
 * Constructing a multiplicative operator is easy.  For example, suppose one
 * wants to construct the multiplicative operator <tt>D = gamma * A * B' * C</tt>.
 * To do so one would do:

 \code
 template<class Scalar>
 void constructD(
    const Scalar                                                   &gamma
    ,const Teuchos::RefCountPtr<const TSFCore::LinearOp<Scalar> >  &A
    ,const Teuchos::RefCountPtr<const TSFCore::LinearOp<Scalar> >  &B
    ,const Teuchos::RefCountPtr<const TSFCore::LinearOp<Scalar> >  &C
    ,Teuchos::RefCountPtr<const TSFCore::LinearOp<Scalar> >        *D
    )
 {
   typedef TSFCore::LinearOpHandle<Scalar> LOP;
   *D = Teuchos::rcp(
     new TSFCore::MultiplicativeLinearOp<Scalar>(
       3, Teuchos::arrayArg<LOP>(LOP(A),LOP(B,TSFCore::TRANS),LOP(C))(), gamma
       )
     );
 }
 \endcode
 *

 *
 * \ingroup TSFCore_ANA_Development_grp
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
		,const LinearOpHandle<Scalar>    Ops[]
		,const Scalar                    &gamma = Teuchos::ScalarTraits<Scalar>::one()
		);

	///
	/** Initialize given a list of linear operators.
	 *
	 * @param  numOps  [in] Number of constituent opeators.
	 * @param  Ops     [in] Array (length <tt>numOps</tt>) of
	 *                 constituent linear operators and their
	 *                 aggregated default definitions of the
	 *                 non-transposed operator.
	 * @param  gamma   [in] Scalar multiplier
	 *
	 * Preconditions:<ul>
	 * <li><tt>numOps > 0</tt>
	 * <li><tt>Ops != NULL</tt>
	 * <li><tt>Ops[k].op().get()!=NULL</tt>, for <tt>k=0...numOps-1</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->numOps()==numOps</tt>
	 * <li><tt>this->getOp(k).op().get()==Ops[k].op().get()</tt>, for <tt>k=0...numOps-1</tt>
	 * <li><tt>this->gamma()==gamma</tt>
	 * </ul>
	 */
	void initialize(
		const int                        numOps
		,const LinearOpHandle<Scalar>    Ops[]
		,const Scalar                    &gamma = Teuchos::ScalarTraits<Scalar>::one()
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
	const LinearOpHandle<Scalar>& getOp(const int k) const;

	///
	/** Set to uninitialized.
	 *
	 * @param  numOps  [in] Number of operators (must be <tt>numOps==this->numOps()</tt>).
	 * @param  Ops     [out] Array (length <tt>numOps</tt>) that if <tt>Ops!=NULL</tt>
	 *                 then <tt>Ops[k]</tt> will be set to <tt>this->getOp(k)</tt>, for
	 *                 <tt>k=0...numOps-1</tt>.
	 * @param  gamma [out] Optional pointer to scalar <tt>gamma</tt>.
	 *               If <tt>gamma!=NULL</tt> then on output <tt>*gamma</tt>
	 *               is set to <tt>this->gamma()</tt> (before call).
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
		,LinearOpHandle<Scalar>    Ops[]    = NULL
		,Scalar                    *gamma   = NULL
		);

	//@}

	/** @name Overridden from OpBase */
	//@{
	///
	/** Returns <tt>this->getOp(0).range()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->numOps()==0</tt>
	 * </ul>
	 */
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > range() const;
	///
	/** Returns <tt>this->getOp(this->numOps()-1).domain()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->numOps()==0</tt>
	 * </ul>
	 */
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > domain() const;
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

	std::vector<LinearOpHandle<Scalar> >   Ops_;
	Scalar                                 gamma_;

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
const LinearOpHandle<Scalar>& MultiplicativeLinearOp<Scalar>::getOp(const int k) const
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
