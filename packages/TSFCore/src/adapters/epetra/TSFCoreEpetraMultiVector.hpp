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

// /////////////////////////////////////////////////////////////////////////////
// TSFCoreEpetraMultiVector.hpp

#ifndef TSFCORE_EPETRA_MULTI_VECTOR_HPP
#define TSFCORE_EPETRA_MULTI_VECTOR_HPP

#include "TSFCoreEpetraTypes.hpp"
#include "TSFCoreMPIMultiVectorBase.hpp"

namespace TSFCore {

/** \brief Concrete <tt>MultiVector</tt> adapter subclass for
 * <tt>Epetra_MultiVector</tt>.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup TSFCore_adapters_Epetra_grp
 */
class EpetraMultiVector : virtual public MPIMultiVectorBase<RTOp_value_type> {
public:
	
	///
	typedef RTOp_value_type Scalar;
	///
	using MultiVector<Scalar>::col;     // Inject *all* functions!
	///
	using MultiVector<Scalar>::subView; // Inject *all* functions!

	/** @name Constructors/Initializers */
	//@{

	/** Construct to initialized.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->range().get() == NULL</tt>
	 * <li> <tt>this->domain().get() == NULL</tt>
	 * </ul>
	 */
	EpetraMultiVector();

	/// Calls <tt>initalize()</tt>.
	EpetraMultiVector(
		const Teuchos::RefCountPtr<Epetra_MultiVector>                        &epetra_multi_vec
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>                  &epetra_range
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>                  &epetra_domain     = Teuchos::null
#else
		,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> > &domain            = Teuchos::null
#endif
		);
	
	/** Initialize given the <tt>Epetra_multi_vec</tt>.
   *
   * ToDo: Finish documentation!
	 */
	void initialize(
		const Teuchos::RefCountPtr<Epetra_MultiVector>                        &epetra_multi_vec
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>                  &epetra_range
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>                  &epetra_domain     = Teuchos::null
#else
		,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> > &domain            = Teuchos::null
#endif
		);
	
  /** Set uninitalized.
   *
   * ToDo: Finish documentation!
   */
	void setUninitialized(
		Teuchos::RefCountPtr<Epetra_MultiVector>                        *epetra_multi_vec = NULL
		,Teuchos::RefCountPtr<const EpetraVectorSpace>                  *epetra_range     = NULL
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
		,Teuchos::RefCountPtr<const EpetraVectorSpace>                  *epetra_domain    = NULL
#else
		,Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> > *domain           = NULL
#endif
		);

	/** Get a smart pointer to non-<tt>const</tt> <tt>Epetra_MultiVector</tt> object.
	 *
	 * Note: Only the numerical values of the Epetra vector should be changed by the
	 * client.  Any other type of change will invalidate <tt>this</tt>.
	 */
 	Teuchos::RefCountPtr<Epetra_MultiVector> epetra_multi_vec();

	/** Get a smart pointer to <tt>const</tt> <tt>Epetra_MultiVector</tt> object.
	 *
	 * This gives access to read the <tt>Epetra_MultiVector</tt> object only
	 * so this should be very safe.
	 */
	Teuchos::RefCountPtr<const Epetra_MultiVector> epetra_multi_vec() const;

	//@}

	/** @name Overridden from EuclideanLinearOpBase */
	//@{
	/// Overridden
	Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> > domainScalarProdVecSpc() const;
#if defined(TSFCORE_EPETRA_USE_EPETRA_MULTI_VECTOR_MULTIPLY) && defined(TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE)
	/// Overridden
	void euclideanApply(
		const ETransp                 M_trans
		,const MultiVector<Scalar>    &X
		,MultiVector<Scalar>          *Y
		,const Scalar                 alpha
		,const Scalar                 beta
		) const;
#endif
	//@}

	/** @name Overridden from MultiVector */
	//@{
	/// Overridden
	Teuchos::RefCountPtr<Vector<Scalar> > col(Index j);
	/// Overridden
	Teuchos::RefCountPtr<MultiVector<Scalar> > subView( const Range1D& col_rng );
	/// Overridden
	Teuchos::RefCountPtr<MultiVector<Scalar> > subView( const int numCols, const int cols[] );
	//@}

	/** @name Overridden from MPIMultiVectorBase */
	//@{
	/// Overridden
	Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> > mpiSpace() const;
	/// Overridden
	void getLocalData( const Scalar **values, Index *leadingDim ) const;
	/// Overridden
	void freeLocalData( const Scalar *values ) const;
	/// Overridden
	void getLocalData( Scalar **values, Index *leadingDim );
	/// Overridden
	void commitLocalData( Scalar *values );
	//@}

private:
	
	Teuchos::RefCountPtr<Epetra_MultiVector>                         epetra_multi_vec_;
	Teuchos::RefCountPtr<const EpetraVectorSpace>                    epetra_range_;
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
	Teuchos::RefCountPtr<const EpetraVectorSpace>                    epetra_domain_;
#else
	Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   domain_;
#endif
	
}; // end class EpetraMultiVector

// ///////////////////////////////////
// Inline members

inline
Teuchos::RefCountPtr<Epetra_MultiVector>
EpetraMultiVector::epetra_multi_vec()
{
	return epetra_multi_vec_;
}

inline
Teuchos::RefCountPtr<const Epetra_MultiVector>
EpetraMultiVector::epetra_multi_vec() const
{
	return epetra_multi_vec_;
}

} // end namespace TSFCore

#endif // TSFCORE_EPETRA_MULTI_VECTOR_HPP
