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
// TSFCoreMPIMultiVectorStdDecl.hpp

#ifndef TSFCORE_MPI_MULTI_VECTOR_STD_DECL_HPP
#define TSFCORE_MPI_MULTI_VECTOR_STD_DECL_HPP

#include "TSFCoreMPIMultiVectorBase.hpp"

namespace TSFCore {

///
/** Efficient concrete implementation subclass for SPMD-MPI-based multi-vectors.
 *
 * This subclass provides a very efficient and very general concrete
 * implementation of a <tt>TSFCore::MultiVector</tt> object.
 *
 * Objects of this type generally should not be constructed directly
 * by a client but instead by using the concrete vector space subclass
 * <tt>TSFCore::MPIVectorSpaceStd</tt> and using the function
 * <tt>TSFCore::MPIVectorSpaceStd::createMembers()</tt>.
 *
 * The storage type can be anything since a
 * <tt>Teuchos::RefCountPtr</tt> is used to pass in the local values
 * pointer into the constructor and <tt>initialize()</tt>.
 *
 * \ingroup TSFCore_adapters_MPI_concrete_std_grp
 */
template<class Scalar>
class MPIMultiVectorStd : virtual public MPIMultiVectorBase<Scalar> {
public:

  ///
  using MPIMultiVectorBase<Scalar>::subView;
  ///
  using MPIMultiVectorBase<Scalar>::col;

  /** @name Constructors/initializers/accessors */
  //@{

	/// Construct to uninitialized
	MPIMultiVectorStd();

  /// Calls <tt>initialize()</tt>
	MPIMultiVectorStd(
    const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >           &mpiRangeSpace
    ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &domainSpace
    ,const Teuchos::RefCountPtr<Scalar>                                     &localValues
    ,const Index                                                            leadingDim
    );

  ///
  /** Initialize.
   *
   * @param  mpiRangeSpace
   *                   [in] Smart pointer to <tt>MPIVectorSpaceBase</tt> object
   *                   that defines the data distribution for <tt>mpiSpace()</tt> and <tt>range()</tt>.
   * @param  domainSpace
   *                   [in] Smart pointer to <tt>VectorSpace</tt> object
   *                   that defines <tt>domain()</tt> space.
   * @param  localValues
	 *                   [in] Smart pointer to beginning of Fortran-style column-major
   *                   array that defines the local localValues in the multi-vector.
   *                   This array must be at least of dimension <tt>mpiRangeSpace->localSubDim()*domainSpace->dim()</tt>
   *                   and <tt>(&*localValues)[ (i-1) + (j-1)*leadingDim ]</tt> gives the local value
   *                   of the one-based entry <tt>(i,j)</tt> where <tt>i=1...mpiSpace()->localSubDim()</tt>
   *                   and <tt>j=1...domainSpace->dim()</tt>.
	 * @param  leadingDim
	 *                   [in] The leading dimension of the multi-vector.
   *
   * Preconditions:<ul>
   * <li><tt>mpiRangeSpace.get()!=NULL</tt>
   * <li><tt>domainSpace.get()!=NULL</tt>
   * <li><tt>localValues.get()!=NULL</tt>
   * <li><tt>leadingDim >= mpiRangeSpace->localSubDim()</tt>
   * </ul>
   */
	void initialize(
    const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >          &mpiRangeSpace
    ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domainSpace
    ,const Teuchos::RefCountPtr<Scalar>                                    &localValues
    ,const Index                                                           leadingDim
    );

  ///
  /** Set to an uninitialized state.
   *
   * Postconditions:<ul>
   * <li><tt>this->mpiSpace().get() == NULL</tt>.
   */
	void uninitialize(
    Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >                  *mpiRangeSpace = NULL
    ,Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >          *domainSpace   = NULL
    ,Teuchos::RefCountPtr<Scalar>                                            *localValues   = NULL
    ,Index                                                                   *leadingDim    = NULL
    );

  //@}

	/** @name Overridden from EuclideanLinearOpBase */
	//@{
	///
	Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> > domainScalarProdVecSpc() const;
  //@}

	/** @name Overridden from MultiVector */
	//@{
	///
	Teuchos::RefCountPtr<Vector<Scalar> > col(Index j);
	///
	Teuchos::RefCountPtr<MultiVector<Scalar> > subView( const Range1D& col_rng );
	//@}

	/** @name Overridden from MPIMultiVectorBase */
	//@{

	///
	Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> > mpiSpace() const;
	///
	void getLocalData( Scalar **localValues, Index *leadingDim );
	///
	void commitLocalData( Scalar *localValues );
	///
	void getLocalData( const Scalar **localValues, Index *leadingDim ) const;
	///
	void freeLocalData( const Scalar *localValues ) const;
	//@}
	
private:
	
	// ///////////////////////////////////////
	// Private data members

  Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    mpiRangeSpace_;
  Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >           domainSpace_;
  Teuchos::RefCountPtr<Scalar>                               localValues_;
  Index                                                      leadingDim_;
	
}; // end class MPIMultiVectorStd

} // end namespace TSFCore

#endif // TSFCORE_MPI_MULTI_VECTOR_STD_DECL_HPP
