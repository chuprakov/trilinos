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
template<class Scalar> class MPIVectorSpaceBase;

///
/** Base class for MPI-based multi-vectors.
 *
 * ToDo: Finish documentation!
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
    const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    &mpiRangeSpace
    ,const Teuchos::RefCountPtr<const VectorSpace<Scalar> >          &domainSpace
    ,const Teuchos::RefCountPtr<Scalar>                              &values
    ,const Index                                                     leadingDim
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
   * @param  values    [in] Smart pointer to beginning of Fortran-style column-major
   *                   array that defines the local values in the multi-vector.
   *                   This array must be at least of dimension <tt>mpiRangeSpace->leadingDim()*domainSpace->dim()</tt>
   *                   and <tt>(&*values)[ (i-1) + (j-1)*leadingDim ]</tt> gives the local value
   *                   of the one-based <tt>(i,j)</tt> entry where <tt>i=1...mpiSpace()->localSubDim()</tt>
   *                   and <tt>j=1...domainSpace->dim()</tt>.
   *
   * Preconditions:<ul>
   * <li><tt>mpiRangeSpace.get()!=NULL</tt>
   * <li><tt>domainSpace.get()!=NULL</tt>
   * <li><tt>values.get()!=NULL</tt>
   * <li><tt>leadingDim >= mpiRangeSpace->localSubDim()</tt>
   * </ul>
   */
	void initialize(
    const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    &mpiRangeSpace
    ,const Teuchos::RefCountPtr<const VectorSpace<Scalar> >          &domainSpace
    ,const Teuchos::RefCountPtr<Scalar>                              &values
    ,const Index                                                     leadingDim
    );

  ///
  /** Set to an uninitialized state.
   *
   * Postconditions:<ul>
   * <li><tt>this->mpiSpace().get() == NULL</tt>.
   */
	void uninitialize(
    Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    *mpiRangeSpace = NULL
    ,Teuchos::RefCountPtr<const VectorSpace<Scalar> >          *domainSpace   = NULL
    ,Teuchos::RefCountPtr<Scalar>                              *values        = NULL
    ,Index                                                     *leadingDim    = NULL
    );

  //@}

	/** @name Overridden from OpBase */
	//@{
	///
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > domain() const;
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
	void getLocalData( const Scalar **values, Index *leadingDim ) const;
	///
	void freeLocalData( const Scalar *values ) const;
	///
	void getLocalData( Scalar **values, Index *leadingDim );
	///
	void commitLocalData( Scalar *values );
	//@}
	
private:
	
	// ///////////////////////////////////////
	// Private data members

  Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    mpiRangeSpace_;
  Teuchos::RefCountPtr<const VectorSpace<Scalar> >           domainSpace_;
  Teuchos::RefCountPtr<Scalar>                               values_;
  Index                                                      leadingDim_;
	
}; // end class MPIMultiVectorStd

} // end namespace TSFCore

#endif // TSFCORE_MPI_MULTI_VECTOR_STD_DECL_HPP
