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
// TSFCoreMPIMultiVectorBaseDecl.hpp

#ifndef TSFCORE_MPI_MULTI_VECTOR_BASE_DECL_HPP
#define TSFCORE_MPI_MULTI_VECTOR_BASE_DECL_HPP

#include <vector>

#include "TSFCoreMultiVector.hpp"
#include "Teuchos_BLAS.hpp"

namespace TSFCore {

///
template<class Scalar> class MPIVectorSpaceBase;

///
/** Base class for MPI-based multi-vectors.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class MPIMultiVectorBase : virtual public MultiVector<Scalar> {
public:

	///
	using MultiVector<Scalar>::apply; // Inject *all* the apply methods!

	///
	MPIMultiVectorBase();

	/** @name Pure virtual methods to be overridden by subclasses */
	//@{

	///
	/** Returns the MPI-based vector space object for the range of <tt>*this</tt> multi-vectorr.
	 */
	virtual Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> > mpiSpace() const = 0;

	///
	/** Returns a <tt>const</tt>  pointer to a fortran-style view of the local multi-vector data.
	 *
	 * @param  values      [out] On output <tt>*values</tt> will point to 
	 *                     the first element in the first colum of the local multi-vector
	 *                     stored as a column-major dense Fortran-style matrix.
	 * @param  leadingDim  [out] On output <tt>*leadingDim</tt> gives the leading dimension
	 *                     of the Fortran-style local multi-vector.
	 *
	 */
	virtual void getLocalData( const Scalar **values, Index *leadingDim ) const = 0;

	///
	/** Free view of local data that was gotten from <tt>getLocalData()</tt>.
	 *
	 * @param  values      [in/out] On input <tt>values</tt> must be the pointer set
	 *                     by <tt>getLocalData()</tt>.
	 */
	virtual void freeLocalData( const Scalar *values ) const = 0;

	///
	/** Returns a non-<tt>const</tt> pointer to a fortran-style view of the local multi-vector data.
	 *
	 * @param  values      [out] On output <tt>*values</tt> will point to 
	 *                     the first element in the first colum of the local multi-vector
	 *                     stored as a column-major dense Fortran-style matrix.
	 * @param  leadingDim  [out] On output <tt>*leadingDim</tt> gives the leading dimension
	 *                     of the Fortran-style local multi-vector.
	 *
	 * The function <tT>commitLocalData()</tt> must be called to
	 * commit changes to the data.
	 */
	virtual void getLocalData( Scalar **values, Index *leadingDim ) = 0;

	///
	/** Commit view of local data that was gotten from <tt>getLocalData()</tt>.
	 *
	 * @param  values      [in/out] On input <tt>*values</tt> must be the pointer set
	 *                     by <tt>getLocalData()</tt>.
	 */
	virtual void commitLocalData( Scalar *values ) = 0;

	//@}

	/** @name Virtual methods with default implementaions */
	//@{


	//@}

	/** @name Overridden from OpBase */
	//@{

	/// Returns <tt>mpiSpace</tt>.
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > range() const;

	//@}

	/** @name Overridden from LinearOp */
	//@{

	///
	/** Wrapps the <tt>Vector</tt> objects in <tt>MultiVector</tt> objects then calls
	 * the <tt>MultiVector</tt> version of <tt>apply()</tt>
	 */
	void apply(
		const ETransp            M_trans
		,const Vector<Scalar>    &x
		,Vector<Scalar>          *y
		,const Scalar            alpha
		,const Scalar            beta
		) const;

	///
	/** Uses GEMM(...) and MPI_Allreduce(...) to implement.
	 *
	 * ToDo: Finish documentation!
	 */
	void apply(
		const ETransp                 M_trans
		,const MultiVector<Scalar>    &X
		,MultiVector<Scalar>          *Y
		,const Scalar                 alpha
		,const Scalar                 beta
		) const;

	//@}

	/** @name Overridden from MultiVector */
	//@{
	///
	void applyOp(
		const RTOpPack::RTOpT<Scalar>   &primary_op
		,const size_t                   num_multi_vecs
		,const MultiVector<Scalar>*     multi_vecs[]
		,const size_t                   num_targ_multi_vecs
		,MultiVector<Scalar>*           targ_multi_vecs[]
		,RTOp_ReductTarget              reduct_objs[]
		,const Index                    primary_first_ele
		,const Index                    primary_sub_dim
		,const Index                    primary_global_offset
		,const Index                    secondary_first_ele
		,const Index                    secondary_sub_dim
		) const;
	///
	void getSubMultiVector(
		const Range1D                       &rowRng
		,const Range1D                      &colRng
		,RTOpPack::SubMultiVectorT<Scalar>  *sub_mv
		) const;
	///
	void freeSubMultiVector( RTOpPack::SubMultiVectorT<Scalar>* sub_mv ) const;
	///
	void getSubMultiVector(
		const Range1D                                &rowRng
		,const Range1D                               &colRng
		,RTOpPack::MutableSubMultiVectorT<Scalar>    *sub_mv
		);
	///
	void commitSubMultiVector( RTOpPack::MutableSubMultiVectorT<Scalar>* sub_mv );
	//@}

protected:

	///
	/** Subclasses should call whenever the structure of the VectorSpace changes.
	 *
	 * This function can be overridden by subclasses but this
	 * particualar function implementation must be called from within
	 * any override.
	 */
	virtual void updateMpiSpace();

	///
	/** Validate and resize the row range.
	 *
	 * This function throws an exception if the input range is invalid
	 */
	Range1D validateRowRange( const Range1D& rowRng ) const;


	///
	/** Validate and resize the column range.
	 *
	 * This function throws an exception if the input range is invalid
	 */
	Range1D validateColRange( const Range1D& rowCol ) const;
	
private:
	
	// ///////////////////////////////////////
	// Private data members
	
	mutable bool in_applyOp_;

	mutable Teuchos::BLAS<int,Scalar> blas_;

	// cached
	Index  globalDim_;
	Index  localOffset_;
	Index  localSubDim_;
	Index  numCols_;

	
}; // end class MPIMultiVectorBase

} // end namespace TSFCore

#endif // TSFCORE_MPI_MULTI_VECTOR_BASE_DECL_HPP
