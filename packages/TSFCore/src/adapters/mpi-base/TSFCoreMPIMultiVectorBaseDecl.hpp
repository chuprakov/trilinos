// /////////////////////////////////////////////////////////////////////////////
// TSFCoreMPIMultiVectorBaseDecl.hpp

#ifndef TSFCORE_MPI_MULTI_VECTOR_BASE_DECL_HPP
#define TSFCORE_MPI_MULTI_VECTOR_BASE_DECL_HPP

#include <vector>

#include "TSFCoreMultiVector.hpp"

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
	MPIMultiVectorBase();

	/** @name Pure virtual methods to be overridden by subclasses */
	//@{

	///
	/** Returns the MPI-based vector space object for the range of <tt>*this</tt> multi-vectorr.
	 */
	virtual Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> > mpiSpace() const = 0;

	///
	/** Returns a pointer to the local multi-vector data.
	 *
	 * @param  values      [out] On output <tt>*values</tt> will point to 
	 *                     the first element in the first colum of the local multi-vector
	 *                     stored as a column-major dense Fortran-style matrix.
	 * @param  leadingDim  [out] On output <tt>*leadingDim</tt> gives the leading dimension
	 *                     of the Fortran-style local multi-vector.
	 */
	virtual void getLocalData( Scalar **values, Index *leadingDim ) = 0;

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
	
private:
	
	// ///////////////////////////////////////
	// Private data members
	
	mutable bool in_applyOp_;
	
}; // end class MPIMultiVectorBase

} // end namespace TSFCore

#endif // TSFCORE_MPI_MULTI_VECTOR_BASE_DECL_HPP
