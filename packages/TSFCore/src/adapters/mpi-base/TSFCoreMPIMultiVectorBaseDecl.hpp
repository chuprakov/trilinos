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
	virtual MemMngPack::ref_count_ptr<const MPIVectorSpaceBase<Scalar> > mpiSpace() const = 0;
	//@}

	/** @name Overridden from OpBase */
	//@{
	/// Returns <tt>mpiSpace</tt>.
	MemMngPack::ref_count_ptr< const VectorSpace<Scalar> > range() const;
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
	//@}
	
private:
	
	// ///////////////////////////////////////
	// Private data members
	
	mutable bool in_applyOp_;
	
}; // end class MPIMultiVectorBase

} // end namespace TSFCore

#endif // TSFCORE_MPI_MULTI_VECTOR_BASE_DECL_HPP
