// /////////////////////////////////////////////////////////////////////
// VectorDenseEncap.h

#ifndef SLAP_VECTOR_DENSE_ENCAP_H
#define SLAP_VECTOR_DENSE_ENCAP_H

#include "SparseLinAlgPackTypes.h"
#include "AbstractLinAlgPack/src/VectorWithOpMutable.h"
#include "LinAlgPack/src/VectorClass.h"

namespace SparseLinAlgPack {

///
/** Extract a constant <tt>LinAlgPack::VectorSlice</tt> view of a <tt>VectorWithOp</tt> object.
 *
 * This class is only to be used with <tt>VectotrWithOp</tt> objects that store all of their
 * elements in the local address space or can easily access all of the vector elements in
 * this process (or thread).  It is generally to be used in serial applications but might
 * also find use in parallel appliations where a vector is replicated across processes.
 *
 * This utility class is defined purly in terms of the abstract interfaces.  It is only to
 * be used as an automatic variable on the stack.  For example, to extract a <tt>VectorSlice</tt>
 * view of an abstract vector and use it to copy to another <tt>VectorSlice</tt> object you could
 * write a function like:
 \code
 void copy(const VectorWithOp& vec_in, VectorSlice* vs_out ) {
     VectorDenseEncap  vs_in(vec_in);
	 *vs_out = vs_in();
 }
 \endcode
 * In the above code, if the underlying <tt>VectorWithOp</tt> object does not have to
 * perform any dynamic memory allocations and copy in the method <tt>VectorWithOp::get_sub_vector()</tt>
 * then the above code will only have a constant time overhead.  However, the above approach will work
 * for any <tt>VectorWithOp</tt> object (no matter how inefficient it may be).
 */
class VectorDenseEncap {
public:

	/// Calls <tt>vec.get_sub_vector(Range1D(),DENSE,&sub_vec)</tt> to get the view.  
	VectorDenseEncap( const VectorWithOp&  vec );
	/// Calls <tt>vec.free_sub_vector(&sub_vec)</tt> to release the view.  
	~VectorDenseEncap();
	/// Returns a reference to a constant view of the dense vector.
	const VectorSlice& operator()() const;

private:

	const VectorWithOp &vec_;
	RTOp_SubVector     sub_vec_;
	VectorSlice        vs_;
	VectorDenseEncap();                                     // Not defined and not to be called!
	VectorDenseEncap(const VectorDenseEncap&);              // ""
	VectorDenseEncap& operator=(const VectorDenseEncap&);   // ""

}; // end class VectorDenseEncap

///
/** Extract a non-const <tt>LinAlgPack::VectorSlice</tt> view of a <tt>VectorWithOpMutable</tt> object.
 *
 * This utility class is defined purly in terms of the abstract interfaces.  It is only to
 * be used as an automatic variable on the stack.  Note that the underlying <tt>VectorWithOpMutable</tt>
 * object is  not guarrenteed to be modified until the destructor for \c this is called.
 */
class VectorDenseMutableEncap {
public:

	/// Calls <tt>vec.get_sub_vector(Range1D(),&sub_vec)</tt> to get the view.  
	VectorDenseMutableEncap( VectorWithOpMutable&  vec );
	/// Calls <tt>vec.commit_sub_vector(&sub_vec)</tt> to release the view.  
	~VectorDenseMutableEncap();
	/// Returns a reference to a constant view of the dense vector.
	VectorSlice& operator()();
	/// Returns a reference to a non-const view of the dense vector.
	const VectorSlice& operator()() const;

private:

	VectorWithOpMutable       &vec_;
	RTOp_MutableSubVector     sub_vec_;
	VectorSlice               vs_;
	VectorDenseMutableEncap();                                            // Not defined and not to be called!
	VectorDenseMutableEncap(const VectorDenseMutableEncap&);              // ""
	VectorDenseMutableEncap& operator=(const VectorDenseMutableEncap&);   // ""

}; // end class VectorDenseMutableEncap

// ///////////////////////////////////////////
// Inline members

// VectorDenseEncap

inline
VectorDenseEncap::VectorDenseEncap( const VectorWithOp&  vec )
	:vec_(vec)
{
	RTOp_sub_vector_null(&sub_vec_);
	vec_.get_sub_vector(Range1D(),VectorWithOp::DENSE,&sub_vec_);
	vs_.bind( VectorSlice(
				  const_cast<value_type*>(sub_vec_.values)
				  ,sub_vec_.sub_dim
				  ,sub_vec_.values_stride
				  )
		);
}

inline
VectorDenseEncap::~VectorDenseEncap()
{
	vec_.free_sub_vector(&sub_vec_);
}

inline
const VectorSlice& VectorDenseEncap::operator()() const
{
	return vs_;
}

// VectorDenseMutableEncap

inline
VectorDenseMutableEncap::VectorDenseMutableEncap( VectorWithOpMutable&  vec )
	:vec_(vec)
{
	RTOp_mutable_sub_vector_null(&sub_vec_);
	vec_.get_sub_vector(Range1D(),&sub_vec_);
	vs_.bind( VectorSlice(
				  sub_vec_.values
				  ,sub_vec_.sub_dim
				  ,sub_vec_.values_stride
				  )
		);
}

inline
VectorDenseMutableEncap::~VectorDenseMutableEncap()
{
	vec_.commit_sub_vector(&sub_vec_);
}

inline
VectorSlice& VectorDenseMutableEncap::operator()()
{
	return vs_;
}

inline
const VectorSlice& VectorDenseMutableEncap::operator()() const
{
	return vs_;
}

} // end namespace SparseLinALgPack

#endif // SLAP_VECTOR_DENSE_ENCAP_H
