// /////////////////////////////////////////////////////
// RTOpPackTypes.hpp

#ifndef RTOPPACK_TYPES_HPP
#define RTOPPACK_TYPES_HPP

#include <stdexcept>

#include "RTOp_config.h"

namespace RTOpPack {

//
// Exceptions
//

///
class UnknownError : public std::logic_error
{public: UnknownError(const std::string& what_arg) : std::logic_error(what_arg) {}};
///
class InvalidUsage : public std::logic_error
{public: InvalidUsage(const std::string& what_arg) : std::logic_error(what_arg) {}};
///
class InvalidNumVecs : public std::logic_error
{public: InvalidNumVecs(const std::string& what_arg) : std::logic_error(what_arg) {}};
///
class InvalidNumTargVecs : public std::logic_error
{public: InvalidNumTargVecs(const std::string& what_arg) : std::logic_error(what_arg) {}};
///
class IncompatibleVecs : public std::logic_error
{public: IncompatibleVecs(const std::string& what_arg) : std::logic_error(what_arg) {}};

//
// Vector subviews
//

///
/** Class for a non-mutable sub-vector.
 *
 * For a sub-vector <tt>vec</tt>, the corresponding entries
 *	in the global vector <tt>x(j)</tt> (one based) are as follows:
 \verbatim

	x( vec.globalOffset() + k ) = v.(k), for k = 1,...,vec.subDim()
  \endverbatim
 * The stride <tt>vec.stride()</tt> may be positive (>0), negative (<0)
 * or even zero (0).  A negative stride <tt>vec.stride() < 0</tt> allows a
 * reverse traversal of the elements.  A zero stride
 * <tt>vec.stride()</tt> allows a sub-vector with all the elements the same.
 *
 * The raw pointer to the start of the memory can be obtained as
 * <tt>&vec(1)</tt>.
 *
 * Warning! the default copy constructor and assignement operators are
 * allowed which results in only pointer copy, not deep copy!  You
 * have been warned!
 */
template<class Scalar>
class MutableSubVectorT {
public:
	///
	MutableSubVectorT() : globalOffset_(0), subDim_(0), values_(NULL), stride_(0) {}
	///
	MutableSubVectorT(RTOp_index_type globalOffset, RTOp_index_type subDim, Scalar *values, ptrdiff_t stride)
		:globalOffset_(globalOffset), subDim_(subDim), values_(values), stride_(stride) 
		{}
	///
	void initialize(RTOp_index_type globalOffset, RTOp_index_type subDim, Scalar *values, ptrdiff_t stride)
		{ globalOffset_=globalOffset; subDim_=subDim; values_=values; stride_=stride;  }
	///
	void set_uninitialized()
		{ globalOffset_ = 0; subDim_=0; values_=NULL; stride_ = 0; }
	///
	void setGlobalOffset(RTOp_index_type globalOffset) { globalOffset_ = globalOffset; } 
	///
	RTOp_index_type   globalOffset() const { return globalOffset_; }
	///
	RTOp_index_type   subDim()       const { return subDim_;  }
	///
	Scalar*           values()       const { return values_;  }
	///
	ptrdiff_t         stride()       const { return stride_;  }
	/// Preconditions: <tt>values()!=NULL && (1 <= i <= subDim())</tt>
	Scalar& operator()(RTOp_index_type i) const { return values_[ stride_*(i-1) ]; }
private:
	RTOp_index_type     globalOffset_;
	RTOp_index_type     subDim_;
	Scalar              *values_;
	ptrdiff_t           stride_;
};

///
/** Class for a non-mutable sub-vector.
 *
 * For a sub-vector <tt>vec</tt>, the corresponding entries
 *	in the global vector <tt>x(j)</tt> (one based) are as follows:
 \verbatim

	x( vec.globalOffset() + k ) = v.(k), for k = 1,...,vec.subDim()
  \endverbatim
 * The stride <tt>vec.stride()</tt> may be positive (>0), negative (<0)
 * or even zero (0).  A negative stride <tt>vec.stride() < 0</tt> allows a
 * reverse traversal of the elements.  A zero stride
 * <tt>vec.stride()</tt> allows a sub-vector with all the elements the same.
 *
 * The raw pointer to the start of the memory can be obtained as
 * <tt>&vec(1)</tt>.
 *
 * Warning! the default copy constructor and assignement operators are
 * allowed which results in only pointer copy, not deep copy!  You
 * have been warned!
 */
template<class Scalar>
class SubVectorT {
public:
	///
	SubVectorT() : globalOffset_(0), subDim_(0), values_(NULL), stride_(0) {}
	///
	SubVectorT(RTOp_index_type globalOffset, RTOp_index_type subDim, const Scalar *values, ptrdiff_t stride)
		:globalOffset_(globalOffset), subDim_(subDim), values_(values), stride_(stride) 
		{}
	///
	SubVectorT( const MutableSubVectorT<Scalar>& sv )
		:globalOffset_(sv.globalOffset()), subDim_(sv.subDim()), values_(sv.values()), stride_(sv.stride()) 
		{}
	///
	void initialize(RTOp_index_type globalOffset, RTOp_index_type subDim, const Scalar *values, ptrdiff_t stride)
		{ globalOffset_=globalOffset; subDim_=subDim; values_=values; stride_=stride;  }
	///
	void set_uninitialized()
		{ globalOffset_ = 0; subDim_=0; values_=NULL; stride_ = 0; }
	///
	void setGlobalOffset(RTOp_index_type globalOffset) { globalOffset_ = globalOffset; } 
	///
	RTOp_index_type   globalOffset() const { return globalOffset_; }
	///
	RTOp_index_type   subDim()       const { return subDim_;  }
	///
	const Scalar*     values()       const { return values_;  }
	///
	ptrdiff_t         stride()       const { return stride_;  }
	/// Preconditions: <tt>values()!=NULL && (1 <= i <= subDim())</tt>
	const Scalar& operator()(RTOp_index_type i) const { return values_[ stride_*(i-1) ]; }
private:
	RTOp_index_type     globalOffset_;
	RTOp_index_type     subDim_;
	const Scalar        *values_;
	ptrdiff_t           stride_;
};

///
/** Class for a (sparse or dense) sub-vector.
 *
 *	Sparse and dense local vectors are supported as follows:
  *
  *	A dense vector <tt>vec</tt> is identified by <tt>vec.subDim() == vec.sub_nz</tt>
  * and <tt>vec.indices() == NULL</tt> in which case
  *	<tt>vec.indicesStride()</tt>, <tt>vec.localOffset()</tt> and <tt>vec.isSorted()</tt>
  * are ignored.  For a dense sub-vector <tt>vec</tt>, the corresponding entries
 *	in the global vector <tt>x(j)</tt> (one based) are as follows:
 \verbatim

	x( vec.globalOffset() + k )
		= vec.values()[ vec.valueStride() * (k-1) ]

	for k = 1,...,vec.subDim()
 \endverbatim
 * The stride member <tt>vec.valueStride()()</tt> may be positive (>0), negative (<0)
 * or even zero (0).  A negative stride <tt>vec.valueStride() < 0</tt> allows a
 * reverse traversal of the elements in <tt>vec.values()[]</tt>.  A zero stride
 * <tt>vec.valueStride()() == 0</tt> allows a vector with all the elements the same.
 *
 *	A sparse vector is identified by <tt>vec.subDim() > vec.sub_nz()</tt>
 * or <tt>vec.indices() != NULL</tt>
 * in which case all the fields in the structure are meaningful.
 *	The corresponding elements in the global vector <tt>x(j)</tt>
 * defined as:
 \verbatim

	x( vec.globalOffset() + vec.localOffset() + vec.indices()[vec.indicesStride()*(k-1)] )
		= vec.values[vec.valueStride()*(k-1)]

	for k = 1,...,vec.sub_nz
 \endverbatim
 * If <tt>vec.sub_nz == 0</tt> then it is allowed for <tt>vec.indices() == NULL</tt>.
 * If <tt>vec.subDim() > vec.sub_nz > 0</tt> then <tt>vec.indices() != NULL</tt> must be true.
 *
  * A sparse sub-vector may be sorted (<tt>vec.isSorted()!=0</tt>) or
  * unsorted (<tt>vec.isSorted()==0</tt>) but the indices <tt>vec.indices()[k]</tt>
  * must be unique.  A sorted vector (<tt>vec.isSorted()!=0</tt>) means that
  * the indices are in ascending order:
  \verbatim

	vec.indices()[vec.indicesStride()*(k-1)] < vec.indices()[vec.indicesStride()*(k)]

	for k = 1,...,vec.sub_nz-1
 \endverbatim
 * The member <tt>vec.localOffset()</tt> is used to shift the values in <tt>vec.indices()[]</tt>
 * to be in range of the local sub-vector.  In other words:
 \verbatim
	
	1 <= vec.localOffset() + vec.indices()[vec.indicesStride()*(k-1)] <= vec.sub_nz

	for k = 1...vec.sub_nz
 \endverbatim
 * The member <tt>vec.valueStride()</tt> may be positive (>0), negative (<0) or zero (0).
 * However, the member <tt>vec.indicesStride()</tt> may be may be positive (>0)
 * or negative (<0) but not zero (0).  Allowing <tt>vec.indicesStride() == 0</tt>
 * would mean that a vector would have <tt>vec.sub_nz</tt> nonzero elements with
 * all the same value and all the same indexes and non-unique indices are
 * not allowed.  Allowing non-unique indexes would make some operations
 * (e.g. dot product) very difficult to implement and therefore can not
 * be allowed.  A sparse vector where <tt>vec.valueStride() == 0</tt> is one
 * where all of the nonzeros have the value <tt>vec.values[0]</tt>.  If
 * <tt>vec.sub_nz == 0</tt> for a sparse vector then it is allowed for
 * <tt>vec.values == NULL</tt> and <tt>vec.indices() == NULL</tt>.
 *
 *	This specification allows a lot of flexibility in determining
 * how the vectors are laid out in memory.  However, allowing vectors to be
 * sparse and unsorted may make many user defined operations
 * considerably harder and expensive to implement.
 *
 * To avoid making mistakes in setting the members of this struct use
 * one of the helper functions <tt>RTOp_sparse_sub_vector_from_dense()</tt>,
 * <tt>RTOp_sparse_sub_vector()</tt> or <tt>RTOp_sub_vector_null()</tt>.
 */
template<class Scalar>
class SparseSubVectorT {
public:
	///
	SparseSubVectorT()
		:globalOffset_(0),subDim_(0),subNz_(0)
		,values_(NULL),valuesStride_(0),indices_(NULL)
		,indicesStride_(0),localOffset_(0),isSorted_(0)
		{}
	///
	SparseSubVectorT(
		RTOp_index_type globalOffset, RTOp_index_type subDim, RTOp_index_type subNz
		,const Scalar values[], ptrdiff_t valuesStride
		,const RTOp_index_type indices[], ptrdiff_t indicesStride
		,ptrdiff_t localOffset, int isSorted
		)
		:globalOffset_(globalOffset),subDim_(subDim),subNz_(subNz)
		,values_(values),valuesStride_(valuesStride),indices_(indices)
		,indicesStride_(indicesStride),localOffset_(localOffset),isSorted_(isSorted)
		{}
	///
	SparseSubVectorT(
		RTOp_index_type globalOffset, RTOp_index_type subDim
		,const Scalar values[], ptrdiff_t valuesStride
		)
		:globalOffset_(globalOffset),subDim_(subDim),subNz_(subDim)
		,values_(values),valuesStride_(valuesStride),indices_(NULL)
		,indicesStride_(0),localOffset_(0),isSorted_(0)
		{}
	///
	SparseSubVectorT( const SubVectorT<Scalar>& sv )
		:globalOffset_(sv.globalOffset()), subDim_(sv.subDim()), subNz_(sv.subDim()), values_(sv.values())
		,valuesStride_(sv.stride()), indices_(NULL),indicesStride_(0),localOffset_(0),isSorted_(0)
		{}
	///
	void initialize(
		RTOp_index_type globalOffset, RTOp_index_type subDim, RTOp_index_type subNz
		,const Scalar values[], ptrdiff_t valuesStride
		,const RTOp_index_type indices[], ptrdiff_t indicesStride
		,ptrdiff_t localOffset, int isSorted
		)
		{
			globalOffset_ = globalOffset; subDim_ = subDim; subNz_ = subNz;
			values_ = values; valuesStride_ = valuesStride; indices_ = indices;
			indicesStride_ = indicesStride; localOffset_ = localOffset; isSorted_ = isSorted;
		}
	///
	void initialize(
		RTOp_index_type globalOffset, RTOp_index_type subDim
		,const Scalar values[], ptrdiff_t valuesStride
		)
		{
			globalOffset_ = globalOffset; subDim_ = subDim; subNz_ = subDim;
			values_ = values; valuesStride_ = valuesStride; indices_ = NULL;
			indicesStride_ = 0; localOffset_ = 0; isSorted_ = 1;
		}
	///
	void set_uninitialized()
		{
			globalOffset_ = 0; subDim_ = 0; subNz_ = 0;
			values_ = NULL; valuesStride_ = 0; indices_ = NULL;
			indicesStride_ = 0; localOffset_ = 0; isSorted_ = 1;
		}
	///
	void setGlobalOffset(RTOp_index_type globalOffset) { globalOffset_ = globalOffset; } 
	/// Offset for the sub-vector into the global vector
	RTOp_index_type                  globalOffset() const { return globalOffset_; } 
	/// Dimension of the sub-vector
	RTOp_index_type                  subDim() const { return subDim_; }
	/// Number of nonzero elements (<tt>subNz == subDim</tt> for dense vectors)
	RTOp_index_type                  subNz() const { return subNz_; }
	/// Array (size min{|<tt>valueStride*subNz</tt>|,1}) for the values in the vector
	const Scalar*                    values() const { return values_; }
	/// Stride between elements in <tt>values[]</tt>
	ptrdiff_t                        valuesStride() const { return valuesStride_; }
	///
	/** Array (size min{|<tt>indicesStride*subNz</tt>|,1} if not <tt>NULL</tt>) for the
	  * indices of the nonzero elements in the vector (sparse vectors only)
	  */
	const RTOp_index_type*           indices() const { return indices_; }
	/// Stride between indices in indices[] (sparse vectors only)
	ptrdiff_t                        indicesStride() const { return indicesStride_; }
	/// Offset of indices[] into local sub-vector (sparse vectors only)
	ptrdiff_t                        localOffset() const { return localOffset_; }
	/// If <tt>isSorted == 0</tt> then the vector is not sorted, otherwise it is sorted (sparse vectors only)
	int                              isSorted() const { return isSorted_; }
private:
	RTOp_index_type                  globalOffset_;
	RTOp_index_type                  subDim_;
	RTOp_index_type                  subNz_;
	const Scalar                     *values_;
	ptrdiff_t                        valuesStride_;
	const RTOp_index_type            *indices_;
	ptrdiff_t                        indicesStride_;
	ptrdiff_t                        localOffset_;
	int                              isSorted_;
};

//
// Templated types
//

///
template<class Scalar>  class MutableSubVectorT;
///
template<class Scalar>  class SubVectorT;
///
template<class Scalar>  class SparseSubVectorT;
///
template<class Scalar>  class ReductTargetT;
///
template<class Scalar>  class RTOpT;

//
// Typedefs
//

///
typedef MutableSubVectorT<RTOp_value_type>  MutableSubVector;
///
typedef SubVectorT<RTOp_value_type>         SubVector;
///
typedef SparseSubVectorT<RTOp_value_type>   SparseSubVector;
///
typedef ReductTargetT<RTOp_value_type>      ReductTarget;
///
typedef RTOpT<RTOp_value_type>              RTOp;

} // namespace RTOpPack

#endif // RTOPPACK_TYPES_HPP

