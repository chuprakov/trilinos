// //////////////////////////////////////////////////////////////
// TSFCoreExplicitVectorView.hpp

#include "TSFCoreVector.hpp"

#ifndef TSFCORE_EXPLICIT_VECTOR_VIEW_HPP
#define TSFCORE_EXPLICIT_VECTOR_VIEW_HPP

namespace TSFCore {

///
/** Create an explicit non-mutable (const) view of a <tt>Vector</tt> object.
 */
template<class Scalar>
class ExplicitVectorView {
public:
	///
	ExplicitVectorView( const Vector<Scalar>& v ) : v_(v) { v_.getSubVector(Range1D(),&sv_); }
	///
	~ExplicitVectorView() { v_.freeSubVector(&sv_); }
	///
	RTOp_index_type   globalOffset() const { return sv_.globalOffset(); }
	///
	RTOp_index_type   subDim()       const { return sv_.subDim();  }
	///
	const Scalar*     values()       const { return sv_.values();  }
	///
	ptrdiff_t         stride()       const { return sv_.stride();  }
	/// Preconditions: <tt>values()!=NULL && (1 <= i <= subDim())</tt>
	const Scalar& operator()(RTOp_index_type i) const { return sv_(i); }
private:
	const Vector<Scalar>          &v_;
	RTOpPack::SubVectorT<Scalar>  sv_;
	// Not defined and not to be called
	ExplicitVectorView();
	ExplicitVectorView(const ExplicitVectorView<Scalar>&);
	ExplicitVectorView<Scalar>& operator==(const ExplicitVectorView<Scalar>&);
};
 
///
/** Create an explicit mutable (non-const) view of a <tt>Vector</tt> object.
 */
template<class Scalar>
class ExplicitMutableVectorView {
public:
	///
	ExplicitMutableVectorView( Vector<Scalar>& v ) : v_(v) { v_.getSubVector(Range1D(),&sv_); }
	///
	~ExplicitMutableVectorView() { v_.commitSubVector(&sv_); }
	///
	RTOp_index_type   globalOffset() const { return sv_.globalOffset(); }
	///
	RTOp_index_type   subDim()       const { return sv_.subDim();  }
	///
	Scalar*           values()       const { return sv_.values();  }
	///
	ptrdiff_t         stride()       const { return sv_.stride();  }
	/// Preconditions: <tt>values()!=NULL && (1 <= i <= subDim())</tt>
	Scalar& operator()(RTOp_index_type i) const { return sv_(i); }
private:
	Vector<Scalar>                       &v_;
	RTOpPack::MutableSubVectorT<Scalar>  sv_;
	// Not defined and not to be called
	ExplicitMutableVectorView();
	ExplicitMutableVectorView(const ExplicitMutableVectorView<Scalar>&);
	ExplicitMutableVectorView<Scalar>& operator==(const ExplicitMutableVectorView<Scalar>&);
};

} // namespace TSFCore

#endif // TSFCORE_EXPLICIT_VECTOR_VIEW_HPP
