// //////////////////////////////////////////////////////////////
// TSFCoreExplicitMultiVectorView.hpp

#include "TSFCoreMultiVector.hpp"

#ifndef TSFCORE_EXPLICIT_MULTI_VECTOR_VIEW_HPP
#define TSFCORE_EXPLICIT_MULTI_VECTOR_VIEW_HPP

namespace TSFCore {

///
/** Create an explicit non-mutable (const) view of a <tt>MultiVector</tt> object.
 */
template<class Scalar>
class ExplicitMultiVectorView {
public:
	///
	ExplicitMultiVectorView(
		const MultiVector<Scalar>& mv, const Range1D &rowRng = Range1D(), const Range1D &colRng = Range1D()
		)
		: mv_(mv) { mv_.getSubMultiVector(rowRng,colRng,&smv_); }
	///
	~ExplicitMultiVectorView() { mv_.freeSubMultiVector(&smv_); }
	///
	const RTOpPack::SubMultiVectorT<Scalar>& smv() const { return smv_; }
	///
	RTOp_index_type   globalOffset() const { return smv_.globalOffset(); }
	///
	RTOp_index_type   subDim()       const { return smv_.subDim();  }
	///
	RTOp_index_type   colOffset()    const { return smv_.colOffset(); }
	///
	RTOp_index_type   numSubCols()   const { return smv_.numSubCols();  }
	///
	const Scalar*     values()       const { return smv_.values();  }
	///
	RTOp_index_type   leadingDim()   const { return smv_.leadingDim();  }
	/// One-based indexing: Preconditions: <tt>values()!=NULL && (1<=i<=subDim()) && (1<=j<=numSubCols())</tt>
	const Scalar& operator()(RTOp_index_type i,RTOp_index_type j) const { return smv_(i,j); }
private:
	const MultiVector<Scalar>          &mv_;
	RTOpPack::SubMultiVectorT<Scalar>  smv_;
	// Not defined and not to be called
	ExplicitMultiVectorView();
	ExplicitMultiVectorView(const ExplicitMultiVectorView<Scalar>&);
	ExplicitMultiVectorView<Scalar>& operator==(const ExplicitMultiVectorView<Scalar>&);
};
 
///
/** Create an explicit mutable (non-const) view of a <tt>MultiVector</tt> object.
 */
template<class Scalar>
class ExplicitMutableMultiVectorView {
public:
	///
	ExplicitMutableMultiVectorView(
		MultiVector<Scalar>& mv, const Range1D &rowRng = Range1D(), const Range1D &colRng = Range1D()
		)
		: mv_(mv) { mv_.getSubMultiVector(rowRng,colRng,&smv_); }
	///
	~ExplicitMutableMultiVectorView() { mv_.commitSubMultiVector(&smv_); }
	///
	const RTOpPack::MutableSubMultiVectorT<Scalar>& smv() const { return smv_; }
	///
	RTOp_index_type   globalOffset() const { return smv_.globalOffset(); }
	///
	RTOp_index_type   subDim()       const { return smv_.subDim();  }
	///
	RTOp_index_type   colOffset()    const { return smv_.colOffset(); }
	///
	RTOp_index_type   numSubCols()   const { return smv_.numSubCols();  }
	///
	Scalar*           values()       const { return smv_.values();  }
	///
	RTOp_index_type   leadingDim()   const { return smv_.leadingDim();  }
	/// One-based indexing: Preconditions: <tt>values()!=NULL && (1<=i<=subDim()) && (1<=j<=numSubCols())</tt>
	Scalar& operator()(RTOp_index_type i,RTOp_index_type j) { return smv_(i,j); }
private:
	MultiVector<Scalar>                       &mv_;
	RTOpPack::MutableSubMultiVectorT<Scalar>  smv_;
	// Not defined and not to be called
	ExplicitMutableMultiVectorView();
	ExplicitMutableMultiVectorView(const ExplicitMutableMultiVectorView<Scalar>&);
	ExplicitMutableMultiVectorView<Scalar>& operator==(const ExplicitMutableMultiVectorView<Scalar>&);
};

} // namespace TSFCore

#endif // TSFCORE_EXPLICIT_MULTI_VECTOR_VIEW_HPP
