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

// //////////////////////////////////////////////////////////////
// TSFCoreExplicitMultiVectorView.hpp

#include "TSFCoreMultiVector.hpp"

#ifndef TSFCORE_EXPLICIT_MULTI_VECTOR_VIEW_HPP
#define TSFCORE_EXPLICIT_MULTI_VECTOR_VIEW_HPP

namespace TSFCore {

///
/** Create an explicit non-mutable (const) view of a <tt>MultiVector</tt> object.
 *
 * \ingroup TSFCore_ANA_Development_grp
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
 *
 * \ingroup TSFCore_ANA_Development_grp
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
