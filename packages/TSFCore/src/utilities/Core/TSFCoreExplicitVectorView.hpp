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
	ExplicitVectorView( const Vector<Scalar>& v, const Range1D &rng = Range1D() ) : v_(v) { v_.getSubVector(rng,&sv_); }
	///
	~ExplicitVectorView() { v_.freeSubVector(&sv_); }
	///
	const RTOpPack::SubVectorT<Scalar>& sv() const { return sv_; }
	///
	RTOp_index_type   globalOffset() const { return sv_.globalOffset(); }
	///
	RTOp_index_type   subDim()       const { return sv_.subDim();  }
	///
	const Scalar*     values()       const { return sv_.values();  }
	///
	ptrdiff_t         stride()       const { return sv_.stride();  }
	/// Zero-based indexing: Preconditions: <tt>values()!=NULL && (0 <= i <= subDim()-1)</tt>
	const Scalar& operator[](RTOp_index_type i) const { return sv_[i]; }
	/// One-based indexing: Preconditions: <tt>values()!=NULL && (1 <= i <= subDim())</tt>
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
	ExplicitMutableVectorView( Vector<Scalar>& v, const Range1D &rng = Range1D() ) : v_(v) { v_.getSubVector(rng,&sv_); }
	///
	~ExplicitMutableVectorView() { v_.commitSubVector(&sv_); }
	///
	const RTOpPack::MutableSubVectorT<Scalar>& sv() const { return sv_; }
	///
	RTOp_index_type   globalOffset() const { return sv_.globalOffset(); }
	///
	RTOp_index_type   subDim()       const { return sv_.subDim();  }
	///
	Scalar*           values()       const { return sv_.values();  }
	///
	ptrdiff_t         stride()       const { return sv_.stride();  }
	/// Zero-based indexing: Preconditions: <tt>values()!=NULL && (0 <= i <= subDim()-1)</tt>
	Scalar& operator[](RTOp_index_type i) const { return sv_[i]; }
	/// One-based indexing: Preconditions: <tt>values()!=NULL && (1 <= i <= subDim())</tt>
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
