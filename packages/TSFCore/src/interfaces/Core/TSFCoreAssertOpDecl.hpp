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

// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreAssertOpDecl.hpp

#ifndef TSFCORE_ASSERT_OP_DECL_HPP
#define TSFCORE_ASSERT_OP_DECL_HPP

#include "TSFCoreTypes.hpp"

namespace TSFCore {

/** \defgroup AssertOp_funcs Assertion functions for linear algebra operations.
  *
  * These functions check the compatibility of the vector spaces for many different types
  * of linear algebra operations and throw <tt>Exceptions::IncompatibleVectorSpaces</tt> 
  * expressions if the vector spaces do not match.  These functions only perform there checks
  * if <tt>TSFCORE_ASSERT_COMPATIBILITY</tt> is defined.  These functions will also throw
  * <tt>std::invalid_argument</tt> if a lhs argument is <tt>NULL</tt>.
  */
//@{

/// v_lhs += op v_rhs
template<class Scalar>
void Vp_V_assert_compatibility(Vector<Scalar>* v_lhs, const Vector<Scalar>& v_rhs);
/// v_rhs1 op v_rhs2
template<class Scalar>
void VopV_assert_compatibility(const Vector<Scalar>& v_rhs1, const Vector<Scalar>& v_rhs2);
/// op(m_lhs) += op op(m_rhs)
template<class Scalar>
void Mp_M_assert_compatibility(
	OpBase<Scalar>* m_lhs, ETransp trans_lhs
	,const OpBase<Scalar>& m_rhs, ETransp trans_rhs );
/// v_lhs += op(m_rhs1) * v_rhs2
template<class Scalar>
void Vp_MtV_assert_compatibility(
	Vector<Scalar>* v_lhs
	,const OpBase<Scalar>& m_rhs1, ETransp trans_rhs1, const Vector<Scalar>& v_rhs2 );
/// op(m_lhs) += op(m_rhs1) * op(m_rhs2)
template<class Scalar>
void Mp_MtM_assert_compatibility(
	OpBase<Scalar>* m_lhs, ETransp trans_lhs
	,const OpBase<Scalar>& m_rhs1, ETransp trans_rhs1
	,const OpBase<Scalar>& m_rhs2, ETransp trans_rhs2 );

//@}

#ifdef _DEBUG
#define TSFCORE_ASSERT_COMPATIBILITY
#endif

#ifndef TSFCORE_ASSERT_COMPATIBILITY

// inline definitions that do nothing

template<class Scalar>
inline
void Vp_V_assert_compatibility(Vector<Scalar>* v_lhs, const Vector<Scalar>& v_rhs)
{} 
template<class Scalar>
inline
void VopV_assert_compatibility(const Vector<Scalar>& v_rhs1, const Vector<Scalar>&  v_rhs2)
{}
template<class Scalar>
inline
void Vp_MtV_assert_compatibility(
	Vector<Scalar>* v_lhs
	,const OpBase<Scalar>& m_rhs1, ETransp trans_rhs1, const Vector<Scalar>& v_rhs2 )
{}
template<class Scalar>
inline
void Mp_MtM_assert_compatibility(
	OpBase<Scalar>* m_lhs, ETransp trans_lhs
	,const OpBase<Scalar>& m_rhs1, ETransp trans_rhs1
	,const OpBase<Scalar>& m_rhs2, ETransp trans_rhs2 )
{}

#endif // TSFCORE_ASSERT_COMPATIBILITY

} // end namespace TSFCore

#endif	// TSFCORE_ASSERT_OP_DECL_HPP
