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
