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
// TSFCoreAssertOp.hpp

#ifndef TSFCORE_ASSERT_OP_HPP
#define TSFCORE_ASSERT_OP_HPP

#include "TSFCoreAssertOpDecl.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreOpBase.hpp"
#include "Teuchos_TestForException.hpp"

// boilerplate code

namespace TSFCore {

template<class Scalar>
struct dump_vec_spaces_t {
public:
	dump_vec_spaces_t(
		const TSFCore::VectorSpace<Scalar>& _vec_space1, const std::string &_vec_space1_name
		,const TSFCore::VectorSpace<Scalar>& _vec_space2, const std::string &_vec_space2_name
		)
		:vec_space1(_vec_space1),vec_space1_name(_vec_space1_name)
		,vec_space2(_vec_space2),vec_space2_name(_vec_space2_name)
		{}
	const TSFCore::VectorSpace<Scalar> &vec_space1;
	const std::string                  *vec_space1_name;
	const TSFCore::VectorSpace<Scalar> &vec_space2;
	const std::string                  *vec_space2_name;
}; // end dum_vec_spaces

template<class Scalar>
inline dump_vec_spaces_t<Scalar> dump_vec_spaces(
	const TSFCore::VectorSpace<Scalar>& vec_space1, const std::string &vec_space1_name
	,const TSFCore::VectorSpace<Scalar>& vec_space2, const std::string &vec_space2_name
	)
{
	return dump_vec_spaces(vec_space1,vec_space1_name,vec_space2,vec_space2_name);
}

// Notice!!!!!!!  Place a breakpoint in following function in order to halt the
// program just before an exception is thrown!

template<class Scalar>
std::ostream& operator<<( std::ostream& o, const dump_vec_spaces_t<Scalar>& d )
{
	o << "Error, " << d.vec_space1_name << " at address " << &d.vec_space1
	  << " of type \'" << typeid(d.vec_space1).name()
	  << "\' with dimension " << d.vec_space1_name << ".dim() = " << d.vec_space1.dim()
	  << " is not compatible with "
	  << d.vec_space2_name  << " at address " << &d.vec_space2
	  << " of type \'" << typeid(d.vec_space2).name()
	  << "\' with dimension " << d.vec_space2_name << ".dim() = " << d.vec_space2.dim();
	return o;
}

enum EM_VS { VS_RANGE, VS_DOMAIN };

template<class Scalar>
const TSFCore::VectorSpace<Scalar>& op(
	const TSFCore::OpBase<Scalar>&     M
	,TSFCore::ETransp            M_trans
	,EM_VS                    M_VS
	)
{
	using TSFCore::NOTRANS;
	using TSFCore::TRANS;
	if(M_trans == NOTRANS && M_VS == VS_RANGE)
		return *M.range();
	if(M_trans == TRANS && M_VS == VS_RANGE)
		return *M.domain();
	if(M_trans == NOTRANS && M_VS == VS_DOMAIN)
		return *M.domain();
	// M_trans == TRANS && M_VS == VS_DOMAIN
	return *M.range();
}

} // end namespace TSFCore

#define TSFCORE_ASSERT_LHS_ARG(FUNC_NAME,LHS_ARG) \
	TEST_FOR_EXCEPTION( \
		(LHS_ARG) == NULL, std::invalid_argument \
		,FUNC_NAME << " : Error!" \
		);

// Notice!!!!!!!  Setting a breakpoint a the line number that is printed by this macro
// and then trying to set the condition !isCompatible does not work (at least not
// in gdb).

#define TSFCORE_ASSERT_VEC_SPACES_NAMES(FUNC_NAME,VS1,VS1_NAME,VS2,VS2_NAME) \
{ \
	const bool isCompatible = (VS1).isCompatible(VS2); \
	TEST_FOR_EXCEPTION( \
		!isCompatible, ::TSFCore::Exceptions::IncompatibleVectorSpaces \
		,FUNC_NAME << " : "	<< ::TSFCore::dump_vec_spaces(VS1,VS1_NAME,VS2,VS2_NAME) \
		) \
}

#define TSFCORE_ASSERT_VEC_SPACES(FUNC_NAME,VS1,VS2)\
TSFCORE_ASSERT_VEC_SPACES_NAMES(FUNC_NAME,VS1,#VS1,VS2,#VS2)

#define TSFCORE_ASSERT_MAT_VEC_SPACES(FUNC_NAME,M,M_T,M_VS,VS) \
{ \
	std::ostringstream M_VS_name; \
	M_VS_name << "(" #M << ( (M_T) == NOTRANS ? "" : "'" ) << ")" \
			   << "." << ( (M_VS) == VS_RANGE ? "range()" : "domain()" ); \
	TSFCORE_ASSERT_VEC_SPACES_NAMES( \
		FUNC_NAME \
		,op(M,M_T,M_VS),M_VS_name.str().c_str() \
		,(VS),#VS \
		) \
}

#define TSFCORE_ASSERT_MAT_MAT_SPACES(FUNC_NAME,M1,M1_T,M1_VS,M2,M2_T,M2_VS) \
{ \
	std::ostringstream M1_VS_name, M2_VS_name; \
	M1_VS_name << "(" #M1 << ( M1_T == NOTRANS ? "" : "'" ) << ")" \
			   << "." << ( M1_VS == VS_RANGE ? "range()" : "domain()" ); \
	M2_VS_name << "(" #M2 << ( M2_T == NOTRANS ? "" : "'" ) << ")" \
			   << "." << ( M2_VS == VS_RANGE ? "range()" : "domain()" ); \
	TSFCORE_ASSERT_VEC_SPACES_NAMES( \
		FUNC_NAME \
		,op(M1,M1_T,M1_VS),M1_VS_name.str().c_str() \
		,op(M2,M2_T,M2_VS),M2_VS_name.str().c_str() \
		) \
}

// function definitions

#ifdef TSFCORE_ASSERT_COMPATIBILITY

template<class Scalar>
void TSFCore::Vp_V_assert_compatibility(Vector<Scalar>* v_lhs, const Vector<Scalar>& v_rhs)
{
	const char func_name[] = "Vp_V_assert_compatibility(v_lhs,v_rhs)";
	TSFCORE_ASSERT_LHS_ARG(func_name,v_lhs)
	TSFCORE_ASSERT_VEC_SPACES("Vp_V_assert_compatibility(v_lhs,v_rhs)",*v_lhs->space(),*v_rhs.space());
}

template<class Scalar>
void TSFCore::VopV_assert_compatibility(const Vector<Scalar>& v_rhs1, const Vector<Scalar>&  v_rhs2)
{
	const char func_name[] = "VopV_assert_compatibility(v_rhs1,v_rhs2)";
	TSFCORE_ASSERT_VEC_SPACES(func_name,*v_rhs1.space(),*v_rhs2.space());
}

template<class Scalar>
void TSFCore::Vp_MtV_assert_compatibility(
	Vector<Scalar>* v_lhs
	,const OpBase<Scalar>& m_rhs1, ETransp trans_rhs1, const Vector<Scalar>& v_rhs2 )
{
	const char func_name[] = "Vp_MtV_assert_compatibility(v_lhs,m_rhs1,trans_rhs1,v_rhs2)";
	TSFCORE_ASSERT_LHS_ARG(func_name,v_lhs)
	TSFCORE_ASSERT_MAT_VEC_SPACES(func_name,m_rhs1,trans_rhs1,VS_RANGE,*v_lhs->space())
	TSFCORE_ASSERT_MAT_VEC_SPACES(func_name,m_rhs1,trans_rhs1,VS_DOMAIN,*v_rhs2.space())
}

template<class Scalar>
void TSFCore::Mp_MtM_assert_compatibility(
	OpBase<Scalar>* m_lhs, ETransp trans_lhs
	,const OpBase<Scalar>& m_rhs1, ETransp trans_rhs1
	,const OpBase<Scalar>& m_rhs2, ETransp trans_rhs2 )
{
	const char func_name[] = "Mp_MtM_assert_compatibility(m_lhs,trans_lhsm_rhs1,trans_rhs1,m_rhs2,trans_rhs2)";
	TSFCORE_ASSERT_LHS_ARG(func_name,m_lhs)
	TSFCORE_ASSERT_MAT_MAT_SPACES(func_name,(*m_lhs),trans_lhs,VS_RANGE,m_rhs1,trans_rhs1,VS_RANGE)
	TSFCORE_ASSERT_MAT_MAT_SPACES(func_name,(*m_lhs),trans_lhs,VS_DOMAIN,m_rhs2,trans_rhs2,VS_DOMAIN)
	TSFCORE_ASSERT_MAT_MAT_SPACES(func_name,m_rhs1,trans_rhs1,VS_DOMAIN,m_rhs2,trans_rhs2,VS_RANGE)
}

#endif // TSFCORE_ASSERT_COMPATIBILITY

#endif // TSFCORE_ASSERT_OP_HPP
