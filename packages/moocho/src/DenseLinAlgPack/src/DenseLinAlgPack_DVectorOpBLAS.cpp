// /////////////////////////////////////////////////////////////////////////////////////////
// VectorOpBLAS.cpp

#include <algorithm>
#include <functional>
#include <math.h> // VC++ 5.0 <cmath> is not CD2 complient yet
	// ToDo: Update math function calls to cmath once you get a compiler that meets the standard.

#include "../include/VectorClass.h"
#include "../include/VectorOp.h"
#include "../include/BLAS_Cpp.h"
#include "../include/LinAlgPackAssertOp.h"

// ToDo: Check into the aliasing problem more completely and take care of partial overlap as
// this is not concidered yet.

// File scope utilites
namespace {

using LinAlgPack::Vector;
using LinAlgPack::VectorSlice;
typedef LinAlgPack::value_type value_type;
typedef VectorSlice::size_type size_type;
typedef VectorSlice::difference_type difference_type;

// Utility for copying vector slices.  Takes care of aliasing etc. but not sizes.
void i_assign(VectorSlice* vs_lhs, const VectorSlice& vs_rhs) {
	switch(vs_lhs->overlap(vs_rhs)) {
		case LinAlgPack::SAME_MEM:	return; // assignment to self so skip the copy
		case LinAlgPack::SOME_OVERLAP:		// create a temp for the copy
		{
			Vector temp(vs_rhs);
			BLAS_Cpp::copy( temp.size(), temp.raw_ptr(), 1, vs_lhs->raw_ptr(), vs_lhs->stride() );
			return;
		}
		default:				// no overlap so just perform the copy
		{
			BLAS_Cpp::copy( vs_rhs.size(),vs_rhs.raw_ptr(), vs_rhs.stride()
				, vs_lhs->raw_ptr(), vs_lhs->stride() );
			return;
		}
	}
}

}	// end namespace

// ////////////////////////////////////////////////////////////////////////////////////////////////
// op= operations

// vs_lhs += alpha
void LinAlgPack::Vp_S(VectorSlice* vs_lhs, value_type alpha) {
	if(vs_lhs->stride() == 1) {
		VectorSlice::value_type
			*itr		= vs_lhs->start_ptr(),
			*itr_end	= itr + vs_lhs->size();
		while(itr != itr_end)
			*itr++ += alpha;
	}
	else {
		VectorSlice::iterator
			itr		= vs_lhs->begin(),
			itr_end	= vs_lhs->end();
		while(itr != itr_end)
			*itr++ += alpha;
	}
}

// vs_lhs *= alpha (BLAS xSCAL)
void LinAlgPack::Vt_S(VectorSlice* vs_lhs, value_type alpha) {
	BLAS_Cpp::scal( vs_lhs->size(), alpha, vs_lhs->raw_ptr(), vs_lhs->stride() );
}

// vs_lhs += alpha * vs_rhs (BLAS xAXPY)
void LinAlgPack::Vp_StV(VectorSlice* vs_lhs, value_type alpha, const VectorSlice& vs_rhs) {
	Vp_V_assert_sizes(vs_lhs->size(), vs_rhs.size());
	BLAS_Cpp::axpy( vs_lhs->size(), alpha, vs_rhs.raw_ptr(), vs_rhs.stride()
		, vs_lhs->raw_ptr(), vs_lhs->stride());
}

// ///////////////////////////////////////////////////////////////////////////////////////////////
// Vector as lhs

// v_lhs = alpha
void LinAlgPack::assign(Vector* v_lhs, value_type alpha)
{
	if(!v_lhs->size())
		throw std::length_error("LinAlgPack::assign(v_lhs,alpha): Vector must be sized.");
	else
		BLAS_Cpp::copy( v_lhs->size(), &alpha, 0, v_lhs->raw_ptr(), v_lhs->stride() );
}
// v_lhs = vs_rhs
void LinAlgPack::assign(Vector* v_lhs, const VectorSlice& vs_rhs) {
	v_lhs->resize(vs_rhs.size());
	i_assign( &(*v_lhs)(), vs_rhs );
}
// v_lhs = vs_rhs1 + vs_rhs2
void LinAlgPack::V_VpV(Vector* v_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2)
{
	assert_vs_sizes(vs_rhs1, vs_rhs2);
	v_lhs->resize(vs_rhs1.size());
	std::transform(vs_rhs1.begin(),vs_rhs1.end(),vs_rhs2.begin(),v_lhs->begin(),std::plus<value_type>());
//	VectorSlice::const_iterator
//		v1 = vs_rhs1.begin(),
//		v2 = vs_rhs2.begin();
//	Vector::iterator
//		v = v_lhs->begin(),
//		v_end = v_lhs->end();
//	while( v != v_end )
//		*v++ = *v1++ + *v2++;
}
// v_lhs = vs_rhs1 - vs_rhs2
void LinAlgPack::V_VmV(Vector* v_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2)
{
	assert_vs_sizes(vs_rhs1, vs_rhs2);
	v_lhs->resize(vs_rhs1.size());
	std::transform(vs_rhs1.begin(),vs_rhs1.end(),vs_rhs2.begin(),v_lhs->begin(),std::minus<value_type>());
}
// v_lhs = - vs_rhs
void LinAlgPack::V_mV(Vector* v_lhs, const VectorSlice& vs_rhs) {
	v_lhs->resize(vs_rhs.size());
	(*v_lhs) = vs_rhs;
	BLAS_Cpp::scal(v_lhs->size(), -1.0, v_lhs->raw_ptr(), 1);
}
// v_lhs = alpha * vs_rhs
void LinAlgPack::V_StV(Vector* v_lhs, value_type alpha, const VectorSlice& vs_rhs) {
	v_lhs->resize(vs_rhs.size());
	(*v_lhs) = vs_rhs;
	BLAS_Cpp::scal( v_lhs->size(), alpha, v_lhs->raw_ptr(), 1);
}

// Elementwise math vector functions. VC++ 5.0 is not allowing use of ptr_fun() with overloaded
// functions so I have to perform the loops straight out.
// ToDo: use ptr_fun() when you get a compiler that works this out.  For now I will just use macros.
#define UNARYOP_VEC(LHS, RHS, FUNC)																			\
	LHS->resize(RHS.size());																				\
	Vector::iterator itr_lhs; VectorSlice::const_iterator itr_rhs;											\
	for(itr_lhs = LHS->begin(), itr_rhs = RHS.begin(); itr_lhs != LHS->end(); ++itr_lhs, ++itr_rhs)			\
	{	*itr_lhs = FUNC(*itr_rhs); }

#define BINARYOP_VEC(LHS, RHS1, RHS2, FUNC)																	\
	LinAlgPack::assert_vs_sizes(RHS1, RHS2); LHS->resize(RHS1.size());										\
	Vector::iterator itr_lhs; VectorSlice::const_iterator itr_rhs1, itr_rhs2;								\
	for(itr_lhs = LHS->begin(), itr_rhs1 = RHS1.begin(), itr_rhs2 = RHS2.begin();							\
		itr_lhs != LHS->end(); ++itr_lhs, ++itr_rhs1, ++itr_rhs2)											\
	{	*itr_lhs = FUNC(*itr_rhs1, *itr_rhs2); }

#define BINARYOP_BIND1ST_VEC(LHS, RHS1, RHS2, FUNC)															\
	LHS->resize(RHS2.size());																				\
	Vector::iterator itr_lhs; VectorSlice::const_iterator itr_rhs2;											\
	for(itr_lhs = LHS->begin(), itr_rhs2 = RHS2.begin();													\
		itr_lhs != LHS->end(); ++itr_lhs, ++itr_rhs2)														\
	{	*itr_lhs = FUNC(RHS1, *itr_rhs2); }

#define BINARYOP_BIND2ND_VEC(LHS, RHS1, RHS2, FUNC)															\
	LHS->resize(RHS1.size());																				\
	Vector::iterator itr_lhs; VectorSlice::const_iterator itr_rhs1;											\
	for(itr_lhs = LHS->begin(), itr_rhs1 = RHS1.begin();													\
		itr_lhs != LHS->end(); ++itr_lhs, ++itr_rhs1)														\
	{	*itr_lhs = FUNC(*itr_rhs1, RHS2); }

// v_lhs = abs(vs_rhs)
void LinAlgPack::abs(Vector* v_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VEC(v_lhs, vs_rhs, ::fabs)
}
// v_lhs = asin(vs_rhs)
void LinAlgPack::asin(Vector* v_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VEC(v_lhs, vs_rhs, ::asin)
}
// v_lhs = acos(vs_rhs)
void LinAlgPack::acos(Vector* v_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VEC(v_lhs, vs_rhs, ::acos)
}
// v_lhs = atan(vs_rhs)
void LinAlgPack::atan(Vector* v_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VEC(v_lhs, vs_rhs, ::atan)
}
// v_lhs = atan(vs_rhs1/vs_rhs2)
void LinAlgPack::atan2(Vector* v_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2)
{
	BINARYOP_VEC(v_lhs, vs_rhs1, vs_rhs2, ::atan2)
}
// v_lhs = atan(vs_rhs/alpha)
void LinAlgPack::atan2(Vector* v_lhs, const VectorSlice& vs_rhs, value_type alpha) {
	BINARYOP_BIND2ND_VEC(v_lhs, vs_rhs, alpha, ::atan2)
}
// v_lhs = atan(alpha/vs_rhs)
void LinAlgPack::atan2(Vector* v_lhs, value_type alpha, const VectorSlice& vs_rhs) {
	BINARYOP_BIND1ST_VEC(v_lhs, alpha, vs_rhs, ::atan2)
}
// v_lhs = cos(vs_rhs)
void LinAlgPack::cos(Vector* v_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VEC(v_lhs, vs_rhs, ::cos)
}
// v_lhs = cosh(vs_rhs)
void LinAlgPack::cosh(Vector* v_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VEC(v_lhs, vs_rhs, ::cosh)
}
// v_lhs = exp(vs_rhs)
void LinAlgPack::exp(Vector* v_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VEC(v_lhs, vs_rhs, ::exp)
}
// v_lhs = max(vs_rhs1,vs_rhs2)
void LinAlgPack::max(Vector* v_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2)
{
	LinAlgPack::assert_vs_sizes(vs_rhs1, vs_rhs2); v_lhs->resize(vs_rhs1.size());
	Vector::iterator itr_lhs; VectorSlice::const_iterator itr_rhs1, itr_rhs2;
	for(itr_lhs = v_lhs->begin(), itr_rhs1 = vs_rhs1.begin(), itr_rhs2 = vs_rhs2.begin();
		itr_lhs != v_lhs->end(); ++itr_lhs, ++itr_rhs1, ++itr_rhs2)
	{	*itr_lhs = std::_MAX(*itr_rhs1, *itr_rhs2); }
}
// v_lhs = max(alpha,vs_rhs)
void LinAlgPack::max(Vector* v_lhs, value_type alpha, const VectorSlice& vs_rhs) {
	v_lhs->resize(vs_rhs.size());
	Vector::iterator itr_lhs; VectorSlice::const_iterator itr_rhs;
	for(itr_lhs = v_lhs->begin(), itr_rhs = vs_rhs.begin(); itr_lhs != v_lhs->end(); ++itr_lhs, ++itr_rhs)
	{	*itr_lhs = std::_MAX(alpha,*itr_rhs); }
}
// v_lhs = min(vs_rhs1,vs_rhs2)
void LinAlgPack::min(Vector* v_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2) 
{
	LinAlgPack::assert_vs_sizes(vs_rhs1, vs_rhs2); v_lhs->resize(vs_rhs1.size());
	Vector::iterator itr_lhs; VectorSlice::const_iterator itr_rhs1, itr_rhs2;
	for(itr_lhs = v_lhs->begin(), itr_rhs1 = vs_rhs1.begin(), itr_rhs2 = vs_rhs2.begin();
		itr_lhs != v_lhs->end(); ++itr_lhs, ++itr_rhs1, ++itr_rhs2)
	{	*itr_lhs = std::_MIN(*itr_rhs1, *itr_rhs2); }
}
// v_lhs = min(alpha,vs_rhs)
void LinAlgPack::min(Vector* v_lhs, value_type alpha, const VectorSlice& vs_rhs) {
	v_lhs->resize(vs_rhs.size());
	Vector::iterator itr_lhs; VectorSlice::const_iterator itr_rhs;
	for(itr_lhs = v_lhs->begin(), itr_rhs = vs_rhs.begin(); itr_lhs != v_lhs->end(); ++itr_lhs, ++itr_rhs)
	{	*itr_lhs = std::_MIN(alpha,*itr_rhs); }
}
// v_lhs = pow(vs_rhs1,vs_rhs2)
void LinAlgPack::pow(Vector* v_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2)
{
	BINARYOP_VEC(v_lhs, vs_rhs1, vs_rhs2, ::pow)
}
// v_lhs = pow(vs_rhs,alpha)
void LinAlgPack::pow(Vector* v_lhs, const VectorSlice& vs_rhs, value_type alpha) {
	BINARYOP_BIND2ND_VEC(v_lhs, vs_rhs, alpha, ::pow)
}
// v_lhs = pow(vs_rhs,n) 
void LinAlgPack::pow(Vector* v_lhs, const VectorSlice& vs_rhs, int n) {
	BINARYOP_BIND2ND_VEC(v_lhs, vs_rhs, n, ::pow)
}
// v_lhs = pow(alpha,vs_rhs)
void LinAlgPack::pow(Vector* v_lhs, value_type alpha, const VectorSlice& vs_rhs) {
	BINARYOP_BIND1ST_VEC(v_lhs, alpha, vs_rhs, ::pow)
}
// v_lhs = sqrt(vs_rhs)
void LinAlgPack::sqrt(Vector* v_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VEC(v_lhs, vs_rhs, ::sqrt)
}
// v_lhs = sin(vs_rhs)
void LinAlgPack::sin(Vector* v_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VEC(v_lhs, vs_rhs, ::sin)
}
// v_lhs = sinh(vs_rhs)
void LinAlgPack::sinh(Vector* v_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VEC(v_lhs, vs_rhs, ::sinh)
}
// v_lhs = tan(vs_rhs)
void LinAlgPack::tan(Vector* v_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VEC(v_lhs, vs_rhs, ::tan)
}
// v_lhs = tanh(vs_rhs)
void LinAlgPack::tanh(Vector* v_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VEC(v_lhs, vs_rhs, ::tanh)
}

// ///////////////////////////////////////////////////////////////////////////////////////////////
// VectorSlice as lhs

// vs_lhs = alpha
void LinAlgPack::assign(VectorSlice* vs_lhs, value_type alpha) {
	if(!vs_lhs->size())
		throw std::length_error("LinAlgPack::assign(vs_lhs,alpha): VectorSlice must be bound and sized.");
	BLAS_Cpp::copy( vs_lhs->size(), &alpha, 0, vs_lhs->raw_ptr(), vs_lhs->stride() );
}
// vs_lhs = vs_rhs
void LinAlgPack::assign(VectorSlice* vs_lhs, const VectorSlice& vs_rhs) {
	Vp_V_assert_sizes( vs_lhs->size(), vs_rhs.size() );
	i_assign(vs_lhs,vs_rhs);
}
// vs_lhs = vs_rhs1 + vs_rhs2
void LinAlgPack::V_VpV(VectorSlice* vs_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2) {
	VopV_assert_sizes( vs_rhs1.size(), vs_rhs2.size() );
	Vp_V_assert_sizes( vs_lhs->size(), vs_rhs1.size() );
	std::transform(vs_rhs1.begin(),vs_rhs1.end(),vs_rhs2.begin(),vs_lhs->begin(),std::plus<value_type>());
}
// vs_lhs = vs_rhs1 - vs_rhs2
void LinAlgPack::V_VmV(VectorSlice* vs_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2)
{
	VopV_assert_sizes( vs_rhs1.size(), vs_rhs2.size() );
	Vp_V_assert_sizes( vs_lhs->size(), vs_rhs1.size() );
	std::transform(vs_rhs1.begin(),vs_rhs1.end(),vs_rhs2.begin(),vs_lhs->begin(),std::minus<value_type>());
}
// vs_lhs = - vs_rhs
void LinAlgPack::V_mV(VectorSlice* vs_lhs, const VectorSlice& vs_rhs)
{
	Vp_V_assert_sizes( vs_lhs->size(), vs_rhs.size() );
	(*vs_lhs) = vs_rhs;
	BLAS_Cpp::scal( vs_lhs->size(), -1.0, vs_lhs->raw_ptr(), vs_lhs->stride() );
}
// vs_lhs = alpha * vs_rhs
void LinAlgPack::V_StV(VectorSlice* vs_lhs, value_type alpha, const VectorSlice& vs_rhs)
{
	Vp_V_assert_sizes( vs_lhs->size(), vs_rhs.size() );
	(*vs_lhs) = vs_rhs;
	BLAS_Cpp::scal( vs_lhs->size(), alpha, vs_lhs->raw_ptr(), vs_lhs->stride() );
}

// Elementwise math vector functions. VC++ 5.0 is not allowing use of ptr_fun() with overloaded
// functions so I have to perform the loops straight out.
// ToDo: use ptr_fun() when you get a compiler that works this out.  For now I will just use macros.
#define UNARYOP_VECSLC(LHS, RHS, FUNC)																		\
	LinAlgPack::Vp_V_assert_sizes( (LHS)->size(), (RHS).size() );											\
	VectorSlice::iterator itr_lhs; VectorSlice::const_iterator itr_rhs;										\
	for(itr_lhs = LHS->begin(), itr_rhs = RHS.begin(); itr_lhs != LHS->end(); ++itr_lhs, ++itr_rhs)			\
	{	*itr_lhs = FUNC(*itr_rhs); }


#define BINARYOP_VECSLC(LHS, RHS1, RHS2, FUNC)																\
	LinAlgPack::VopV_assert_sizes( (RHS1).size(), (RHS2).size() );											\
	LinAlgPack::Vp_V_assert_sizes( (LHS)->size(), (RHS1).size() );											\
	VectorSlice::iterator itr_lhs; VectorSlice::const_iterator itr_rhs1, itr_rhs2;							\
	for(itr_lhs = LHS->begin(), itr_rhs1 = RHS1.begin(), itr_rhs2 = RHS2.begin();							\
		itr_lhs != LHS->end(); ++itr_lhs, ++itr_rhs1, ++itr_rhs2)											\
	{	*itr_lhs = FUNC(*itr_rhs1, *itr_rhs2); }

#define BINARYOP_BIND1ST_VECSLC(LHS, RHS1, RHS2, FUNC)														\
	LinAlgPack::Vp_V_assert_sizes( (LHS)->size(), (RHS2).size() );											\
	VectorSlice::iterator itr_lhs; VectorSlice::const_iterator itr_rhs2;									\
	for(itr_lhs = LHS->begin(), itr_rhs2 = RHS2.begin();													\
		itr_lhs != LHS->end(); ++itr_lhs, ++itr_rhs2)														\
	{	*itr_lhs = FUNC(RHS1, *itr_rhs2); }

#define BINARYOP_BIND2ND_VECSLC(LHS, RHS1, RHS2, FUNC)														\
	LinAlgPack::Vp_V_assert_sizes( (LHS)->size(), (RHS1).size());											\
	VectorSlice::iterator itr_lhs; VectorSlice::const_iterator itr_rhs1;									\
	for(itr_lhs = LHS->begin(), itr_rhs1 = RHS1.begin();													\
		itr_lhs != LHS->end(); ++itr_lhs, ++itr_rhs1)														\
	{	*itr_lhs = FUNC(*itr_rhs1, RHS2); }

// vs_lhs = abs(vs_rhs)
void LinAlgPack::abs(VectorSlice* vs_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VECSLC(vs_lhs, vs_rhs, ::fabs)
}
// vs_lhs = asin(vs_rhs)
void LinAlgPack::asin(VectorSlice* vs_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VECSLC(vs_lhs, vs_rhs, ::asin)
}
// vs_lhs = acos(vs_rhs)
void LinAlgPack::acos(VectorSlice* vs_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VECSLC(vs_lhs, vs_rhs, ::acos)
}
// vs_lhs = atan(vs_rhs)
void LinAlgPack::atan(VectorSlice* vs_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VECSLC(vs_lhs, vs_rhs, ::atan)
}
// vs_lhs = atan(vs_rhs1/vs_rhs2)
void LinAlgPack::atan2(VectorSlice* vs_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2)
{
	BINARYOP_VECSLC(vs_lhs, vs_rhs1, vs_rhs2, ::atan2)
}
// vs_lhs = atan(vs_rhs/alpha)
void LinAlgPack::atan2(VectorSlice* vs_lhs, const VectorSlice& vs_rhs, value_type alpha)
{
	BINARYOP_BIND2ND_VECSLC(vs_lhs, vs_rhs, alpha, ::atan2)
}
// vs_lhs = atan(alpha/vs_rhs)
void LinAlgPack::atan2(VectorSlice* vs_lhs, value_type alpha, const VectorSlice& vs_rhs)
{
	BINARYOP_BIND1ST_VECSLC(vs_lhs, alpha, vs_rhs, ::atan2)
}
// vs_lhs = cos(vs_rhs)
void LinAlgPack::cos(VectorSlice* vs_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VECSLC(vs_lhs, vs_rhs, ::cos)
}
// vs_lhs = cosh(vs_rhs)
void LinAlgPack::cosh(VectorSlice* vs_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VECSLC(vs_lhs, vs_rhs, ::cosh)
}
// vs_lhs = exp(vs_rhs)
void LinAlgPack::exp(VectorSlice* vs_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VECSLC(vs_lhs, vs_rhs, ::exp)
}
// vs_lhs = max(vs_rhs1,vs_rhs2)
void LinAlgPack::max(VectorSlice* vs_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2)
{
	LinAlgPack::VopV_assert_sizes( vs_rhs1.size(), vs_rhs2.size() );
	LinAlgPack::Vp_V_assert_sizes( vs_lhs->size(), vs_rhs1.size() );
	VectorSlice::iterator itr_lhs; VectorSlice::const_iterator itr_rhs1, itr_rhs2;
	for(itr_lhs = vs_lhs->begin(), itr_rhs1 = vs_rhs1.begin(), itr_rhs2 = vs_rhs2.begin();
		itr_lhs != vs_lhs->end(); ++itr_lhs, ++itr_rhs1, ++itr_rhs2)
	{	*itr_lhs = std::_MAX(*itr_rhs1, *itr_rhs2); }
}
// vs_lhs = max(alpha,vs_rhs)
void LinAlgPack::max(VectorSlice* vs_lhs, value_type alpha, const VectorSlice& vs_rhs)
{
	LinAlgPack::Vp_V_assert_sizes( vs_lhs->size(), vs_rhs.size() );
	VectorSlice::iterator itr_lhs; VectorSlice::const_iterator itr_rhs;
	for(itr_lhs = vs_lhs->begin(), itr_rhs = vs_rhs.begin(); itr_lhs != vs_lhs->end(); ++itr_lhs, ++itr_rhs)
	{	*itr_lhs = std::_MAX(alpha,*itr_rhs); }
}
// vs_lhs = min(vs_rhs1,vs_rhs2)
void LinAlgPack::min(VectorSlice* vs_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2)
{
	LinAlgPack::VopV_assert_sizes( vs_rhs1.size(), vs_rhs2.size() );
	LinAlgPack::Vp_V_assert_sizes( vs_lhs->size(), vs_rhs1.size() );
	VectorSlice::iterator itr_lhs; VectorSlice::const_iterator itr_rhs1, itr_rhs2;
	for(itr_lhs = vs_lhs->begin(), itr_rhs1 = vs_rhs1.begin(), itr_rhs2 = vs_rhs2.begin();
		itr_lhs != vs_lhs->end(); ++itr_lhs, ++itr_rhs1, ++itr_rhs2)
	{	*itr_lhs = std::_MIN(*itr_rhs1, *itr_rhs2); }
}
// vs_lhs = min(alpha,vs_rhs)
void LinAlgPack::min(VectorSlice* vs_lhs, value_type alpha, const VectorSlice& vs_rhs)
{
	LinAlgPack::Vp_V_assert_sizes( vs_lhs->size(), vs_rhs.size() );
	VectorSlice::iterator itr_lhs; VectorSlice::const_iterator itr_rhs;
	for(itr_lhs = vs_lhs->begin(), itr_rhs = vs_rhs.begin(); itr_lhs != vs_lhs->end(); ++itr_lhs, ++itr_rhs)
	{	*itr_lhs = std::_MIN(alpha,*itr_rhs); }
}
// vs_lhs = pow(vs_rhs1,vs_rhs2)
void LinAlgPack::pow(VectorSlice* vs_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2) 
{
	BINARYOP_VECSLC(vs_lhs, vs_rhs1, vs_rhs2, ::pow)
}
// vs_lhs = pow(vs_rhs,alpha)
void LinAlgPack::pow(VectorSlice* vs_lhs, const VectorSlice& vs_rhs, value_type alpha)
{
	BINARYOP_BIND2ND_VECSLC(vs_lhs, vs_rhs, alpha, ::pow)
}
// vs_lhs = pow(vs_rhs,n) 
void LinAlgPack::pow(VectorSlice* vs_lhs, const VectorSlice& vs_rhs, int n) {
	BINARYOP_BIND2ND_VECSLC(vs_lhs, vs_rhs, n, ::pow)
}
// vs_lhs = pow(alpha,vs_rhs)
void LinAlgPack::pow(VectorSlice* vs_lhs, value_type alpha, const VectorSlice& vs_rhs)
{
	BINARYOP_BIND1ST_VECSLC(vs_lhs, alpha, vs_rhs, ::pow)
}
// vs_lhs = sqrt(vs_rhs)
void LinAlgPack::sqrt(VectorSlice* vs_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VECSLC(vs_lhs, vs_rhs, ::sqrt)
}
// vs_lhs = sin(vs_rhs)
void LinAlgPack::sin(VectorSlice* vs_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VECSLC(vs_lhs, vs_rhs, ::sin)
}
// vs_lhs = sinh(vs_rhs)
void LinAlgPack::sinh(VectorSlice* vs_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VECSLC(vs_lhs, vs_rhs, ::sinh)
}
// vs_lhs = tan(vs_rhs)
void LinAlgPack::tan(VectorSlice* vs_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VECSLC(vs_lhs, vs_rhs, ::tan)
}
// vs_lhs = tanh(vs_rhs)
void LinAlgPack::tanh(VectorSlice* vs_lhs, const VectorSlice& vs_rhs) {
	UNARYOP_VECSLC(vs_lhs, vs_rhs, ::tanh)
}

// ////////////////////////////////////////////////////////////////////////////////////////
// Scalar returning functions

// result = trans(vs_rhs1) * vs_rhs2
LinAlgPack::value_type LinAlgPack::dot(const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2)
{
	VopV_assert_sizes( vs_rhs1.size(), vs_rhs2.size() );
	return BLAS_Cpp::dot( vs_rhs1.size(), vs_rhs1.raw_ptr(), vs_rhs1.stride()
		, vs_rhs2.raw_ptr(), vs_rhs2.stride() );
}
// result = max(vs_rhs)
LinAlgPack::value_type LinAlgPack::max(const VectorSlice& vs_rhs)
{	return *std::max_element(vs_rhs.begin(),vs_rhs.end()); }
// result = min(vs_rhs)
LinAlgPack::value_type LinAlgPack::min(const VectorSlice& vs_rhs)
{	return *std::min_element(vs_rhs.begin(),vs_rhs.end()); }
// result = ||vs_rhs||1
LinAlgPack::value_type LinAlgPack::norm_1(const VectorSlice& vs_rhs)
{	return BLAS_Cpp::asum( vs_rhs.size(), vs_rhs.raw_ptr(), vs_rhs.stride() ); }
// result = ||vs_rhs||2
LinAlgPack::value_type LinAlgPack::norm_2(const VectorSlice& vs_rhs) 
{	return BLAS_Cpp::nrm2( vs_rhs.size(), vs_rhs.raw_ptr(), vs_rhs.stride() ); }
// result = ||vs_rhs||infinity
LinAlgPack::value_type LinAlgPack::norm_inf(const VectorSlice& vs_rhs) {
//	return BLAS_Cpp::iamax( vs_rhs.size(), vs_rhs.raw_ptr(), vs_rhs.stride() );
//	For some reason iamax() is not working properly?
	value_type max_ele = 0.0;
	for(VectorSlice::const_iterator itr = vs_rhs.begin(); itr != vs_rhs.end(); ) {
		value_type ele = ::fabs(*itr++);
		if(ele > max_ele) max_ele = ele;
	}
	return max_ele;
}

// ////////////////////////////////////////////////////////////////////////////////////////
// Misc operations

// swap(vs1, vs2)
void LinAlgPack::swap( VectorSlice* vs1, VectorSlice* vs2 ) {
	if( vs1->overlap(*vs2) == SAME_MEM ) return;
	VopV_assert_sizes( vs1->size(), vs2->size() );
	BLAS_Cpp::swap( vs1->size(), vs1->raw_ptr(), vs1->stride(), vs2->raw_ptr(), vs2->stride() );
}
