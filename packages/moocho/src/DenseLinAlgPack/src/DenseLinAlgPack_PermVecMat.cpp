// ///////////////////////////////////////////////////////////////////////////////////////
// PermVecMat.cpp

#include "../include/PermVecMat.h"
#include "../include/GenMatrixClass.h"
#include "../include/IVector.h"

// Local assert function
namespace {
inline void i_assert_perm_size(size_t size1, size_t size2)
{
#ifdef LINALGPACK_CHECK_RHS_SIZES
	if(size1 != size2)
		throw std::length_error("The size of the permutation vector is not correct");
#endif
}

}

void LinAlgPack::identity_perm(IVector* perm) {
	if(!perm->size())
		throw std::length_error("LinAlgPack::identity_perm(): perm must be sized");
	IVector::iterator	itr_perm		= perm->begin();
	for(size_type i = 1; i <= perm->size(); ++i)
		*itr_perm++ = i;
}

void LinAlgPack::inv_perm(const IVector& perm, IVector* inv_perm) {
	inv_perm->resize(perm.size());
	for(size_type i = 1; i <= perm.size(); ++i)
		(*inv_perm)(perm(i)) = i;
}

void LinAlgPack::perm_ele(const IVector& perm, VectorSlice* vs)
{
	i_assert_perm_size(vs->size(),perm.size());
	Vector tmp_v(vs->size());
	Vector::iterator		v_itr		= tmp_v.begin(),
							v_itr_end	= tmp_v.end();
	IVector::const_iterator perm_itr = perm.begin();
	// Copy the elements in the permed order into the temp vector
	for(; v_itr != v_itr_end; ++v_itr, ++perm_itr)
		*v_itr = (*vs)(*perm_itr);
		
	// Copy the permed vector back
	(*vs) = tmp_v;
}

void LinAlgPack::perm_ele(const VectorSlice& x, const IVector& perm, VectorSlice* y)
{
	i_assert_perm_size(x.size(),perm.size());
	i_assert_perm_size(y->size(),perm.size());

	IVector::const_iterator
		perm_itr	= perm.begin();
	VectorSlice::iterator
		y_itr		= y->begin(),
		y_end		= y->end();
	while(y_itr != y_end)
		*y_itr++ = x(*perm_itr++);
}

void LinAlgPack::inv_perm_ele(const VectorSlice& y, const IVector& perm, VectorSlice* x)
{
	i_assert_perm_size(y.size(),perm.size());
	i_assert_perm_size(x->size(),perm.size());

	VectorSlice::const_iterator
		y_itr		= y.begin(),
		y_end		= y.end();
	IVector::const_iterator
		perm_itr	= perm.begin();
	while(y_itr != y_end)
		(*x)(*perm_itr++) = *y_itr++;
}

void LinAlgPack::perm_rows(const IVector& row_perm, GenMatrixSlice* gms)
{
	i_assert_perm_size(gms->rows(),row_perm.size());
	GenMatrix tmp_gm(gms->rows(),gms->cols());
	GenMatrixSlice::size_type rows = gms->rows(), i;
	// Copy the rows in the correct order into the temp matrix.
	for(i = 1; i <= rows; ++i)
		tmp_gm.row(i) = gms->row(row_perm(i));
	// Copy the permed matrix back
	(*gms) = tmp_gm;
}

void LinAlgPack::perm_cols(const IVector& col_perm, GenMatrixSlice* gms)
{
	i_assert_perm_size(gms->cols(),col_perm.size());
	GenMatrix tmp_gm(gms->rows(),gms->cols());
	GenMatrixSlice::size_type cols = gms->cols(), i;
	// Copy the columns in the correct order into the temp matrix.
	for(i = 1; i <= cols; ++i)
		tmp_gm.col(i) = gms->col(col_perm(i));
	// Copy the permed matrix back
	(*gms) = tmp_gm;
}

void LinAlgPack::perm_rows_cols(const IVector& row_perm, const IVector& col_perm
	, GenMatrixSlice* gms)
{
	i_assert_perm_size(gms->rows(),row_perm.size());
	i_assert_perm_size(gms->cols(),col_perm.size());
	GenMatrix tmp_gm(gms->rows(),gms->cols());
	GenMatrixSlice::size_type rows = gms->rows(), cols = gms->cols(), i;
	// Copy the rows in the correct order into the temp matrix.
	for(i = 1; i <= rows; ++i)
		tmp_gm.row(i) = gms->row(row_perm(i));
	// Copy the columns in the correct order back into matrix.
	for(i = 1; i <= cols; ++i)
		gms->col(i) = tmp_gm.col(col_perm(i));
}
