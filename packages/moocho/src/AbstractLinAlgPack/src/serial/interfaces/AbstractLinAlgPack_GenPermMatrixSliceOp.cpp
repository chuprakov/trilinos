// /////////////////////////////////////////////////////////
// GenPermMatrixSliceOp.cpp

#include <assert.h>

#include "SparseLinAlgPack/include/GenPermMatrixSliceOp.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/LinAlgPackAssertOp.h"

void SparseLinAlgPack::V_StMtV(
	  SpVector* y, value_type a, const GenPermMatrixSlice& P
	, BLAS_Cpp::Transp P_trans, const VectorSlice& x )
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	namespace GPMSIP = GenPermMatrixSliceIteratorPack;
	using LinAlgPack::MtV_assert_sizes;
	MtV_assert_sizes( P.rows(), P.cols(), P_trans, x.size() );

	y->resize( BLAS_Cpp::rows( P.rows(), P.cols(), P_trans ), P.nz() );

	typedef SpVector::element_type ele_t;

	if( P_trans == no_trans ) {
		for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
			const size_type
				i = itr->row_i(),
				j = itr->col_j();
			y->add_element( ele_t( i, a * x(j) ) );
		}
	}
	else {
		for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
			const size_type
				j = itr->row_i(),
				i = itr->col_j();
			y->add_element( ele_t( i, a * x(j) ) );
		}
	}
	if(		( P_trans == no_trans	&& P.ordered_by() == GPMSIP::BY_ROW )
		||	( P_trans == trans		&& P.ordered_by() == GPMSIP::BY_COL )	)
	{
		y->assume_sorted(true);
	}
}

void SparseLinAlgPack::Vp_StMtV(
	  SpVector* y, value_type a, const GenPermMatrixSlice& P
	, BLAS_Cpp::Transp P_trans, const VectorSlice& x )
{
	using LinAlgPack::Vp_MtV_assert_sizes;
	Vp_MtV_assert_sizes( y->size(), P.rows(), P.cols(), P_trans, x.size() );

	typedef SpVector::element_type ele_t;

	if( P_trans == BLAS_Cpp::no_trans ) {
		for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
			const size_type
				i = itr->row_i(),
				j = itr->col_j();
			y->add_element( ele_t( i, a * x(j) ) );
		}
	}
	else {
		for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
			const size_type
				j = itr->row_i(),
				i = itr->col_j();
			y->add_element( ele_t( i, a * x(j) ) );
		}
	}
}

void SparseLinAlgPack::Vp_StMtV(
	  VectorSlice* y, value_type a, const GenPermMatrixSlice& P
	, BLAS_Cpp::Transp P_trans, const VectorSlice& x, value_type b = 1.0)
{
	using LinAlgPack::Vp_MtV_assert_sizes;
	Vp_MtV_assert_sizes( y->size(), P.rows(), P.cols(), P_trans, x.size() );
	// y = b*y
	if( b == 0.0 )
		*y = 0.0;
	else
		Vt_S(y,b);	
	// y += a*op(P)*x
	if( P_trans == BLAS_Cpp::no_trans ) {
		for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
			const size_type
				i = itr->row_i(),
				j = itr->col_j();
			(*y)(i) += a * x(j);
		}
	}
	else {
		for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
			const size_type
				j = itr->row_i(),
				i = itr->col_j();
			(*y)(i) += a * x(j);
		}
	}
}

void SparseLinAlgPack::Vp_StMtV(
	  VectorSlice* y, value_type a, const GenPermMatrixSlice& P
	, BLAS_Cpp::Transp P_trans, const SpVectorSlice& x, value_type b = 1.0)
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	namespace GPMSIP = GenPermMatrixSliceIteratorPack;
	using LinAlgPack::Vp_MtV_assert_sizes;
	
	Vp_MtV_assert_sizes( y->size(), P.rows(), P.cols(), P_trans, x.size() );
	// y = b*y
	if( b == 0.0 )
		*y = 0.0;
	else
		Vt_S(y,b);	
	// y += a*op(P)*x
	if( x.is_sorted() ) {
		const SpVectorSlice::difference_type x_off = x.offset();
		if( P_trans == no_trans && P.ordered_by() == GPMSIP::BY_COL ) {
			assert(0);	// ToDo: implement this!
		}
		else if( P_trans == trans && P.ordered_by() == GPMSIP::BY_ROW ) {
			GenPermMatrixSlice::const_iterator
				P_itr = P.begin(),
				P_end = P.end();
			SpVectorSlice::const_iterator
				x_itr = x.begin(),
				x_end = x.end();
			while( P_itr != P_end && x_itr != x_end ) {
				const size_type
					j = P_itr->row_i(),
					i = P_itr->col_j();
				if( j < x_itr->indice() + x_off ) {
					++P_itr;
					continue;
				}
				else if( j > x_itr->indice() + x_off ) {
					++x_itr;
					continue;
				}
				else {	// they are equal
					(*y)(i) += a * x_itr->value();
					++P_itr;
					++x_itr;
				}
			}
		}
	}
	else {
		// Since things do not match up we will have to create a
		// temporary dense copy of x to operate on.
		assert(0);	// ToDo: Implement this!
	}
}
