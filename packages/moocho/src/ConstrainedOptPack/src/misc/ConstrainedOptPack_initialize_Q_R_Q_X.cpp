// ///////////////////////////////////////////////////////////////
// initialize_Q_R_Q_X.cpp

#include "ConstrainedOptimizationPack/include/initialize_Q_R_Q_X.h"
#include "SparseLinAlgPack/include/GenPermMatrixSlice.h"

void ConstrainedOptimizationPack::initialize_Q_R_Q_X(
	size_type            n_R
	,size_type           n_X
	,const size_type     i_x_free[]
	,const size_type     i_x_fixed[]
	,bool                test_setup
	,size_type           Q_R_row_i[]
	,size_type           Q_R_col_j[]
	,GenPermMatrixSlice  *Q_R
	,size_type           Q_X_row_i[]
	,size_type           Q_X_col_j[]
	,GenPermMatrixSlice  *Q_X
	)
{
	namespace GPMSIP = SparseLinAlgPack::GenPermMatrixSliceIteratorPack;
	const size_type
		n = n_R + n_X;
	// Validate i_x_free[] and i_x_fixed[] are correct!
	// ToDo: Implement this!

	// Setup Q_R
	if( n_X == 0 && i_x_free == NULL ) {
		// consider i_x_free as identity
		size_type
			*row_i_itr = Q_R_row_i,
			*col_j_itr = Q_R_col_j;
		for( size_type i = 1; i <= n_R; ++i, ++row_i_itr, ++col_j_itr ) {
			*row_i_itr = i;
			*col_j_itr = i;
		}
	}
	else if( n_R > 0 ) {
		const size_type
			*i_x_R = i_x_free;
		size_type
			*row_i_itr = Q_R_row_i,
			*col_j_itr = Q_R_col_j;
		for( size_type i = 1; i <= n_R; ++i, ++i_x_R, ++row_i_itr, ++col_j_itr ) {
			*row_i_itr = *i_x_R;
			*col_j_itr = i;
		}
	}
	Q_R->initialize_and_sort(
		n,n_R,n_R,0,0,GPMSIP::BY_ROW
		,Q_R_row_i,Q_R_col_j,test_setup);
	// Setup Q_X
	if( n_R == 0 && i_x_fixed == NULL ) {
		// consider i_x_fixed as identity
		size_type
			*row_i_itr = Q_X_row_i,
			*col_j_itr = Q_X_col_j;
		for( size_type i = 1; i <= n_X; ++i, ++row_i_itr, ++col_j_itr ) {
			*row_i_itr = i;
			*col_j_itr = i;
		}
	}
	else if( n_X > 0 ) {
		const size_type
			*i_x_X = i_x_fixed;
		size_type
			*row_i_itr = Q_X_row_i,
			*col_j_itr = Q_X_col_j;
		for( size_type i = 1; i <= n_X; ++i, ++i_x_X, ++row_i_itr, ++col_j_itr ) {
			*row_i_itr = *i_x_X;
			*col_j_itr = i;
		}					
	}
	Q_X->initialize_and_sort(
		n,n_X,n_X,0,0,GPMSIP::BY_ROW
		,Q_X_row_i,Q_X_col_j,test_setup);
}
