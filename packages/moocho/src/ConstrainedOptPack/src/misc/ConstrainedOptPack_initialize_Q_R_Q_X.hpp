// ///////////////////////////////////////////////////////////////
// initialize_Q_R_Q_X.h

#ifndef INITIALIZE_Q_R_Q_X_H
#define INITIALIZE_Q_R_Q_X_H

#include "ConstrainedOptimizationPackTypes.h"

namespace ConstrainedOptimizationPack {

///
/** Initialize #GenPermMatrixSlice# mapping matrices for #Q_R# and #Q_X#.
 *
 *
 * @param  n_R        [in] Number of free variables
 * @param  n_X        [in] Number of fixed variables
 * @param  i_x_free   [in] array (length n_R) of indices of free variables.
 *                    If n_R == 0 then i_x_free can be NULL.  If n_X == 0
 *                    then i_x_free can be NULL in which case it is considered
 *                    identitiy.
 * @param  i_x_fixed  [in] array (length n_X) of indices of fixed variables.
 *                    If n_X == 0 then i_x_fixed can be NULL.  If n_R == 0
 *                    then i_x_fixed can be NULL in which case it is considered
 *                    identitiy.
 * @param  test_setup [in] If true then i_x_free[] and i_x_fixed[] will be
 *                    validated and if not okay then an exception will be
 *                    thown.
 * @param  Q_R_row_i  [out] array (length n_R) of row indices for Q_R.
 *                    If n_R == 0 then Q_R_row_i can be NULL.
 * @param  Q_R_col_j  [out] array (length n_R) of column indices for Q_R.
 *                    If n_R == 0 then Q_R_col_j can be NULL.
 * @param  Q_R        [out] GenPermMatixSlice object initialized with 
 *                    Q_R_row_i and Q_R_col_j.  If n_R == 0 then Q_R
 *                    will be initialized to (n_X by 0).
 * @param  Q_X_row_i  [out] array (length n_X) of row indices for Q_X.
 *                    If n_X == 0 then Q_X_row_i can be NULL.
 * @param  Q_X_col_j  [out] array (length n_X) of column indices for Q_X
 *                    If n_X == 0 then Q_X_col_j can be NULL.
 * @param  Q_X        [out] GenPermMatixSlice object initialized with 
 *                    Q_X_row_i and Q_X_col_j  If n_X == 0 then Q_X
 *                    will be initialized to (n_X by 0).
 */
void initialize_Q_R_Q_X(
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
	);

}  // end namespace ConstrainedOptimizationPack

#endif INITIALIZE_Q_R_Q_X_H
