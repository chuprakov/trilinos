// /////////////////////////////////////////////////////////////////
// InvCholUpdate.h

#ifndef INV_CHOL_UPDATE_H
#define INV_CHOL_UPDATE_H

#include "LinAlgPackTypes.h"

namespace LinAlgPack {

///
/** Perform an update of the Cholesky factor of a symetric matrix (SM = L * L^T)
  * and while maintaining the triangular structure.
  *
  * This function permforms an update of L^T or L^-1 depending on whether
  * we are updating the cholesky factor or its inverse.  This function
  * implements algorithm A3.4.1 in Dennis and Schabel.  The update is:
  * (J_new^T = L_old^T + u * v^T) where J_new is rotated back to triangular form.
  *
  * Preconditions: \begin{itemize}
  * \item UpTriM.rows() == UpTriM.cols() == u.size() == v.size() (throw std::length_error)
  * \end{itemize}
  */
void update_chol_factor(GenMatrixSlice* UpTriM, VectorSlice* u
	, const VectorSlice& v);

///
/** Perform a jacobi rotation or a matrix about row i.
  *
  * Preconditions: \begin{itemize}
  * \item UpTriM.rows() == UpTriM.cols() (throw std::length_error)
  * \end{itemize}
  */
void jacobi_rotate(GenMatrixSlice* UpTriM, size_type row_i, value_type alpha
	, value_type beta); 

}	// end namespace LinAlgPack

#endif	// INV_CHOL_UPDATE_H
