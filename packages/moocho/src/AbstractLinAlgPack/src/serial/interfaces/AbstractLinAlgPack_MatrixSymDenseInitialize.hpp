// //////////////////////////////////////////////////////////
// MatrixSymDenseInitialize.h

#ifndef MATRIX_SYM_DENSE_INITIALIZE_H
#define MATRIX_SYM_DENSE_INITIALIZE_H

#include "Matrix.h"

namespace SparseLinAlgPack {

///
/** Interface for initializing a matrix with a dense symmetric matrix.
  */
class MatrixSymDenseInitialize : public virtual Matrix {
public:

	///
	/** Initialize with a symmetric dense matrix.
	 * Through this interface there are absolutly no preconditions
	 * as the the state of #this# after this function executes.
	 * The implementation can use #M# to initialize itself any way
	 * it would like.
	 */
	virtual void initialize( const sym_gms& M ) = 0;

};	// end class MatrixSymDenseInitialize

}	// end namespace SparseLinAlgPack 

#endif	// MATRIX_SYM_DENSE_INITIALIZE_H
