// /////////////////////////////////////////////////////////////////////////////
// DecompositionSystemCoordinate.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	
// disable VC 5.0 warnings about truncated identifier names (templates).
#pragma warning(disable : 4503)	

#include "../include/DecompositionSystemCoordinate.h"
#include "../include/IdentZeroVertConcatMatrixSubclass.h"
#include "SparseSolverPack/test/TestBasisSystem.h"
#include "SparseLinAlgPack/include/GenMatrixSubclass.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "Misc/include/update_success.h"

namespace ConstrainedOptimizationPack {

void DecompositionSystemCoordinate::update_decomp(MatrixWithOp* A, MatrixWithOp* Z, MatrixWithOp* Y
		, MatrixWithOp* U, MatrixWithOp* V)
{

	// Validate the matrix types and get concrete references.

	// Z
	validate_Z(Z);

	// Y
	IdentZeroVertConcatMatrixSubclass*
		_Y = dynamic_cast<IdentZeroVertConcatMatrixSubclass*>(Y);
	if(!_Y)
		throw InvalidMatrixType( "DecompositionSystemCoordinate::update_decomp(...):  The concrete type "
									" of the Y matrix must be a subclass of IdentZeroVertConcatMatrixSubclass" );
	IdentZeroVertConcatMatrix &cY = _Y->m();

	// V
	GenMatrixSubclass*
		_V = dynamic_cast<GenMatrixSubclass*>(V);
	if(!_V)
		throw InvalidMatrixType( "DecompositionSystemCoordinate::update_decomp(...):  The concrete type "
									" of the V matrix must be a subclass of GenMatrixSubclass" );
	GenMatrix &cV = _V->m();
	
	bool new_basis = factor_a_new_basis();

	// Factor the basis (the type of A is validated by basys_sys().
	factor(A);

	if(new_basis) {
		// Setup the access to C, N, E and F.
		// Here the type of the E matrix is determined by the basys_sys() and done so externally to this
		// class.

		if( !access_matrices(BasisSystem::C) ) {
			// If part_[BasisSystem::C] == 0 then we need to allocate all new matrices.
			
			bool allocate[BasisSystem::NUM_ACCESS_MATRICES];
			
			std::fill_n( allocate, (int)BasisSystem::NUM_ACCESS_MATRICES, true );
			
			allocate[BasisSystem::E] = false;	// E is U so don't allocate it
			
			basis_sys().create_access_matrices( allocate, &access_matrices(BasisSystem::C) );

			// Use U for E to be set up by basis_sys()
			access_matrices(BasisSystem::E) = U;
		}

		// Setup access to C, N and perhaps E and F if they exist
		// If V and therfore E is not of the proper type then this will
		// throw an exception.
		try {
			basis_sys().setup_access_matrices( &access_matrices(BasisSystem::C) );
		}
		catch(const BasisSystem::InvalidMatrixType& excpt) {
			throw InvalidMatrixType( excpt.what() );
		}
	}

	// Test the basis system
	if( check_results() && get_basis_sys_tester().get() ) {
		if(!basis_sys_tester().check_basis_system(
			basis_sys()
			, const_cast<const MatrixWithOp**>(&access_matrices(BasisSystem::C))
			, *A, trase(), out() ))
		{
			throw std::runtime_error( "DecompositionSystemCoordinate::update_decomp(...) : "
				"Error, BasisSystem object does not check out." );
		}
	}

	// Update the matrices Z, Y, U and V

	// * Z, V
	update_Z_and_V(Z,&cV);

	// * Y = [I, 0]'

	cY.resize( n(), r(), true );

	// * U is the same thing as E and was already set up by basis_sys()

}

void DecompositionSystemCoordinate::solve_transAtY(const VectorSlice& b, BLAS_Cpp::Transp trans
	, VectorSlice* x) const
{
	solve_C(b,trans,x);
}

void DecompositionSystemCoordinate::solve_transAtY(const GenMatrixSlice& B, BLAS_Cpp::Transp trans
	, GenMatrixSlice* X) const
{
	solve_C(B,trans,X);
}

void DecompositionSystemCoordinate::delete_access_matrices() {
	access_matrices(BasisSystem::E) = 0;	// E which is U is not ours to destroy
	basis_sys().destroy_access_matrices( &access_matrices(BasisSystem::C) );
}
		
}	// end namespace ConstrainedOptimizationPack
