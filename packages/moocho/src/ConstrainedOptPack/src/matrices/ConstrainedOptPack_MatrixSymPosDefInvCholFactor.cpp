// //////////////////////////////////////////////////////////////////////////////////
// MatrixSymPosDefInvCholFactor.cpp

#include "../include/MatrixSymPosDefInvCholFactor.h"
#include "../include/SymInvCholMatrixOp.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/GenMatrixOp.h"
#include "LinAlgPack/include/GenMatrixOut.h"

namespace LinAlgOpPack {

using SparseLinAlgPack::Vp_StV;
using SparseLinAlgPack::Vp_StMtV;
using SparseLinAlgPack::Mp_StM;
using ConstrainedOptimizationPack::Vp_StMtV;

}	// end namespace LinAlgOpPack

namespace ConstrainedOptimizationPack {

// Overridden from Matrix 

size_type MatrixSymPosDefInvCholFactor::cols() const
{
	return rows();
}

// Overridden from MatrixWithOp

MatrixWithOp& MatrixSymPosDefInvCholFactor::operator=(const MatrixWithOp& m)
{
	return MatrixWithOpConcreteEncap<SymInvCholMatrix>::operator=(m);
}

std::ostream& MatrixSymPosDefInvCholFactor::output(std::ostream& out) const
{	return out << "\n*** Inverse Cholesky factor:\n" << m().UInv(); }

// Level-2 BLAS

void MatrixSymPosDefInvCholFactor::Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const VectorSlice& vs_rhs2, value_type beta) const
{
	ConstrainedOptimizationPack::Vp_StMtV(vs_lhs,alpha,m(),trans_rhs1,vs_rhs2,beta);
}

void MatrixSymPosDefInvCholFactor::Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2, value_type beta) const
{
	using LinAlgOpPack::assign;
	Vector vs_rhs2;
	assign(&vs_rhs2,sv_rhs2);
	ConstrainedOptimizationPack::Vp_StMtV(vs_lhs,alpha,m(),trans_rhs1,vs_rhs2,beta);
}

value_type MatrixSymPosDefInvCholFactor::transVtMtV(const VectorSlice& vs_rhs1, BLAS_Cpp::Transp trans_rhs2
	, const VectorSlice& vs_rhs3) const
{
	return ConstrainedOptimizationPack::transVtMtV(vs_rhs1,m(),vs_rhs3);
}

value_type MatrixSymPosDefInvCholFactor::transVtMtV(const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
	, const SpVectorSlice& sv_rhs3) const
{
	using LinAlgOpPack::assign;
	Vector vs_rhs1, vs_rhs3;
	assign(&vs_rhs1,sv_rhs1);
	assign(&vs_rhs3,sv_rhs3);
	return ConstrainedOptimizationPack::transVtMtV(vs_rhs1,m(),vs_rhs3);
}

// Overridden from MatrixFactorized

void MatrixSymPosDefInvCholFactor::V_InvMtV(Vector* v_lhs, BLAS_Cpp::Transp trans_rhs1
	, const VectorSlice& vs_rhs2) const
{
	ConstrainedOptimizationPack::V_InvMtV(v_lhs,m(),vs_rhs2);
}

void MatrixSymPosDefInvCholFactor::V_InvMtV(VectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
	, const VectorSlice& vs_rhs2) const
{
	ConstrainedOptimizationPack::V_InvMtV(vs_lhs,m(),vs_rhs2);
}

void MatrixSymPosDefInvCholFactor::V_InvMtV(Vector* v_lhs, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2) const
{
	ConstrainedOptimizationPack::V_InvMtV(v_lhs,m(),sv_rhs2);
}

void MatrixSymPosDefInvCholFactor::V_InvMtV(VectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2) const
{
	ConstrainedOptimizationPack::V_InvMtV(vs_lhs,m(),sv_rhs2);
}

value_type MatrixSymPosDefInvCholFactor::transVtInvMtV(const VectorSlice& vs_rhs1
	, BLAS_Cpp::Transp trans_rhs2, const VectorSlice& vs_rhs3) const
{
	return ConstrainedOptimizationPack::transVtInvMtV(vs_rhs1,m(),vs_rhs3);
}

value_type MatrixSymPosDefInvCholFactor::transVtInvMtV(const SpVectorSlice& sv_rhs1
	, BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3) const
{
	return ConstrainedOptimizationPack::transVtInvMtV(sv_rhs1,m(),sv_rhs3);}

// Overridden from MatrixSymFactorized

void MatrixSymPosDefInvCholFactor::M_StMtInvMtM(
	  sym_gms* S, value_type a, const MatrixWithOp& B
	, BLAS_Cpp::Transp B_trans, EMatrixDummyArg dummy_arg ) const
{
//	// Uncomment to use the defalut implementation (for debugging)
//	MatrixSymFactorized::M_StMtInvMtM(S,a,B,B_trans,dummy_arg); return;

	using BLAS_Cpp::trans;
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans_not;
	using BLAS_Cpp::upper;
	using BLAS_Cpp::nonunit;
	using SparseLinAlgPack::M_StInvMtM;
	using LinAlgPack::tri;
	using LinAlgPack::syrk;
	using LinAlgPack::M_StInvMtM;
	using LinAlgOpPack::M_StMtM;
	using LinAlgOpPack::assign;

	LinAlgPack::MtM_assert_sizes( rows(), cols(), no_trans
		, B.rows(), B.cols(), trans_not(B_trans) );
	LinAlgPack::Mp_MtM_assert_sizes( S->rows(), S->cols(), no_trans
		, B.rows(), B.cols(), B_trans
		, B.rows(), B.cols(), trans_not(B_trans) );
	//
	// S = a * op(B) * inv(M) * op(B)'
	// 
	// M = L * L'
	// inv(M) = inv(L * L') = inv(L') * inv(L) = UInv * UInv'
	// 
	// S = a * op(B) * UInv * UInv' * op(B)'
	// 
	// T = op(B)'
	// 
	// T = UInv' * T (inplace with BLAS)
	// 
	// S = a * T' * T
	// 

	// T = op(B)'
	GenMatrix T;
	assign( &T, B, trans_not(B_trans) );
	// T = UInv' * T (inplace with BLAS)
	M_StMtM( &T(), 1.0, tri(m().UInv(),upper,nonunit), trans, T(), no_trans );
	// S = a * T' * T
	syrk( trans, a, T(), 0.0, S );
}

// Overridden from MatrixSymSecantUpdateable

void MatrixSymPosDefInvCholFactor::init_identity(size_type n, value_type alpha)
{
	if( alpha <= 0.0 ) {
		std::ostringstream omsg;
		omsg	<< "MatrixSymPosDefInvCholFactor::init_identity(...) : Error, alpha = " << alpha
				<< " <= 0.0 and therefore this is not a positive definite matrix.";
		throw UpdateSkippedException( omsg.str() );	
	}
	m().resize(n);
	m().UInv() = 0.0;
	m().UInv().diag() = 1.0 / ::sqrt( alpha );
}

void MatrixSymPosDefInvCholFactor::init_diagonal( const VectorSlice& diag )
{
	VectorSlice::const_iterator
		min_ele_ptr = std::min_element( diag.begin(), diag.end() );
	if( *min_ele_ptr <= 0.0 ) {
		std::ostringstream omsg;
		omsg	<< "MatrixSymPosDefInvCholFactor::init_diagonal(...) : Error, "
				<< "diag(" << min_ele_ptr - diag.begin() + 1 << " ) = "
				<< (*min_ele_ptr) << " <= 0.0.\n"
				<< "Therefore this is not a positive definite matrix.";
		throw UpdateSkippedException( omsg.str() );	
	}
	const size_type n = diag.size();
	m().resize(n);
	m().UInv() = 0.0;

	VectorSlice::const_iterator
		diag_itr = diag.begin();
	VectorSlice::iterator
		inv_fact_diag_itr = m().UInv().diag().begin();

	while( diag_itr != diag.end() )
		*inv_fact_diag_itr++ = 1.0 / ::sqrt( *diag_itr++ );
}

void MatrixSymPosDefInvCholFactor::secant_update(VectorSlice* s, VectorSlice* y, VectorSlice* _Bs)
{
	using LinAlgOpPack::V_MtV;
	try {
		if(!_Bs) {
			Vector Bs;
			V_MtV( &Bs, *this, BLAS_Cpp::no_trans, *s );
			ConstrainedOptimizationPack::BFGS_update(&m(),s,y,&Bs());
		}
		else {
			ConstrainedOptimizationPack::BFGS_update(&m(),s,y,_Bs);
		}
	}
	catch(const BFGSUpdateSkippedException& excpt) {
		throw UpdateSkippedException( excpt.what() );
	}
}

// Overridden from MatrixExtractInvCholFactor

void MatrixSymPosDefInvCholFactor::extract_inv_chol( tri_ele_gms* InvChol ) const
{
	LinAlgPack::assign( InvChol, LinAlgPack::tri_ele( m().UInv(), BLAS_Cpp::upper ) );
}

}	// end namespace ConstrainedOptimizationPack
