// //////////////////////////////////////////////////////////////
// MatrixWithOp.cpp

#include <assert.h>

#include <typeinfo>
#include <stdexcept>

#include "AbstractLinAlgPack/include/MatrixWithOp.h"
#include "AbstractLinAlgPack/include/MatrixWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/SpVectorClass.h"
#include "AbstractLinAlgPack/include/SpVectorView.h"
#include "AbstractLinAlgPack/include/EtaVector.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "ThrowException.h"

namespace AbstractLinAlgPack {

MatrixWithOp::mat_ptr_t
MatrixWithOp::sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
	if( 
		( ( row_rng.lbound() == 1 && row_rng.ubound() == this->rows() )
		  || row_rng.full_range() )
		&&
		( ( col_rng.lbound() == 1 && col_rng.ubound() == this->cols() )
		  || row_rng.full_range() )
		)
		return mat_ptr_t(this,false); // don't clean up memory
	return NULL; // requested a view that was not the entire matrix!
}

MatrixWithOp& MatrixWithOp::zero_out()
{
	THROW_EXCEPTION(
		true, std::logic_error, "MatrixWithOp::zero_out(): "
		"Error, this method as not been defined by the subclass \'"
		<<typeid(*this).name()<<"\'" );
}

MatrixWithOp& MatrixWithOp::operator=(const MatrixWithOp& M)
{
	assert(0); // ToDo: Implement!
	return *this;
}

std::ostream& MatrixWithOp::output(std::ostream& out) const
{
	const size_type m = this->rows(), n = this->cols();
	VectorSpace::vec_mut_ptr_t
		row_vec = space_rows().create_member(); // dim() == n
	out << m << " " << n << std::endl;
	for( size_type i = 1; i <= m; ++i ) {
		LinAlgOpPack::V_StMtV( row_vec.get(), 1.0, *this, BLAS_Cpp::trans, EtaVector(i,m)() );
		row_vec->output(out,false,true);
	}
	return out;
}

// Level-1 BLAS

// rhs matrix argument

void MatrixWithOp::Mp_StM(
	MatrixWithOp* m_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs) const
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;

#ifdef _DEBUG
	THROW_EXCEPTION(
		!m_lhs->space_rows().is_compatible(
			trans_rhs == no_trans ? this->space_rows() : this->space_cols() )
		|| !m_lhs->space_cols().is_compatible(
			trans_rhs == no_trans ? this->space_cols() : this->space_rows() )
		, IncompatibleMatrices
		,"MatrixWithOp::Mp_StM(m_lhs,...): Error, m_lhs of type \'"<<typeid(*m_lhs).name()<<"\' "
		<<"is not compatible with this of type \'"<<typeid(*m_lhs).name()<<"\'" );
#endif
	MatrixWithOpMutable
		*m_mut_lhs = dynamic_cast<MatrixWithOpMutable*>(m_lhs);
	if(!m_mut_lhs)
		return m_lhs->Mp_StM(alpha,*this,trans_rhs);
		
	const size_type
		rows = BLAS_Cpp::rows( m_mut_lhs->rows(), m_mut_lhs->cols(), trans_rhs ),
		cols = BLAS_Cpp::cols( m_mut_lhs->rows(), m_mut_lhs->cols(), trans_rhs );
	for( size_type j = 1; j <= cols; ++j )
		AbstractLinAlgPack::Vp_StMtV( m_mut_lhs->col(j).get(), alpha, *this, trans_rhs, EtaVector(j,cols)() );
}

void MatrixWithOp::Mp_StMtP(
	MatrixWithOp* m_lhs, value_type alpha
	, BLAS_Cpp::Transp M_trans
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	) const
{
	assert(0); // ToDo: Implement!
}

void MatrixWithOp::Mp_StPtM(
	MatrixWithOp* m_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	, BLAS_Cpp::Transp M_trans
	) const
{
	assert(0); // ToDo: Implement!
}

void MatrixWithOp::Mp_StPtMtP(
	MatrixWithOp* m_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_trans
	, const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
	) const
{
	assert(0); // ToDo: Implement!
}

// lhs matrix argument

void MatrixWithOp::Mp_StM(
	value_type alpha,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs)
{
	assert(0); // ToDo: Implement!
}

void MatrixWithOp::Mp_StMtP(
	value_type alpha
	,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
	,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	)
{
	assert(0); // ToDo: Implement!
}

void MatrixWithOp::Mp_StPtM(
	value_type alpha
	,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
	)
{
	assert(0); // ToDo: Implement!
}

void MatrixWithOp::Mp_StPtMtP(
	value_type alpha
	,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
	,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
	)
{
	assert(0); // ToDo: Implement!
}

// Level-2 BLAS

void MatrixWithOp::Vp_StMtV(
	VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2, value_type beta) const
{
	Vp_MtV_assert_sizes(v_lhs->dim(), this->rows(), this->cols(), trans_rhs1, sv_rhs2.dim() );
	VectorSpace::vec_mut_ptr_t
		v_rhs2 = (trans_rhs1 == BLAS_Cpp::no_trans
				  ? this->space_rows()
				  : this->space_cols()
			).create_member();
	v_rhs2->set_sub_vector(sub_vec_view(sv_rhs2));
	this->Vp_StMtV(v_lhs,alpha,trans_rhs1,*v_rhs2,beta);
}

void MatrixWithOp::Vp_StPtMtV(
	VectorWithOpMutable* v_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_rhs2_trans
	, const VectorWithOp& vs_rhs3, value_type beta) const
{
	assert(0); // ToDo: Implement!
}

void MatrixWithOp::Vp_StPtMtV(
	VectorWithOpMutable* v_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_rhs2_trans
	, const SpVectorSlice& sv_rhs3, value_type beta) const
{
	assert(0); // ToDo: Implement!
}

value_type MatrixWithOp::transVtMtV(
	const VectorWithOp& vs_rhs1, BLAS_Cpp::Transp trans_rhs2
	, const VectorWithOp& vs_rhs3) const
{
	assert(0); // ToDo: Implement!
	return 0.0;
}

value_type MatrixWithOp::transVtMtV(
	const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
	, const SpVectorSlice& sv_rhs3) const
{
	assert(0); // ToDo: Implement!
	return 0.0;
}

void MatrixWithOp::syr2k(
	BLAS_Cpp::Transp M_trans, value_type alpha
	, const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
	, const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
	, value_type beta, MatrixSymWithOp* sym_lhs ) const
{
	assert(0); // ToDo: Implement!
}

// Level-3 BLAS

void MatrixWithOp::Mp_StMtM(
	MatrixWithOp* m_lhs, value_type alpha
	, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOp& mwo_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
	assert(0); // ToDo: Implement!
}

void MatrixWithOp::Mp_StMtM(
	MatrixWithOp* m_lhs, value_type alpha
	, const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const
{
	assert(0); // ToDo: Implement!
}

void MatrixWithOp::Mp_StMtM(
	value_type alpha
	,const MatrixWithOp& mvw_rhs1, BLAS_Cpp::Transp trans_rhs1
	,const MatrixWithOp& mwo_rhs2,BLAS_Cpp::Transp trans_rhs2
	,value_type beta )
{
	assert(0); // ToDo: Implement!
}

void MatrixWithOp::syrk(
	BLAS_Cpp::Transp M_trans, value_type alpha
	, value_type beta, MatrixSymWithOp* sym_lhs ) const
{
	assert(0); // ToDo: Implement!
}

// overridden from MatrixBase

size_type MatrixWithOp::rows() const
{
	return this->space_cols().dim();
}

size_type MatrixWithOp::cols() const
{
	return this->space_rows().dim();
}

} // end namespace AbstractLinAlgPack
