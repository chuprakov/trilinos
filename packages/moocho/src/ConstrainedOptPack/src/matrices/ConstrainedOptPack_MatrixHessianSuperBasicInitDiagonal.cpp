// ////////////////////////////////////////////////////////
// MatrixHessianSuperBasicInitDiagonal.cpp

#include "ConstrainedOptimizationPack/include/MatrixHessianSuperBasicInitDiagonal.h"
#include "LinAlgPack/include/VectorClass.h"
#include "Misc/include/dynamic_cast_verbose.h"

namespace ConstrainedOptimizationPack {

MatrixHessianSuperBasicInitDiagonal::MatrixHessianSuperBasicInitDiagonal()
	: B_RR_init_(NULL)
{}

void MatrixHessianSuperBasicInitDiagonal::initialize(
	size_type            n
	,size_type           n_R
	,const size_type     i_x_free[]
	,const size_type     i_x_fixed[]
	,const EBounds       bnd_fixed[]
	,const B_RR_ptr_t&   B_RR_ptr
	,const B_RX_ptr_t&   B_RX_ptr
	,BLAS_Cpp::Transp    B_RX_trans
	,const B_XX_ptr_t&   B_XX_ptr
	)
{
	using DynamicCastHelperPack::dyn_cast;

	// Validate the B_RR supports this interface
#ifdef _WINDOWS
	B_RR_init_ = &dynamic_cast<MatrixSymInitDiagonal&>(
		const_cast<MatrixSymWithOpFactorized&>(*B_RR_ptr)
		);
#else
	B_RR_init_ = &dyn_cast<MatrixSymInitDiagonal>(
		const_cast<MatrixSymWithOpFactorized&>(*B_RR_ptr)
		);
#endif

	MatrixHessianSuperBasic::initialize(
		n,n_R,i_x_free,i_x_fixed,bnd_fixed
		,B_RR_ptr,B_RX_ptr,B_RX_trans,B_XX_ptr 
		);
}

// Overridden from MatrixSymInitDiagonal

void MatrixHessianSuperBasicInitDiagonal::init_identity(
	size_type n, value_type alpha )
{
	assert_initialized();
	B_RR_init_->init_identity(n,alpha);
	MatrixHessianSuperBasic::initialize(
		n,n,NULL,NULL,NULL
		,this->B_RR_ptr()
		,this->B_RX_ptr(),this->B_RX_trans()
		,this->B_XX_ptr()
		);
}

void MatrixHessianSuperBasicInitDiagonal::init_diagonal(
	const VectorSlice& diag )
{
	assert_initialized();
	B_RR_init_->init_diagonal(diag);
	MatrixHessianSuperBasic::initialize(
		diag.size(),diag.size(),NULL,NULL,NULL
		,this->B_RR_ptr()
		,this->B_RX_ptr(),this->B_RX_trans()
		,this->B_XX_ptr()
		);
}

} // end namespace ConstrainedOptimizationPack
