// ////////////////////////////////////////////////////////////////////
// MatrixOpNonsingThyra.cpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.
//

#include "MatrixOpNonsingThyra.hpp"
#include "VectorMutableThyra.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

// Constructors / Initializers

MatrixOpNonsingThyra::MatrixOpNonsingThyra()
{}

MatrixOpNonsingThyra::MatrixOpNonsingThyra(
	const Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<value_type> >  &thyra_linear_op_ns
	,BLAS_Cpp::Transp                                                            thyra_linear_op_trans
	)
{
	this->initialize(thyra_linear_op_ns,thyra_linear_op_trans);
}

void MatrixOpNonsingThyra::initialize(
	const Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<value_type> >  &thyra_linear_op_ns
	,BLAS_Cpp::Transp                                                            thyra_linear_op_trans
	)
{
	namespace mmp = MemMngPack;
	TEST_FOR_EXCEPTION(
		thyra_linear_op_ns.get()==NULL, std::invalid_argument
		,"MatrixOpNonsingThyra::initialize(thyra_linear_op_ns): Error!"
		);
	MatrixOpThyra::initialize(thyra_linear_op_ns,thyra_linear_op_trans);
}

Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<value_type> > 
MatrixOpNonsingThyra::set_uninitialized()
{
	Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<value_type> >
		tmp_thyra_linear_op_ns = thyra_linear_op_ns();
	MatrixOpThyra::set_uninitialized();
	return tmp_thyra_linear_op_ns;
}

Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<value_type> >
MatrixOpNonsingThyra::thyra_linear_op_ns() const
{
	return Teuchos::rcp_dynamic_cast<const Thyra::LinearOpWithSolveBase<value_type> >(this->thyra_linear_op());
}

// Overridden from MatrixOp (needed to remove ambiguities)

MatrixOp::mat_mut_ptr_t
MatrixOpNonsingThyra::clone()
{
	return this->MatrixOpThyra::clone();
}

// Overridden from MatrixNonsing

void MatrixOpNonsingThyra::V_InvMtV(
	VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
	,const Vector& v_rhs2
	) const
{
	using Teuchos::dyn_cast;
	using BLAS_Cpp::trans_trans;
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		v_lhs==NULL, std::invalid_argument
		,"MatrixOpThyra::Vp_StMtV(...): Error!"
		);
#endif
	*v_lhs = 0.0; // Must initialize before sending to solve(...)!
	VectorMutableThyra &v_thyra_lhs = dyn_cast<VectorMutableThyra>(*v_lhs);
	Teuchos::RefCountPtr<Thyra::VectorBase<value_type> > thyra_vec_lhs = v_thyra_lhs.set_uninitialized();
	Thyra::solve(
    *thyra_linear_op_ns()
		,trans_trans(trans_rhs1,thyra_linear_op_trans())==BLAS_Cpp::no_trans ? Thyra::NOTRANS : Thyra::TRANS  // M_trans
		,*dyn_cast<const VectorMutableThyra>(v_rhs2).thyra_vec()                                              // y
		,thyra_vec_lhs.get()                                                                                  // x
		);
	v_thyra_lhs.initialize(thyra_vec_lhs);
}

void MatrixOpNonsingThyra::M_StInvMtM(
	MatrixOp* m_lhs, value_type alpha
	,BLAS_Cpp::Transp trans_rhs1
	,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
	) const
{
	MatrixNonsing::M_StInvMtM(m_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2); // ToDo: Specialize!
}

// Overridden from MatrixOpNonsing

MatrixOpNonsing::mat_mwons_ptr_t
MatrixOpNonsingThyra::clone_mwons() const
{
	return Teuchos::null; // ToDo: Add a clone function to Thyra::LinearOpWithSolveBase???
	//return Teuchos::rcp(new MatrixOpNonsingThyra(thyra_linear_op_ns()->clone_lows()));
}

} // end namespace AbstractLinAlgPack
