// ///////////////////////////////////////////////////////////////////
// ExampleBasisSystem.cpp
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

#include "NLPInterfacePack_ExampleBasisSystem.hpp"
#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack_MatrixComposite.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_TestForException.hpp"

namespace NLPInterfacePack {
 
ExampleBasisSystem::ExampleBasisSystem(
	const VectorSpace::space_ptr_t       &space_x
	,const Range1D                       &var_dep
	,const Range1D                       &var_indep
	)
	:BasisSystemComposite(
		space_x
		,var_dep
		,var_indep
		,space_x->sub_space(var_dep)
		,Teuchos::rcp(
			new Teuchos::AbstractFactoryStd<MatrixOpNonsing,MatrixSymDiagStd>())       // C
		,Teuchos::rcp(
			new Teuchos::AbstractFactoryStd<MatrixSymOp,MatrixSymDiagStd>())           // D'*D
		,Teuchos::rcp(
			new Teuchos::AbstractFactoryStd<MatrixSymOpNonsing,MatrixSymDiagStd>())    // S
		,Teuchos::rcp(
			new Teuchos::AbstractFactoryStd<MatrixOp,MatrixSymDiagStd>())              // D
		)
{}
	
void ExampleBasisSystem::initialize(
	const VectorSpace::space_ptr_t       &space_x
	,const Range1D                       &var_dep
	,const Range1D                       &var_indep
	)
{
	namespace mmp = MemMngPack;
	TEST_FOR_EXCEPTION(
		space_x.get() == NULL, std::invalid_argument
		,"ExampleBasisSystem::initialize(...) : Error, space_x must be specified!"
		);
	BasisSystemComposite::initialize(
		space_x
		,var_dep
		,var_indep
		,space_x->sub_space(var_dep)
		,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixOpNonsing,MatrixSymDiagStd>())      // C
		,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixSymOp,MatrixSymDiagStd>())          // D'*D
		,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixSymOpNonsing,MatrixSymDiagStd>())   // S
		,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixOp,MatrixSymDiagStd>())             // D
		);
}

void ExampleBasisSystem::update_D(
	const MatrixOpNonsing   &C
	,const MatrixOp         &N
	,MatrixOp               *D
	,EMatRelations          mat_rel
	) const
{
	using Teuchos::dyn_cast;

	TEST_FOR_EXCEPTION(
		D == NULL, std::logic_error
		,"ExampleBasisSystem::update_D(...): Error!" );

	const MatrixSymDiagStd
		&C_aggr = dyn_cast<const MatrixSymDiagStd>(C),
		&N_aggr = dyn_cast<const MatrixSymDiagStd>(N);
	MatrixSymDiagStd
		&D_sym_diag = dyn_cast<MatrixSymDiagStd>(*D);
	if( D_sym_diag.rows() != C.rows() )
		D_sym_diag.initialize(
			this->space_x()->sub_space(this->var_dep())->create_member()
			);
	AbstractLinAlgPack::ele_wise_divide(                           // D_diag = - N_diag ./ C_diag
		-1.0, N_aggr.diag(), C_aggr.diag(), &D_sym_diag.diag() );  // ...
}

} // end namespace NLPInterfacePack
