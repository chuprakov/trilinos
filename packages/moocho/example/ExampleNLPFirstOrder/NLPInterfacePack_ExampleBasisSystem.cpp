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

#include "ExampleBasisSystem.h"
#include "AbstractLinAlgPack/include/MatrixSymDiagonalStd.h"
#include "AbstractLinAlgPack/include/MatrixCompositeStd.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "RTOpPack/include/RTOpCppC.h"
#include "AbstractFactoryStd.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"

namespace NLPInterfacePack {
 
ExampleBasisSystem::ExampleBasisSystem(
	const VectorSpace::space_ptr_t       &space_x
	,const Range1D                       &var_dep
	,const Range1D                       &var_indep
	)
	:BasisSystemCompositeStd(
		space_x
		,var_dep
		,var_indep
		,space_x->sub_space(var_dep)
		,MemMngPack::rcp(
			new MemMngPack::AbstractFactoryStd<MatrixWithOpNonsingular,MatrixSymDiagonalStd>())       // C
		,MemMngPack::rcp(
			new MemMngPack::AbstractFactoryStd<MatrixSymWithOp,MatrixSymDiagonalStd>())               // D'*D
		,MemMngPack::rcp(
			new MemMngPack::AbstractFactoryStd<MatrixSymWithOpNonsingular,MatrixSymDiagonalStd>())    // S
		,MemMngPack::rcp(
			new MemMngPack::AbstractFactoryStd<MatrixWithOp,MatrixSymDiagonalStd>())                  // D
		)
{}
	
void ExampleBasisSystem::initialize(
	const VectorSpace::space_ptr_t       &space_x
	,const Range1D                       &var_dep
	,const Range1D                       &var_indep
	)
{
	namespace mmp = MemMngPack;
	THROW_EXCEPTION(
		space_x.get() == NULL, std::invalid_argument
		,"ExampleBasisSystem::initialize(...) : Error, space_x must be specified!"
		);
	BasisSystemCompositeStd::initialize(
		space_x
		,var_dep
		,var_indep
		,space_x->sub_space(var_dep)
		,mmp::rcp(new mmp::AbstractFactoryStd<MatrixWithOpNonsingular,MatrixSymDiagonalStd>())      // C
		,mmp::rcp(new mmp::AbstractFactoryStd<MatrixSymWithOp,MatrixSymDiagonalStd>())              // D'*D
		,mmp::rcp(new mmp::AbstractFactoryStd<MatrixSymWithOpNonsingular,MatrixSymDiagonalStd>())   // S
		,mmp::rcp(new mmp::AbstractFactoryStd<MatrixWithOp,MatrixSymDiagonalStd>())                 // D
		);
}

void ExampleBasisSystem::update_D(
	const MatrixWithOpNonsingular&  C
	,const MatrixWithOp&            N
	,MatrixWithOp*                  D
	,EMatRelations                  mat_rel
	) const
{
	using DynamicCastHelperPack::dyn_cast;

	THROW_EXCEPTION(
		D == NULL, std::logic_error
		,"ExampleBasisSystem::update_D(...): Error!" );

	const MatrixSymDiagonalStd
		&C_aggr = dyn_cast<const MatrixSymDiagonalStd>(C),
		&N_aggr = dyn_cast<const MatrixSymDiagonalStd>(N);
	MatrixSymDiagonalStd
		&D_sym_diag = dyn_cast<MatrixSymDiagonalStd>(*D);
	if( D_sym_diag.rows() != C.rows() )
		D_sym_diag.initialize(
			this->space_x()->sub_space(this->var_dep())->create_member()
			);
	AbstractLinAlgPack::ele_wise_divide(                           // D_diag = - N_diag ./ C_diag
		-1.0, N_aggr.diag(), C_aggr.diag(), &D_sym_diag.diag() );  // ...
}

} // end namespace NLPInterfacePack
