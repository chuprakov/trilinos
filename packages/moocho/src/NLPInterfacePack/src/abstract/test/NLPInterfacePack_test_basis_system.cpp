// /////////////////////////////////////////////////////////////////////
// test_basis_system.cpp
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

#include "test_basis_system.h"
#include "NLPInterfacePack/src/NLPFirstOrderInfo.h"
#include "AbstractLinAlgPack/src/BasisSystem.h"
#include "AbstractLinAlgPack/src/BasisSystemTester.h"
#include "AbstractLinAlgPack/src/BasisSystemTesterSetOptions.h"
#include "AbstractLinAlgPack/src/MatrixWithOpNonsingular.h"

bool NLPInterfacePack::test_basis_system(
 	NLPFirstOrderInfo*                             nlp
	,BasisSystem*                                  basis_sys
	,OptionsFromStreamPack::OptionsFromStream*     options
	,std::ostream*                                 out
	)
{
	namespace rcp = MemMngPack;

	const index_type
		n = nlp->n(),
		m = nlp->m(),
		mI = nlp->mI();

	// Create the matrices Gc and Gh
	NLPFirstOrderInfo::mat_fcty_ptr_t::element_type::obj_ptr_t
		Gc = ( m  ? nlp->factory_Gc()->create() : rcp::null ),
		Gh = ( mI ? nlp->factory_Gh()->create() : rcp::null );
	
	// Compute the matrices at xinit
	const VectorWithOp
		&xo = nlp->xinit();
	if(m)
		nlp->set_Gc(Gc.get());
	if(mI)
		nlp->set_Gh(Gh.get());
	if(m)
		nlp->calc_Gc(xo);
	if(mI)
		nlp->calc_Gh(xo,m==0);

	// Create the matrices C and D
	BasisSystem::mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t
		C = ( m ? basis_sys->factory_C()->create() : rcp::null);
	BasisSystem::mat_fcty_ptr_t::element_type::obj_ptr_t
		D = ( m && n > m && basis_sys->factory_C().get() ? basis_sys->factory_C()->create() : rcp::null);
	BasisSystem::mat_fcty_ptr_t::element_type::obj_ptr_t
		GcUP = ( m && n > m && basis_sys->factory_GcUP().get()  ? basis_sys->factory_GcUP()->create() : rcp::null);
	BasisSystem::mat_fcty_ptr_t::element_type::obj_ptr_t
		GhUP = ( mI && n > m && basis_sys->factory_GhUP().get() ? basis_sys->factory_GhUP()->create() : rcp::null);

	// Initialize C and D with basis_sys
	basis_sys->update_basis(
		Gc.get()
		,Gh.get()
		,C.get()
		,D.get()
		,GcUP.get()
		,GhUP.get()
		);

	// Test the basis and basis system objects.
	BasisSystemTester
		basis_sys_tester;
	if(options) {
		BasisSystemTesterSetOptions
			opt_setter(&basis_sys_tester);
		opt_setter.set_options(*options);
	}
	const bool result = basis_sys_tester.test_basis_system(
		*basis_sys
		,Gc.get()
		,Gh.get()
		,C.get()
		,NULL    // Create the N matrix internally
		,D.get()
		,GcUP.get()
		,GhUP.get()
		,out
		);

	return result;
}
