// ////////////////////////////////////////////////////////////////////////
// test_nlp_first_order_info.cpp
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

#include <assert.h>

#include "test_nlp_first_order_info.h"
#include "NLPInterfacePack/include/CalcFiniteDiffProd.h"
#include "NLPInterfacePack/include/CalcFiniteDiffProdSetOptions.h"
#include "NLPTester.h"
#include "NLPTesterSetOptions.h"
#include "NLPFirstDerivativesTester.h"
#include "NLPFirstDerivativesTesterSetOptions.h"
#include "NLPInterfacePack/include/NLPFirstOrderInfo.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "AbstractLinAlgPack/include/VectorSpaceTester.h"
#include "AbstractLinAlgPack/include/VectorSpaceTesterSetOptions.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorWithOpOut.h"
#include "AbstractLinAlgPack/include/MatrixSpace.h"
#include "AbstractLinAlgPack/include/MatrixWithOp.h"
#include "AbstractLinAlgPack/include/MatrixWithOpOut.h"
#include "Range1D.h"
#include "update_success.h"

bool NLPInterfacePack::test_nlp_first_order_info(
	NLPFirstOrderInfo*                          nlp
	,OptionsFromStreamPack::OptionsFromStream*    options
	,std::ostream*                                out
	)
{
	namespace rcp = ReferenceCountingPack;
	using TestingHelperPack::update_success;

	bool result;
	bool success = true;

	if(out)
		*out << "\n**************************************"
			 << "\n*** test_nlp_first_order_info(...) ***"
			 << "\n**************************************\n";
	
	nlp->initialize();

	// Test the Vector spaces
	if(out)
		*out << "\nTesting the vector spaces ...\n";
	
	VectorSpaceTester vec_space_tester;
	if(options) {
		VectorSpaceTesterSetOptions
			opt_setter(&vec_space_tester);
		opt_setter.set_options(*options);
	}

	if(out)
		*out << "\nTesting nlp->space_x() ...\n";
	result = vec_space_tester.check_vector_space(*nlp->space_x(),out);
	if(out) {
		if(result)
			*out << "nlp->space_x() checks out!\n";
		else
			*out << "nlp->space_x() check failed!\n";
	}
	update_success( result, &success );

	if( nlp->m() ) {
		if(out)
			*out << "\nTesting nlp->space_c() ...\n";
		result = vec_space_tester.check_vector_space(*nlp->space_c(),out);
		if(out) {
			if(result)
				*out << "nlp->space_c() checks out!\n";
			else
				*out << "nlp->space_c() check failed!\n";
		}
		update_success( result, &success );
	}

	if( nlp->mI() ) {
		if(out)
			*out << "\nTesting nlp->space_h() ...\n";
		result = vec_space_tester.check_vector_space(*nlp->space_h(),out);
		if(out) {
			if(result)
				*out << "nlp->space_h() checks out!\n";
			else
				*out << "nlp->space_h() check failed!\n";
		}
		update_success( result, &success );
	}
	
	// Test the NLP interface first!

	NLPTester nlp_tester;
	if(options) {
		NLPTesterSetOptions
			nlp_tester_opt_setter(&nlp_tester);
		nlp_tester_opt_setter.set_options(*options);
	}
	const bool print_all_warnings = nlp_tester.print_all();

	result = nlp_tester.test_interface(
		nlp, nlp->xinit(), print_all_warnings, out );
	update_success( result, &success );
	
	// Test the NLPFirstOrderInfo interface now!

	const size_type
		n  = nlp->n(),
		m  = nlp->m(),
		mI = nlp->mI();
	VectorSpace::vec_mut_ptr_t
		c   = m  ? nlp->space_c()->create_member() : rcp::null,
		h   = mI ? nlp->space_h()->create_member() : rcp::null,
		Gf  =      nlp->space_x()->create_member();
	NLPFirstOrderInfo::mat_space_ptr_t::element_type::mat_ptr_t
		Gc  = m  ? nlp->space_Gc()->create_member() : rcp::null,
		Gh  = mI ? nlp->space_Gh()->create_member() : rcp::null;

	if(m) {
		if(out)
			*out << "\nCalling nlp->calc_Gc(...) at nlp->xinit() ...\n";
		nlp->set_Gc( Gc.get() );
		nlp->calc_Gc( nlp->xinit(), true );
		if(nlp_tester.print_all())
			*out << "\nGc =\n" << *Gc;
	}
	if(mI) {
		if(out)
			*out << "\nCalling nlp->calc_Gh(...) at nlp->xinit() ...\n";
		nlp->set_Gh( Gh.get() );
		nlp->calc_Gh( nlp->xinit(), m == 0 );
		if(nlp_tester.print_all())
			*out << "\nGh =\n" << *Gh;
	}

	if(out)
		*out << "\nCalling nlp->calc_Gf(...) at nlp->xinit() ...\n";
	nlp->set_Gf( Gf.get() );
	nlp->calc_Gf( nlp->xinit(), m == 0 && mI == 0 );
	if(nlp_tester.print_all())
		*out << "\nGf =\n" << *Gf;

	CalcFiniteDiffProd
		calc_fd_prod;
	if(options) {
		CalcFiniteDiffProdSetOptions
			options_setter( &calc_fd_prod );
		options_setter.set_options(*options);
	}
	NLPFirstDerivativesTester
		nlp_first_derivatives_tester(rcp::rcp(&calc_fd_prod,false));
	if(options) {
		NLPFirstDerivativesTesterSetOptions
			nlp_tester_opt_setter(&nlp_first_derivatives_tester);
		nlp_tester_opt_setter.set_options(*options);
	}
	result = nlp_first_derivatives_tester.finite_diff_check(
		nlp, nlp->xinit()
		,nlp->num_bounded_x() ? &nlp->xl() : NULL
		,nlp->num_bounded_x() ? &nlp->xu() : NULL
		,Gc.get(), Gh.get(), Gf.get()
		,print_all_warnings, out
		);
	update_success( result, &success );

	return success;
}
