// ////////////////////////////////////////////////////////////////////////
// test_nlp_first_order_direct.cpp

#include <assert.h>

#include "test_nlp_first_order_direct.h"
#include "NLPTester.h"
#include "NLPTesterSetOptions.h"
#include "NLPFirstOrderDirectTester.h"
#include "NLPFirstOrderDirectTesterSetOptions.h"
#include "NLPInterfacePack/include/NLPFirstOrderDirect.h"
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

bool NLPInterfacePack::test_nlp_first_order_direct(
	NLPFirstOrderDirect*                          nlp
	,OptionsFromStreamPack::OptionsFromStream*    options
	,std::ostream*                                out
	)
{
	using TestingHelperPack::update_success;

	bool result;
	bool success = true;

	if(out)
		*out << "\n****************************************"
			 << "\n*** test_nlp_first_order_direct(...) ***"
			 << "\n****************************************\n";
	
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
	if(out)
		*out << "\nTesting nlp->space_x()->sub_space(nlp->var_dep()) ...\n";
	result = vec_space_tester.check_vector_space(
		*nlp->space_x()->sub_space(nlp->var_dep()),out);
	if(out) {
		if(result)
			*out << "nlp->space_x()->sub_space(nlp->var_dep()) checks out!\n";
		else
			*out << "nlp->space_x()->sub_space(nlp->var_dep()) check failed!\n";
	}
	update_success( result, &success );
	if(out)
		*out << "\nTesting nlp->space_x()->sub_space(nlp->var_indep()) ...\n";
	result = vec_space_tester.check_vector_space(
		*nlp->space_x()->sub_space(nlp->var_indep()),out);
	if(out) {
		if(result)
			*out << "nlp->space_x()->sub_space(nlp->var_indep()) checks out!\n";
		else
			*out << "nlp->space_x()->sub_space(nlp->var_indep()) check failed!\n";
	}
	update_success( result, &success );
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
	if(out)
		*out << "\nTesting nlp->space_c()->sub_space(nlp->con_decomp()) ...\n";
	result = vec_space_tester.check_vector_space(
		*nlp->space_c()->sub_space(nlp->con_decomp()),out);
	if(out) {
		if(result)
			*out << "nlp->space_c()->sub_space(nlp->con_decomp()) checks out!\n";
		else
			*out << "nlp->space_c()->sub_space(nlp->con_decomp()) check failed!\n";
	}
	update_success( result, &success );
	if( nlp->con_decomp().size() < nlp->m() ) {
		if(out)
			*out << "\nTesting nlp->space_c()->sub_space(nlp->con_undecomp()) ...\n";
		result = vec_space_tester.check_vector_space(
			*nlp->space_c()->sub_space(nlp->con_undecomp()),out);
		if(out) {
			if(result)
				*out << "nlp->space_c()->sub_space(nlp->con_undecomp()) checks out!\n";
			else
				*out << "nlp->space_c()->sub_space(nlp->con_undecomp()) check failed!\n";
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
	
	// Test the NLPFirstOrderDirect interface now!

	if(out)
		*out << "\nCalling nlp->calc_point(...) at nlp->xinit() ...\n";
	const size_type
		n  = nlp->n(),
		m  = nlp->m(),
		mI = nlp->mI();
	const Range1D
		var_dep      = nlp->var_dep(),
		var_indep    = nlp->var_indep(),
		con_decomp   = nlp->con_decomp(),
		con_undecomp = nlp->con_undecomp();
	VectorSpace::vec_mut_ptr_t
		c   =      nlp->space_c()->create_member(),
		h   = mI ? nlp->space_h()->create_member() : NULL,
		Gf  =      nlp->space_x()->create_member(),
		py  =      nlp->space_x()->sub_space(var_dep)->create_member(),
		rGf =      nlp->space_x()->sub_space(var_indep)->create_member();
	NLPFirstOrderDirect::mat_space_ptr_t::element_type::mat_ptr_t
		GcU = con_decomp.size() < m ? nlp->space_GcU()->create_member() : NULL,
		Gh  = mI                    ? nlp->space_Gh()->create_member()  : NULL,
		D   =                         nlp->space_D()->create_member(),
		V   = con_decomp.size() < m ? nlp->space_V()->create_member()   : NULL,
		P   = mI                    ? nlp->space_P()->create_member()   : NULL;
	nlp->calc_point(
		nlp->xinit(), NULL, c.get(), true, h.get(), Gf.get(), py.get(), rGf.get()
		,GcU.get(), Gh.get(), D.get(), V.get(), P.get() );
	if(out) {
		*out << "\n||Gf||inf  = " << Gf->norm_inf();
		if(nlp_tester.print_all())
			*out << "\nGf =\n" << *Gf;
		*out << "\n||py||inf  = " << py->norm_inf();
		if(nlp_tester.print_all())
			*out << "\npy =\n" << *py;
		*out << "\n||rGf||inf  = " << rGf->norm_inf();
		if(nlp_tester.print_all())
			*out << "\nrGf =\n" << *rGf;
		if(nlp_tester.print_all())
			*out << "\nD =\n" << *D;
		if( con_decomp.size() < m ) {
			assert(0); // ToDo: Print GcU and V
		}
		if( mI ) {
			assert(0); // ToDo: Print Gh and P
		}
	}

	NLPFirstOrderDirectTester nlp_first_order_direct_tester;
	if(options) {
		NLPFirstOrderDirectTesterSetOptions
			nlp_tester_opt_setter(&nlp_first_order_direct_tester);
		nlp_tester_opt_setter.set_options(*options);
	}
	// ToDo: Set options from stream!
	result = nlp_first_order_direct_tester.finite_diff_check(
		nlp, nlp->xinit(), nlp->num_bounded_x() ? &nlp->xl() : NULL
		,nlp->num_bounded_x() ? &nlp->xu() : NULL, 0.0, c.get(), h.get()
		,Gf.get(),py.get(),rGf.get(),GcU.get(),Gh.get(),D.get(),V.get(),P.get()
		,print_all_warnings, out );
	update_success( result, &success );

	return success;
}
