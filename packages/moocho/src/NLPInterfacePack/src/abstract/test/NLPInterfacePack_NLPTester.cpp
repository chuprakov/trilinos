// ///////////////////////////////////////////////////////////
// NLPTester.cpp
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

#include <iomanip>
#include <ostream>

#include "NLPTester.hpp"
#include "NLPInterfacePack/src/NLP.hpp"
#include "AbstractLinAlgPack/src/VectorSpace.hpp"
#include "AbstractLinAlgPack/src/VectorWithOpMutable.hpp"
#include "AbstractLinAlgPack/src/VectorWithOpOut.hpp"
#include "AbstractLinAlgPack/src/VectorStdOps.hpp"
#include "AbstractLinAlgPack/src/VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack/src/assert_print_nan_inf.hpp"
#include "update_success.hpp"

namespace NLPInterfacePack {

NLPTester::NLPTester(
	bool     print_all
	,bool    throw_exception
	)
	:print_all_(print_all), throw_exception_(throw_exception)
{}

bool NLPTester::test_interface(
	NLP                           *nlp
	,const VectorWithOp           &xo
	,bool                         print_all_warnings
	,std::ostream                 *out
	)
{
	using TestingHelperPack::update_success;
	using AbstractLinAlgPack::assert_print_nan_inf;

	bool result;
	bool success = true;

	if(out) {
		*out << std::boolalpha
			 << std::endl
			 << "**************************************\n"
			 << "*** NLPTester::test_interface(...) ***\n"
			 << "**************************************\n";
	}

	try {
	
		// Check the input x
		assert_print_nan_inf(xo,"xo",true,out); 

		// Initialize the NLP if it has not been already and force in bounds
		if(out)
			*out << "\nnlp->force_xinit_in_bounds(true)";
		nlp->force_xinit_in_bounds();
		if(out)
			*out << "\nnlp->initialize(true)\n";
		nlp->initialize(true);
		
		const size_type
			n = nlp->n(),
			m = nlp->m(),
			mI = nlp->mI();
		if(out)
			*out << "\n*** Dimensions of the NLP"
				 << "\nnlp->n()  = " << n
				 << "\nnlp->m()  = " << m
				 << "\nnlp->mI() = " << mI
				 << std::endl;
		if( n < m ) {
			if(*out)
				*out << "Error! n = " << n << " < m = " << m << " is not allowed!\n";
			THROW_EXCEPTION(
				throw_exception_, std::logic_error
				,"NLPTester::test_interface(...): Error! n = " << n << " < m = " << m << " is not allowed!"
				);
		}

		// Validate the vector spaces
		if(out)
			*out << "\n*** Validate the dimensions of the vector spaces";
		
		result = nlp->space_x()->dim() == nlp->n();
		update_success( result, &success );
		if(out)
			*out << "\ncheck: nlp->space_x()->dim() = " << nlp->space_x()->dim()
				 << " == nlp->n() = " << nlp->n() << ": " << result;

		if( nlp->m() ) {
			result = nlp->space_c()->dim() == nlp->m();
			update_success( result, &success );
			if(out)
				*out << "\ncheck: nlp->space_c()->dim() = " << nlp->space_c()->dim()
					 << " == nlp->m() = " << nlp->m() << ": " << result;
		}
		else {
			result = nlp->space_c().get() == NULL;
			update_success( result, &success );
			if(out)
				*out << "\ncheck: nlp->space_c().get() = " << nlp->space_c().get()
					 << " == NULL: " << result;
		}

		if( nlp->mI() ) {
			result = nlp->space_h()->dim() == nlp->mI();
			update_success( result, &success );
			if(out)
				*out << "\ncheck: nlp->space_h()->dim() = " << nlp->space_h()->dim()
					 << " == nlp->mI() = " << nlp->mI() << ": " << result;
		}
		else {
			result = nlp->space_h().get() == NULL;
			update_success( result, &success );
			if(out)
				*out << "\ncheck: nlp->space_h().get() = " << nlp->space_h().get()
					 << " == NULL: " << result;
		}

		// Validate the bounds on the unknowns.

		const VectorWithOp &xinit = nlp->xinit();
		if(out)
			*out << "\n||nlp->xinit()||inf = " << xinit.norm_inf() << std::endl;
		if(out && print_all())
			*out << "\nnlp->xinit() =\n" << xinit;

		assert_print_nan_inf(xinit,"xinit",true,out); 

		if(out)
			*out << "\n*** Validate that the initial starting point is in bounds ...\n";
		const VectorWithOp
			&xl = nlp->xl(),
			&xu = nlp->xu();
		if(out && print_all())
			*out << "\nnlp->xl() =\n" << xl
				 << "\nnlp->xu() =\n" << xu;
		assert_print_nan_inf(xl,"xl",true,out); 
		assert_print_nan_inf(xu,"xu",true,out); 

		// Validate that xl <= xinit <= xu.
		VectorSpace::vec_mut_ptr_t
			d = nlp->space_x()->create_member();
		*d = 1.0;
		std::pair<value_type,value_type>
			u = AbstractLinAlgPack::max_near_feas_step(
				xinit, *d, nlp->xl(), nlp->xu(), 0.0
				);
		result = u.first >= 0.0;
		update_success( result, &success );
		if(out) {
			*out << "\ncheck: xl <= x <= xu : " << result;
			if(result)
				*out << "\nxinit is in bounds with { max |u| | xl <= x + u <= xu } -> "
					 << ( u.first > -u.second ? u.first : u.second  ) << std::endl;
		}
			
		size_type 
			num_bounded_x = AbstractLinAlgPack::num_bounded(
				nlp->xl(), nlp->xu(), NLP::infinite_bound() );
		result = (num_bounded_x == nlp->num_bounded_x());
		update_success( result, &success );
		if(out)
			*out << "\ncheck: num_bounded(nlp->xl(),nlp->xu()) = " << num_bounded_x
				 << " == nlp->num_bounded_x() = " << nlp->num_bounded_x()
				 << ": " << result << std::endl;

		// Validate bounds on the general inequalities
		if( nlp->mI() ) {
			
			assert(0); // ToDo: Implement this once we have an NLP with general inequalities

		}

		// Get the initial Lagrange multipliers
		if(out)
			*out << "\nGetting the initial estimates for the Lagrange mutipliers ...\n";
		VectorSpace::vec_mut_ptr_t
			lambda, lambdaI, nu;
		nlp->get_init_lagrange_mult(
			(  nlp->m()
			   ? (lambda  = nlp->space_c()->create_member()).get() 
			   : (VectorWithOpMutable*)NULL )
			,( nlp->mI()
			   ? (lambdaI = nlp->space_h()->create_member()).get()
			   : (VectorWithOpMutable*)NULL )
			,( nlp->num_bounded_x()
			   ? (nu = nlp->space_x()->create_member()).get()
			   : (VectorWithOpMutable*)NULL )
			);

		if(out) {
			if(lambda.get())
				*out << "\n||lambda||inf  = " << lambda->norm_inf();
			if(lambdaI.get())
				*out << "\n||lambdaI||inf = " << lambda->norm_inf()
					 << "\nlambdaI.nz()   = " << lambda->nz();
			if(nu.get())
				*out << "\n||nu||inf      = " << nu->norm_inf()
					 << "\nnu.nz()        = " << nu->nz();
			*out << std::endl;
			if(print_all()) {
				if(lambda.get())
					*out << "\nlambda =\n" << *lambda;
				if(lambdaI.get())
					*out << "\nlambdaI =\n" << *lambdaI;
				if(nu.get())
					*out << "\nnu =\n" << *nu;
			}
		}
		if(lambda.get())
			assert_print_nan_inf(*lambda,"lambda",true,out); 
		if(lambdaI.get())
			assert_print_nan_inf(*lambdaI,"lambdaI",true,out); 
		if(nu.get())
			assert_print_nan_inf(*nu,"nu",true,out); 

		// Save the current reference that are set to be set back at the end
		value_type            *f_saved = NULL;
		VectorWithOpMutable   *c_saved = NULL;
		VectorWithOpMutable   *h_saved = NULL;
		f_saved = nlp->get_f();
		if( nlp->m() )  c_saved = nlp->get_c();
		if( nlp->mI() ) h_saved = nlp->get_h();

		// Create calcualtion quantities
		value_type                   f;
		VectorSpace::vec_mut_ptr_t   c;
		VectorSpace::vec_mut_ptr_t   h;
		if( nlp->m() )
			c = nlp->space_c()->create_member();
		if( nlp->mI() )
			h = nlp->space_h()->create_member();

		// Set the calculation quantities
		nlp->set_f(&f);
		if( nlp->m() )  nlp->set_c(c.get());
		if( nlp->mI() ) nlp->set_h(h.get());

		// Calculate the quantities at xo
		if(out)
			*out << "\n*** Evaluate the point xo ...\n";
		nlp->set_multi_calc(true);
		nlp->calc_f(xo,true);
		if(nlp->m())  nlp->calc_c(xo,false);
		if(nlp->mI()) nlp->calc_h(xo,false);

		if(out) {
			*out << "\nf(xo) = " << f;
			if(nlp->m())
				*out << "\n||c(xo)||inf = " << nlp->c().norm_inf();
			if(nlp->mI())
				*out << "\n||h(xo)||inf = " << nlp->h().norm_inf();
			*out << std::endl;
			if(print_all()) {
				if(nlp->m())
					*out << "\nc(xo) =\n" << nlp->c();
				if(nlp->mI())
					*out << "\nh(xo) =\n" << nlp->h();
			}
		}

		if(c.get())
			assert_print_nan_inf(*c,"c(xo)",true,out); 
		if(h.get())
			assert_print_nan_inf(*h,"h(xo)",true,out); 

		// Report the final solution!
		if(out)
			*out << "\n*** Report this point to the NLP as suboptimal ...\n";
		nlp->report_final_solution(	xo, lambda.get(), lambdaI.get(), nu.get(), false );

		// Print the number of evaluations!
		if(out) {
			*out << "\n*** Print the number of evaluations ...\n";
			*out << "\nnlp->num_f_evals() = " << nlp->num_f_evals();
			if(nlp->m())
				*out << "\nnlp->num_c_evals() = " << nlp->num_c_evals();
			if(nlp->mI())
				*out << "\nnlp->num_h_evals() = " << nlp->num_h_evals();
			*out << std::endl;
		}

		// Set the original quantities back
		nlp->set_f(f_saved);
		if(nlp->m())  nlp->set_c(c_saved);
		if(nlp->mI()) nlp->set_h(h_saved);

	}
	catch(const std::exception& except) {
		if(out)
			*out << "Caught a std::exception: " << except.what() << std::endl;
		success = false;
		if(throw_exception())
			throw;
	}
	catch(...) {
		if(out)
			*out << "Caught an unknown exception!\n";
		success = false;
		if(throw_exception())
			throw;
	}

	return success;
}

} // namespace NLPInterfacePack
