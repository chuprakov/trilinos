// ////////////////////////////////////////////////////////////////////////////
// MeritFunc_PenaltyParamUpdateGuts_AddedStep.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>
#include <typeinfo>

#include "../../include/std/MeritFunc_PenaltyParamUpdateGuts_AddedStep.h"
#include "../../include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/MeritFuncPenaltyParam.h"
#include "ConstrainedOptimizationPack/include/MeritFuncNLPDirecDeriv.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOut.h"

namespace {

typedef ReducedSpaceSQPPack::value_type value_type;
inline value_type max(value_type v1, value_type v2)
{	return (v1 > v2) ? v1 : v2; }

}

namespace ReducedSpaceSQPPack {

MeritFunc_PenaltyParamUpdateGuts_AddedStep::MeritFunc_PenaltyParamUpdateGuts_AddedStep(
		  const merit_func_ptr_t& merit_func, value_type small_mu
		, value_type mult_factor, value_type kkt_near_sol )
	: merit_func_(merit_func), near_solution_(false)
		, small_mu_(small_mu), mult_factor_(mult_factor), kkt_near_sol_(kkt_near_sol)
{}

bool MeritFunc_PenaltyParamUpdateGuts_AddedStep::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type
	, poss_type assoc_step_poss)
{
	using LinAlgPack::norm_inf;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();
	
	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	MeritFuncPenaltyParam
		*param = dynamic_cast<MeritFuncPenaltyParam*>(&merit_func());
	if( !param ) {
		std::ostringstream omsg;
		omsg
			<< "MeritFunc_PenaltyParamUpdateGuts_AddedStep::do_step(...), Error "
			<< "The class " << typeid(&merit_func()).name() << " does not support the "
			<< "MeritFuncPenaltyParam iterface\n";
		out << omsg.str();
		throw std::logic_error( omsg.str() );
	}

	MeritFuncNLPDirecDeriv
		*direc_deriv = dynamic_cast<MeritFuncNLPDirecDeriv*>(&merit_func());
	if( !direc_deriv ) {
		std::ostringstream omsg;
		omsg
			<< "MeritFunc_PenaltyParamUpdateGuts_AddedStep::do_step(...), Error "
			<< "The class " << typeid(&merit_func()).name() << " does not support the "
			<< "MeritFuncNLPDirecDeriv iterface\n";
		out << omsg.str();
		throw std::logic_error( omsg.str() );
	}

	value_type	new_mu = 0.0;

	if( s.mu().updated_k(0) ) {
		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			out << "\nmu_k is already updated by someone else so we will use it.\n";
		}
		new_mu = s.mu().get_k(0);
	}
	else {
		value_type min_mu = 0.0;
		if ( this->min_mu(s,&min_mu) ) {

			// Update the penalty parameter as defined in the fortran rSQP code (EXACT2())
		
			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
				out << "\nUpdate the penalty parameter...\n";
			}

			value_type
					mu_km1 = param->mu(),
					mult_fact = (1.0 + mult_factor_);

			if(near_solution_) {
				if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
					out << "\nNear solution, forcing mu_k >= mu_km1...\n";
				}
				new_mu = max( max( mu_km1, mult_fact * min_mu ), small_mu_ );
			}
			else {
				if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
					out << "\nNot near solution, allowing reduction in mu ...\n";
				}
				new_mu =	max(
								  (3.0 * mu_km1 + min_mu) / 4.0	
								, max( mult_fact * min_mu, small_mu_ )
								); 
				value_type kkt_error = s.opt_kkt_err().get_k(0) + s.feas_kkt_err().get_k(0);
				
				if(kkt_error <= kkt_near_sol_) {
					if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
						out << "\nkkt_error = " << kkt_error << " <= kkt_near_sol = "
								<< kkt_near_sol_ << std::endl
							<< "Switching to forcing mu_k >= mu_km1 in the future\n";
					}
					near_solution_ = true;
				}
			}

			s.mu().set_k(0) = new_mu;
		}
		else {
			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
				out << "\nDon't have the info to update penalty parameter so just use the last updated...\n";
			}
			new_mu = param->mu();
		}
	}

	// Set the penalty parameter
	param->mu( new_mu );

	// In addition also compute the directional derivative
	direc_deriv->calc_deriv( s.Gf().get_k(0)(), s.c().get_k(0)(), s.d().get_k(0)() );

	if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
		out << "\nmu = " << s.mu().get_k(0) << "\n";
	}

	return true;
}

void MeritFunc_PenaltyParamUpdateGuts_AddedStep::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Update the penalty parameter for the merit function to ensure\n"
		<< L << "*** a descent direction a directional derivatieve.\n"
		<< L << "*** phi is a merit function object that uses the penalty parameter mu.\n"
		<< L << "default: near_solution = false\n"
		<< L << "         small_mu = " << small_mu_ << std::endl
		<< L << "         mult_factor = " << mult_factor_ << std::endl
		<< L << "         kkt_near_sol = " << kkt_near_sol_ << std::endl
		<< L << "if mu_k is already updated then\n"
		<< L << "    *** Just use mu_k already computed.\n"
		<< L << "else\n"
		;
	            print_min_mu_step( out, L + "    " ); 
	out
		<< L << "    if update_mu == true then\n"
		<< L << "        mu_last = phi.mu()\n"
		<< L << "        mult_fact = 1.0 + mult_factor\n"
		<< L << "        if near_solution == true\n"
		<< L << "            mu_k = max( max( mu_last, mult_fact*min_mu ), small_mu )\n"
		<< L << "        else\n"
		<< L << "            mu_k = max(   ( 3.0 * mu_last + min_mu ) / 4.0\n"
		<< L << "                        , max( mult_fact * min_mu , small_mu )     )\n"
		<< L << "            kkt_error = opt_kkt_err_k + feas_kkt_err_k\n"
		<< L << "            if kkt_error <= kkt_near_sol then\n"
		<< L << "                near_solution = true\n"
		<< L << "            end\n"
		<< L << "        end\n"
		<< L << "    else\n"
		<< L << "        mu_k = phi.mu()\n"
		<< L << "    end\n"
		<< L << "end\n"
		<< L << "phi.mu(mu_k)\n"
		<< L << "phi.calc_deriv(Gf_k,c_k,d_k)\n";
}

// Overridden from MeritFunc_PenaltyParamUpdate_AddedStep

void MeritFunc_PenaltyParamUpdateGuts_AddedStep::small_mu( value_type small_mu )
{
	small_mu_ = small_mu;
}

value_type MeritFunc_PenaltyParamUpdateGuts_AddedStep::small_mu() const
{
	return small_mu_;
}

void MeritFunc_PenaltyParamUpdateGuts_AddedStep::mult_factor( value_type mult_factor )
{
	mult_factor_ = mult_factor;
}

value_type MeritFunc_PenaltyParamUpdateGuts_AddedStep::mult_factor() const
{
	return mult_factor_;
}

void MeritFunc_PenaltyParamUpdateGuts_AddedStep::kkt_near_sol( value_type kkt_near_sol )
{
	kkt_near_sol_ = kkt_near_sol;
}

value_type MeritFunc_PenaltyParamUpdateGuts_AddedStep::kkt_near_sol() const
{
	return kkt_near_sol_;
}


}	// end namespace ReducedSpaceSQPPack