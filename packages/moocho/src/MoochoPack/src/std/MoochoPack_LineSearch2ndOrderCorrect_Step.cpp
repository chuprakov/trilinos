// ////////////////////////////////////////////////////////////////////////////
// LineSearch2ndOrderCorrect_Step.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>
#include <typeinfo>

#include "../../include/std/LineSearch2ndOrderCorrect_Step.h"
#include "../../include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/print_vector_change_stats.h"
#include "ConstrainedOptimizationPack/include/MeritFuncCalc1DQuadratic.h"
#include "ConstrainedOptimizationPack/include/MeritFuncCalcNLP.h"
#include "ConstrainedOptimizationPack/include/MeritFuncCalcNLE.h"
#include "ConstrainedOptimizationPack/include/MeritFuncNLESqrResid.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/max_near_feas_step.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"
#include "LinAlgPack/include/LinAlgOpPack.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

namespace ReducedSpaceSQPPack {

LineSearch2ndOrderCorrect_Step::LineSearch2ndOrderCorrect_Step(
		  const direct_ls_sqp_ptr_t&		direct_ls_sqp
		, const merit_func_ptr_t&			merit_func
		, const direct_ls_newton_ptr_t&		direct_ls_newton
		, value_type						eta
		, ENewtonOutputLevel				newton_olevel
		, value_type						constr_norm_threshold
		, int								after_k_iter
		, EForcedConstrReduction			forced_constr_reduction
		, value_type						max_step_ratio
		, int								max_newton_iter			)
	:
		  direct_ls_sqp_(direct_ls_sqp)
		, merit_func_(merit_func)
		, direct_ls_newton_(direct_ls_newton)
		, eta_(eta)
		, newton_olevel_(newton_olevel)
		, constr_norm_threshold_(constr_norm_threshold)
		, after_k_iter_(after_k_iter)
		, forced_constr_reduction_(forced_constr_reduction)
		, max_step_ratio_(max_step_ratio)
		, max_newton_iter_(max_newton_iter)
		, considering_correction_(false)
{}

bool LineSearch2ndOrderCorrect_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using std::setw;

	using LinAlgPack::dot;
	using LinAlgPack::norm_inf;
	using LinAlgPack::V_VpV;
	using LinAlgPack::V_VmV;
	using LinAlgPack::Vp_StV;
	using LinAlgPack::Vt_S;

	using LinAlgOpPack::Vp_V;
	using LinAlgOpPack::V_MtV;

	using SparseLinAlgPack::max_near_feas_step;

	using ConstrainedOptimizationPack::print_vector_change_stats;

	typedef LineSearch2ndOrderCorrect_Step	this_t;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();
	NLP			&nlp	= algo.nlp();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();
	out << std::boolalpha;

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// /////////////////////////////////////////
	// Set references to iteration quantities
	//
	// Set k+1 first then go back to get k to ensure
	// we have backward storage.
	
	Vector
		&x_kp1 = s.x().set_k(+1).v();
	value_type
		&f_kp1 = s.f().set_k(+1);
	Vector
		&c_kp1 = s.c().set_k(+1).v();

	const value_type
		&f_k = s.f().get_k(0);
	const Vector
		&c_k = s.c().get_k(0).v();
	const Vector
		&x_k = s.x().get_k(0).v();
	const Vector
		&d_k = s.d().get_k(0).v();
	value_type
		&alpha_k = s.alpha().get_k(0);

	// /////////////////////////////////////
	// Compute Dphi_k, phi_kp1 and phi_k

	// Dphi_k
	const value_type
		Dphi_k = merit_func().deriv();
	if( Dphi_k >= 0 ) {
		throw LineSearchFailure( "LineSearch2ndOrderCorrect_Step::do_step(...) : " 
			"Error, d_k is not a descent direction for the merit function " );
	}

	// ph_kp1
	value_type
		&phi_kp1 = s.phi().set_k(+1) = merit_func().value( f_kp1, c_kp1 );

	// Must compute phi(x) at the base point x_k since the penalty parameter may have changed.
	const value_type
		&phi_k = s.phi().set_k(0) = merit_func().value( f_k, c_k );

	// //////////////////////////////////////
	// Setup the calculation merit function

	// Here f_kp1, and c_kp1 are updated at the same time the
	// line search is being performed.
	nlp.set_f( &f_kp1 );
	nlp.set_c( &c_kp1 );
	MeritFuncCalcNLP
		phi_calc( &merit_func(), &nlp );

	// ///////////////////////////////////////////////
	// Concider 2nd order correction if near solution

	if( !considering_correction_ ) {
		const value_type
			nrm_c_x  = s.c().get_k(0).norm_inf();
		if( nrm_c_x <= constr_norm_threshold() && s.k() >= after_k_iter() ) {
			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
				out	<< "\nConsider the 2nd order correction for x_kp1 = x_k + d_k + w from now on:\n"
					<< "||c_k|| = " << nrm_c_x << " <= constr_norm_threshold = "
						<< constr_norm_threshold() << std::endl
					<< "k = " << s.k() << " <= after_k_iter = "
						<< after_k_iter() << std::endl;
			}
			considering_correction_ = true;
		}
	}

	// //////////////////////////////
	// See if we can take a full step

	bool chose_point = false;

	const value_type frac_phi = phi_k + eta() * Dphi_k;
	const bool armijo_test = phi_kp1 <= frac_phi;
	if( armijo_test ) {
		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			out	<< "\nAccepting full step x_kp1 = x_k + d_k\n";
		}
		chose_point = true;	// The point meets the Armijo test.
	}

	// This is storage for function and gradient evaluations for
	// the trial newton points and must be remembered for latter
	value_type f_xdww;
	Vector     c_xdww;
	Vector w,		// Full correction after completed computation.
		   xdww;	// Will be set to xdw + sum( w(newton_i), newton_i = 1... )
					//     where w(itr) is the local corrections for the current
					//		newton iteration.
	bool use_correction = false;

	bool considered_correction = ( considering_correction_ && !chose_point );
	if( considered_correction ) {

		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			out	<< "\nConsidering whether to compute a 2nd order correction for\n"
					<< "x_kp1 = x_k + alpha_k * d_k + alpha_k^2 * w ...\n";
		}

		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			const value_type obj_descent = dot( s.Gf().get_k(0)(), d_k() );
			out	<< "\nGf_k' * d_k = " << obj_descent << std::endl;
			if( obj_descent >= 0.0 ) {
				out	<< "\nWarning, this may not work well with Gf_k'*d_k >= 0.0\n";
			}
		}

		// Merit function for newton line searches
		ConstrainedOptimizationPack::MeritFuncNLESqrResid
			phi_c;

		Vector
			xdw = x_kp1;	// Will be set to x + d + sum(w(i),i=1..itr-1)
							//     where w(i) are previous local corrections
		value_type
			phi_c_x    = phi_c.value( c_k() ),
			phi_c_xd   = phi_c.value( c_kp1() ),
			phi_c_xdw  = phi_c_xd,		// No correction is computed yet so w = 0
			phi_c_xdww = phi_c_xdw,
			nrm_d	   = norm_inf( d_k() );

		// Merit function for newton line searches
		nlp.set_f( &(f_xdww = f_kp1) );
		nlp.set_c( &(c_xdww = c_kp1) );
		ConstrainedOptimizationPack::MeritFuncCalcNLE
			phi_c_calc( &phi_c, &nlp );

		Vector wy(s.decomp_sys().r());	// Range space wy (see latter).

		if( phi_c_xd < phi_c_x ) {
			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
				out	<< "\nphi_c(c(x_k+d_k)) = " << phi_c_xd
						<< " < phi_c(c(x_k)) = " << phi_c_x << std::endl
					<< "No need for a 2nd order correciton, perform regular line search ... \n";
			}
			use_correction = false;
		}
		else {
			// Try to compute a second order correction term.

			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
				out	<< "\nphi_c(c(x_k+d_k)) = " << phi_c_xdw
						<< " >= phi_c(c(x_k)) = " << phi_c_x << std::endl
					<< "Lets try to compute a second order correction w ... \n";
			}
			
			// Print header for summary information
			int owidth = 22;
			int prec = 8;
			out	<< std::setprecision(prec);
			if( newton_olevel() == PRINT_NEWTON_SUMMARY_INFO ) {
				out	<< "\nStarting Newton iterations\n\n"
					<< "\nphi_c_x   = "	<< phi_c_x 
					<< "\nphi_c_xd  = "	<< phi_c_xd
					<< "\n||d_k||nf = "	<< phi_c_xd << "\n\n"
					<< setw(5)			<< "it"
					<< setw(owidth)		<< "||w||inf"
					<< setw(owidth)		<< "u"
					<< setw(owidth)		<< "step_ratio"
					<< setw(5)			<< "lsit"
					<< setw(owidth)		<< "a"
					<< setw(owidth)		<< "phi_c_xdww"
					<< setw(owidth)		<< "phi_c_xdww-phi_c_x"
					<< setw(owidth)		<< "phi_c_xdww-phi_c_xd\n"
					<< setw(5)			<< "----"
					<< setw(owidth)		<< "-------------------"
					<< setw(owidth)		<< "-------------------"
					<< setw(owidth)		<< "-------------------"
					<< setw(5)			<< "----"
					<< setw(owidth)		<< "-------------------"
					<< setw(owidth)		<< "-------------------"
					<< setw(owidth)		<< "-------------------"
					<< setw(owidth)		<< "-------------------\n";
			}

			int newton_i;
			for( newton_i = 1; newton_i <= max_newton_iter(); ++newton_i ) {
				
				if( (int)newton_olevel() >= (int)this_t::PRINT_NEWTON_STEPS ) {
					out << "\n**** newton_i = " << newton_i << std::endl;
				}

				// ToDo: The calculation of the correction needs to be
				// delegated to somewhere else so that we can compute it
				// different ways.

				// Compute the local second order correction w which in the end
				// the full correction will be sum(w(i)).
				//
				// Compute w s.t. Gc'*w + c(xdw) = 0
				//
				// To find such a w:
				//
				// Gc'*w + c(xdw)
				//	=> Gc'* (Z*wz + Y*wy) + c(xdw) = 0
				//
				// Set wz = 0 then solve:
				//
				// wy = -inv(Gc'*Y) * c(xdw)
				// w = Y*wy

				// wy = -inv(Gc'*Y) * c(xdw)
				// Note: c(xdw) was already computed when phi_c_calc(xdw) was computed.
				s.decomp_sys().solve_transAtY( nlp.c(), BLAS_Cpp::no_trans, &wy() );
				Vt_S( &wy(), -1.0 );

				// w = Y*wy
				V_MtV( &w, s.Y().get_k(0), BLAS_Cpp::no_trans, wy() );

				// End code to delagete to some where else

				value_type
					nrm_w = norm_inf(w());				

				if( (int)newton_olevel() >= (int)this_t::PRINT_NEWTON_STEPS ) {
					out << "\n||w||inf = " << nrm_w << std::endl;
				}

				if( (int)newton_olevel() >= (int)this_t::PRINT_NEWTON_VECTORS ) {
					out << "\nw = " << w();
				}

				// ////////////////////////////////
				// Cutting back w

				value_type a = 1.0;	// This is the alpha for your linesearch

				// Cut back w to be in the relaxed bounds.
				std::pair<value_type,value_type>
					u_steps = max_near_feas_step( s.x().get_k(0)(), w()
						, algo.nlp().xl(), algo.nlp().xu()
						, algo.algo_cntr().max_var_bounds_viol() );
				const value_type u = u_steps.first;

				if( u < a ) {
					if( (int)newton_olevel() >= (int)this_t::PRINT_NEWTON_STEPS ) {
						out << "\nCutting back w = (a=u) * w to be within relaxed bounds:\n"
							<< "u = " << u << std::endl;
					}
					a = u;
				}

				// Cut back step so x+d+sum(w(i)) is not too far from x+d
				value_type
					step_ratio = nrm_w / nrm_d;
				if( (int)newton_olevel() >= (int)this_t::PRINT_NEWTON_STEPS ) {
					out << "\nstep_ratio = ||w||inf/||d||inf = " << step_ratio
							<< std::endl;
				}
				if( a * step_ratio > max_step_ratio() ) {
					const value_type aa = a*(max_step_ratio()/step_ratio);
					if( (int)newton_olevel() >= (int)this_t::PRINT_NEWTON_STEPS ) {
						out << "\na*step_ratio = " << a*step_ratio
								<< " > max_step_ratio = " << max_step_ratio() << std::endl
							<< "Cutting back a = a*max_step_ratio/step_ratio = "
								<< aa << std::endl;
					}
					a = aa;
				}

				// /////////////////////////////////////////////////
				// Perform a line search along xdww = xdw + a * w
				
				if( (int)newton_olevel() >= (int)this_t::PRINT_NEWTON_STEPS ) {
					out << "\nPerform linesearch along xdww = xdw + a*w\n"
						<< "starting from a = " << a << " ...\n";
				}

				xdww = xdw();									// xdww = xdw + a*w
				Vp_StV( &xdww(), a, w() );
				phi_c.calc_deriv(nlp.c());	// Set the directional derivative at c(xdw)
				phi_c_xdww = phi_c_calc( xdww() );	// phi_c_xdww = phi(xdww)
				const VectorSlice xdw_w[2] = { xdw(), w() };
				MeritFuncCalc1DQuadratic
					phi_c_calc_1d( phi_c_calc, 1 , xdw_w, &xdww() );
				const bool
					ls_okay = direct_ls_newton().do_line_search(phi_c_calc_1d,phi_c_xdw
						,&a,&phi_c_xdww
						, (int)newton_olevel() >= (int)this_t::PRINT_NEWTON_STEPS 
							? &out : 0												);
				// Note that the last value c(x) computed but the line search is for
				// xdw + a*w.

				// Print line for summary output
				if( newton_olevel() == PRINT_NEWTON_SUMMARY_INFO ) {
					out	<< setw(5)			<< newton_i
						<< setw(owidth)		<< nrm_w
						<< setw(owidth)		<< u
						<< setw(owidth)		<< step_ratio
						<< setw(5)			<< direct_ls_newton().num_iterations()
						<< setw(owidth)		<< a
						<< setw(owidth)		<< phi_c_xdww
						<< setw(owidth)		<< (phi_c_xdww-phi_c_x)
						<< setw(owidth)		<< (phi_c_xdww-phi_c_xd) << std::endl;
				}

				if(!ls_okay) {
					if( (int)newton_olevel() >= (int)this_t::PRINT_NEWTON_SUMMARY_INFO ) {
						out << "\nMaximum number of linesearch iterations has been exceeded\n"
							<< "so forget about computing a correction ...\n";
					}
					use_correction = false;
					break;
				}

				// See if this point is okay
				bool good_correction = false;
				switch( forced_constr_reduction() ) {
					case CONSTR_LESS_X_D: {
						good_correction = ( phi_c_xdww < phi_c_xd );
						if( good_correction
							&& (int)newton_olevel() >= (int)this_t::PRINT_NEWTON_SUMMARY_INFO )
						{
							out << "\nphi_c(c(x_k+d_k+w)) = " << phi_c_xdww
									<< " < phi_c(c(x_k+d_k)) = " << phi_c_xd << std::endl;
						}
						break;
					}
					case CONSTR_LESS_X: {
						good_correction = ( phi_c_xdww < phi_c_x );
						if( good_correction
							&& (int)newton_olevel() >= (int)this_t::PRINT_NEWTON_SUMMARY_INFO )
						{
							out << "\nphi_c(c(x_k+d_k+w)) = " << phi_c_xdww
									<< " < phi_c(c(x_k)) = " << phi_c_x << std::endl;
						}
						break;
					}
					default:
						assert(0);
				}

				if(good_correction) {
					if( (int)newton_olevel() >= (int)this_t::PRINT_NEWTON_SUMMARY_INFO ) {
						out << "\nAccept this point and compute our full correction w ... \n";
					}
					// Compute the full correction and do a curved linesearch
					// w = xdww - x_kp1
					V_VmV( &w(), xdww(), x_kp1() );
					if( (int)newton_olevel() >= (int)this_t::PRINT_NEWTON_STEPS ) {
						out << "\n||w||inf = " << (a * nrm_w) << std::endl;
					}
					if( (int)newton_olevel() >= (int)this_t::PRINT_NEWTON_VECTORS ) {
						out << "\nw = " << w();
					}
					use_correction = true;
					break;
				}

				// Else perform another newton iteration.
				xdw       = xdww;
				phi_c_xdw = phi_c_xdww;

			}	// end for
			if( !use_correction ) {
				if( forced_constr_reduction() == CONSTR_LESS_X_D ) {
					if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
						out	<< "\nDam! This is really bad!\n"
							<< "We where only looking for point phi_c(c(x_k+d_k+w)"
							<< " < phi_c(c(x_k+k_k) and we could not find it\n"
							<< " in the aloted number of newton iterations!\n"
							<< "Perhaps the Gc_k did not give us a descent direction?\n"
							<< "Just perform a standard line search from here ...\n";
					}
					use_correction = false;
				}
				else {
					if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
						out	<< "\nWe where looking for point phi_c(c(x_k+d_k+w))"
							<< " < phi_c(c(x_k)) and we could not find it.\n";
					}
					if( phi_c_xdww < phi_c_xd ) {
						if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
							out	<< "However, we did find a point less than phi_c(c(x_k+d_k))\n"
								<< "so lets use the correction anyway.\n";
						}
						// Compute the full correction and do a curved linesearch
						// w = xdww - x_kp1
						V_VmV( &w(), xdww(), x_kp1() );
						use_correction = true;
					}
					else {
						if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
							out	<< "Dam! We did not even find a point less than phi_c(c(x_k+d_k))\n"
								<< "just perform a standard line search along x_k + alpha_k*d_k.\n";
						}
						use_correction = false;
					}
				}
			}
		}	// end else from if phi_c_xdw > phi_c_x
	} // end considered_correction

	// //////////////////////////
	// Set up for the line search

	if( considered_correction ) {
		if( use_correction ) {
			// We are using the correction so setup the full step for the
			// NLP linesearch to come.
			Vp_V( &x_kp1(), w() );	// Set point to x_kp1 = x_k + d_k + w
			f_kp1 = nlp.f();		// Here f and c where computed at x_k+d_k+w
			c_kp1 = nlp.c()();
			phi_kp1 = merit_func().value( f_kp1, c_kp1 );
		}
		else {
			// Just pretend the second order correction never happened
			// and we don't need to do anything.
		}
		// Set back the base point
		nlp.set_f( &f_kp1 );
		nlp.set_c( &c_kp1 );
	}

	// //////////////////////
	// Do the line search

	if( !chose_point ) {
		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			if( use_correction ) {
				out	<< "\nPerform a curved linesearch along:\n"
					<< "x_kp1 = x_k + alpha_k * d_k + alpha_k^2 * w ...\n";
			}
			else {
				out	<< "\nPerform a standard linesearch along:\n"
					<< "x_kp1 = x_k + alpha_k * d_k ...\n";
			}
		}
		const VectorSlice xdw[3] = { x_k(), d_k(), w() };
		MeritFuncCalc1DQuadratic
			phi_calc_1d( phi_calc, (use_correction?2:1) , xdw, &x_kp1() );
		if( !direct_ls_sqp().do_line_search( phi_calc_1d, phi_k, &alpha_k, &phi_kp1
			, static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ?
				&out : static_cast<std::ostream*>(0)	)		)
		{
			// If the line search failed but the value of the merit function is less than
			// the base point value then just accept it and move on.  This many be a
			// bad thing to do.

			const value_type
				scaled_ared		= (s.phi().get_k(0) - s.phi().get_k(+1))/s.phi().get_k(0),
				keep_on_frac	= 1.0e-10;	// Make adjustable?
			bool keep_on = scaled_ared < keep_on_frac;

			if( (int)olevel >= (int)PRINT_BASIC_ALGORITHM_INFO )
			{
				out
					<< "\nThe maximum number of linesearch iterations has been exceeded "
					<< "(k = " << algo.state().k() << ")\n"
					<< "(phi_k - phi_kp1)/phi_k = " << scaled_ared;
//				if(keep_on) {
//					out
//						<< " < " << keep_on_frac
//						<< "\nso we will accept to step and move on.\n";
//				}
//				else {
//					out
//						<< " > " << keep_on_frac
//						<< "\nso we will reject the step and declare a line search failure.\n";
//				}
			}
//
//			if( keep_on ) return true;
			
			throw LineSearchFailure( "LineSearch2ndOrderCorrect_Step::do_step(): "
									 "Error, Line search failure" );
		}
	}

	if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
		out	<< "\nalpha_k      = "	<< alpha_k				<< std::endl
			<< "\n||x_kp1||inf = "	<< norm_inf( x_kp1 )	<< std::endl
			<< "\nf_kp1        = "	<< f_kp1				<< std::endl
			<< "\n||c_kp1||inf = "	<< norm_inf(c_kp1)		<< std::endl
			<< "\nphi_kp1      = "	<< phi_kp1				<< std::endl;
	}

	if( (int)olevel >= (int)PRINT_VECTORS ) {
		out << "\nx_kp1 =\n"	<< x_kp1
			<< "\nc_kp1 =\n"	<< c_kp1;
	}

	return true;
}

void LineSearch2ndOrderCorrect_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out	<< L << "*** Calculate a second order correction when near solution.\n"
		<< L << "*** If we can compute a correction w then perform a curved\n"
		<< L << "*** line search along x_kp1 = x_k + alpha_k * d_k + alpha_k^2 * w.\n"
		<< L << "default: eta                     = " << eta() << std::endl
		<< L << "         constr_norm_threshold   = " << constr_norm_threshold() << std::endl
		<< L << "         after_k_iter            = " << after_k_iter() << std::endl
		<< L << "         forced_constr_reduction = CONSTR_LESS_X_D\n"
		<< L << "         max_step_ratio          = " << max_step_ratio() << std::endl
		<< L << "         max_newton_iter         = " << max_newton_iter() << std::endl
		<< L << "         considering_correction  = false\n"
		<< L << "begin definition of NLP merit function phi.value(f(x),c(x)):\n";

	merit_func().print_merit_func( out, L + "    " );
	
	out	<< L << "end definition\n"
		<< L << "Dphi_k = phi.deriv()\n"
		<< L << "if Dphi_k >= 0 then\n"
		<< L << "    throw line_search_failure\n"
		<< L << "end\n"
		<< L << "phi_kp1 = phi_k.value(f_kp1,c_kp1)\n"
		<< L << "phi_k = phi.value(f_k,c_k)\n"
		<< L << "if considering_correction == false then\n"
		<< L << "    if (norm_inf_c_k < constr_norm_threshold) and (k >= after_k_iter) then\n"
		<< L << "        considering_correction = true\n"
		<< L << "    end\n"
		<< L << "end\n"
		<< L << "chose_point = false\n"
		<< L << "if phi_kp1 < phi_k + eta * Dphi_k then\n"
		<< L << "    chose_point = true\n"
		<< L << "else\n"
		<< L << "if (considering_correction == true) and (chose_point == false) then\n"
		<< L << "    considered_correction = true\n"
		<< L << "    begin definition of c(x) merit function phi_c.value(c(x)):\n";

	ConstrainedOptimizationPack::MeritFuncNLESqrResid().print_merit_func(
		out, L + "        " );
	
	out	<< L << "    end definition\n"
		<< L << "    xdw = x_kp1;\n"
		<< L << "    phi_c_x = phi_c.value(c_k);\n"
		<< L << "    phi_c_xd = phi_c.value(c_kp1);\n"
		<< L << "    phi_c_xdw = phi_c_xd;\n"
		<< L << "    phi_c_xdww = phi_c_xdw;\n"
		<< L << "    if phi_c_xd < phi_c_x then\n"
		<< L << "        *** There is no need to perform a correction.\n"
		<< L << "        use_correction = false;\n"
		<< L << "    else\n"
		<< L << "        *** Lets try to compute a correction by performing\n"
		<< L << "        *** a series of newton steps to compute local steps w\n"
		<< L << "        for newton_i = 1...max_newton_itr\n"
		<< L << "            wy = -inv(Gc'*Y)*c(xdw);\n"
		<< L << "            w = Y*wy;\n"
		<< L << "            Find the largest positive step u where the slightly\n"
		<< L << "            relaxed variable bounds:\n"
		<< L << "                xl - delta <= xdw + u * w <= xu + delta\n"
		<< L << "            are strictly satisfied (where delta = max_var_bounds_viol).\n"
		<< L << "            a = min(1,u);\n"
		<< L << "            step_ratio = norm(w,inf)/norm(d,inf);\n"
		<< L << "            a = min(a,max_step_ratio/step_ratio);\n"
		<< L << "            Perform line search for phi_c.value(c(xdww = xdw+a*w))\n"
		<< L << "            starting from a and compute:\n"
		<< L << "                a,\n"
		<< L << "                xdww = xdw + a * w,\n"
		<< L << "                phi_c_xdww = phi_c.value(c(xdww))\n"
		<< L << "            print summary information;\n"
		<< L << "            if line search failed then\n"
		<< L << "                use_correction = false;\n"
		<< L << "                exit for loop;\n"
		<< L << "            end\n"
		<< L << "            *** Determine if this is sufficent reduction in c(x) error\n"
		<< L << "            if forced_constr_reduction == CONSTR_LESS_X_D then\n"
		<< L << "                good_correction = (phi_c.value(c(xdww))\n"
		<< L << "                                        < phi_c.value(c(x_k+d_k)));\n"
		<< L << "            else if forced_constr_reduction == CONSTR_LESS_X then\n"
		<< L << "                good_correction = (phi_c.value(c(xdww))\n"
		<< L << "                                        < phi_c.value(c(x_k)));\n"
		<< L << "            end\n"
		<< L << "            if good_correction == true then\n"
		<< L << "                w = xdww - (x_k+d_k);\n"
		<< L << "                use_correction = true;\n"
		<< L << "                exit for loop;\n"
		<< L << "            end\n"
		<< L << "            *** This is not sufficient reduction is c(x) error yet.\n"
		<< L << "            xdw = xdww;\n"
		<< L << "            phi_c_xdw = phi_c_xdww;\n"
		<< L << "        end\n"
		<< L << "        if use_correction == false then\n"
		<< L << "            if forced_constr_reduction == CONSTR_LESS_X_D then\n"
		<< L << "               *** Dam! We could not find a point phi_c_xdww < phi_c_xd.\n"
		<< L << "               *** Perhaps Gc_k does not give a descent direction for phi_c\n"
		<< L << "            else if forced_constr_reduction == CONSTR_LESS_X then\n"
		<< L << "               if phi_c_dww < phi_c_xd then\n"
		<< L << "                   *** Accept this correction anyway.\n"
		<< L << "                   use_correction = true\n"
		<< L << "               else\n"
		<< L << "                   *** Dam! we could not find any reduction in phi_c so\n"
		<< L << "                   *** Perhaps Gc_k does not give a descent direction for phi_c\n"
		<< L << "            end\n"
		<< L << "        end\n"
		<< L << "    end\n"
		<< L << "end\n"
		<< L << "if chose_point == false then\n"
		<< L << "    if use_correction == true then\n"
		<< L << "        Perform line search along x_kp1 = x_k + alpha_k * d_k + alpha_k^2 * w\n"
		<< L << "    else\n"
		<< L << "        Perform line search along x_kp1 = x_k + alpha_k * d_k\n"
		<< L << "    end\n"
		<< L << "    begin direct line search : \"" << typeid(direct_ls_sqp()).name() << "\"\n";

	direct_ls_sqp().print_algorithm( out, L + "        " );

	out
		<< L << "    end direct line search\n"
		<< L << "    if maximum number of linesearch iterations are exceeded then\n"
		<< L << "        throw line_search_failure\n"
		<< L << "    end\n"
		<< L << "end\n";
}

}	// end namespace ReducedSpaceSQPPack