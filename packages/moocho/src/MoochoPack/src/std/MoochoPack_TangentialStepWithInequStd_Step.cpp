// ////////////////////////////////////////////////////////////////////////////
// IndepDirecWithBoundsStd_Step.cpp
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

#include <ostream>
#include <sstream>

#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackExceptions.h"
#include "ReducedSpaceSQPPack/include/std/IndepDirecWithBoundsStd_Step.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackExceptions.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/ComputeMinMult.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "ConstrainedOptimizationPack/include/ZVarReductMatrix.h"
#include "SparseLinAlgPack/include/MatrixWithOpFactorized.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/include/SpVectorOut.h"
#include "SparseLinAlgPack/include/max_near_feas_step.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"
#include "Misc/include/dynamic_cast_verbose.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

namespace ReducedSpaceSQPPack {

IndepDirecWithBoundsStd_Step::IndepDirecWithBoundsStd_Step(
	const qp_solver_ptr_t       &qp_solver
	,const qp_tester_ptr_t      &qp_tester
	,value_type                 warm_start_frac
	,EQPTesting                 qp_testing
	,bool                       primal_feasible_point_error
	,bool                       dual_feasible_point_error
	)
	:qp_solver_(qp_solver)
	,qp_tester_(qp_tester)
	,warm_start_frac_(warm_start_frac)
	,qp_testing_(qp_testing)
	,primal_feasible_point_error_(primal_feasible_point_error)
	,dual_feasible_point_error_(dual_feasible_point_error)
{}

bool IndepDirecWithBoundsStd_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using DynamicCastHelperPack::dyn_cast;
	using ::fabs;
	using LinAlgPack::norm_inf;
	using LinAlgPack::Vt_S;
	using LinAlgPack::V_VpV;
	using LinAlgPack::Vp_StV;
	using LinAlgOpPack::Vp_V;
	using LinAlgOpPack::V_StV;
	using LinAlgOpPack::V_MtV;
	using SparseLinAlgPack::norm_inf;
	using ConstrainedOptimizationPack::min_abs;
	using SparseLinAlgPack::max_near_feas_step;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();
	const bool check_results = algo.algo_cntr().check_results();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// problem dimensions
	const size_type
		n = algo.nlp().n(),
		m = algo.nlp().m(),
		r = algo.nlp().r();

	// Comupte qp_grad which is an approximation to rGf + Z' * HL * Y * py

	// qp_grad = rGf
	VectorSlice
		qp_grad = s.qp_grad().set_k(0).v() = s.rGf().get_k(0)();

	// qp_grad += zeta * w
	if(s.w().updated_k(0)) {
		if(s.zeta().updated_k(0))
			Vp_StV( &qp_grad, s.zeta().get_k(0), s.w().get_k(0)() );
		else
			Vp_V( &qp_grad, s.w().get_k(0)() );
	}

	// Set the bounds for:
	// 		dl <= Z*pz + Y*py <= du		->		dl - Ypy <= Z*pz <= du - Ypz

	const SparseBounds
		&d_bounds = d_bounds_(s).get_k(0);

	// bl = dl - Ypy_k

	const Vector
		&Ypy_k	= s.Ypy().get_k(0).v();
	const SpVectorSlice
		&dl = d_bounds.l,
		&du = d_bounds.u;

	SpVector bl;
	bl.uninitialized_resize( dl.size(), dl.nz(), dl.nz() );
	bl.assume_sorted(true);
	if(dl.nz()) {
		SpVectorSlice::const_iterator
			dl_itr		= dl.begin(),
			dl_itr_end	= dl.end();
		SpVector::iterator
			bl_itr		= bl.begin();
		for(; dl_itr != dl_itr_end; ++dl_itr, ++bl_itr) {
			const SpVectorSlice::element_type::indice_type i =  dl_itr->indice() + dl.offset();
			bl_itr->initialize( i, dl_itr->value() - Ypy_k(i) );
		}
	}

	// bu = du - Ypy_k

	SpVector bu;
	bu.uninitialized_resize( du.size(), du.nz(), du.nz() );
	bu.assume_sorted(true);
	if(du.nz()) {
		SpVectorSlice::const_iterator
			du_itr		= du.begin(),
			du_itr_end	= du.end();
		SpVector::iterator
			bu_itr		= bu.begin();
		for(; du_itr != du_itr_end; ++du_itr, ++bu_itr) {
			const SpVectorSlice::element_type::indice_type i =  du_itr->indice() + du.offset();
			bu_itr->initialize( i, du_itr->value() - Ypy_k(i) );
		}
	}

	// Print out the QP bounds for the constraints
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nqp_grad = \n" << qp_grad();
		out << "\nbl = \n" << bl();
		out << "\nbu = \n" << bu();
	}

	if(algo.algo_cntr().check_results()) {
		bl.assert_valid_and_sorted();
		bu.assert_valid_and_sorted();
	}

	//
	// Determine if we should perform a warm start or not.
	//
	bool do_warm_start = false;
	if( act_set_stats_(s).updated_k(-1) ) {
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out	<< "\nDetermining if the QP should use a warm start ...\n";
		}
		// We need to see if we should preform a warm start for the next iteration
		ActSetStats &stats = act_set_stats_(s).get_k(-1);
		const size_type
			num_active = stats.num_active(),
			num_adds   = stats.num_adds(),
			num_drops  = stats.num_drops();
		const value_type
			frac_same
			= ( num_adds == ActSetStats::NOT_KNOWN || num_active == 0
				? 0.0
				: std::_MAX(((double)(num_active)-num_adds-num_drops) / num_active, 0.0 ) );
		do_warm_start = ( num_active > 0 && frac_same >= warm_start_frac() );
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out << "\nnum_active = " << num_active;
			if( num_active ) {
				out	<< "\nmax(num_active-num_adds-num_drops,0)/(num_active) = "
					<< "max("<<num_active<<"-"<<num_adds<<"-"<<num_drops<<",0)/("<<num_active<<") = "
					<< frac_same;
				if( do_warm_start )
					out << " >= ";
				else
					out << " < ";
				out << "warm_start_frac = " << warm_start_frac();
			}
			if( do_warm_start )
				out << "\nUse a warm start this time!\n";
			else
				out << "\nDon't use a warm start this time!\n";
		}
	}

	// Use active set from last iteration as an estimate for current active set
	// if we are to use a warm start.
	// 
	// ToDo: If the selection of dependent and independent variables changes
	// , then you will have to adjust this or not perform a warm start at all!
	if( do_warm_start ) {
		const SpVector &nu_km1 = s.nu().get_k(-1);
		s.nu().set_k(0) = nu_km1;
	}
	else {
		s.nu().set_k(0).resize(n,n-r+1);	// no non-zero elements initially
	}

	// *** ToDo: Issolate out the following code into a function call.

	//
	// Setup the reduced QP subproblem
	//
	// The call to the QP is setup for the more flexible call to the QPSolverRelaxed
	// interface to deal with the two independent variabilities of using simple
	// bounds for pz or not and including extra equality constraints or not.
	// If this method of calling the QP solver were not used then 4 separate
	// calls to solve_qp(...) would have to be included to handle the four possible
	// QP formulations.
	//

	Vector &pz_k = s.pz().set_k(0).v();
	pz_k.resize(n-r);
	const MatrixWithOp
		&Z_k = s.Z().get_k(0);

	// The numeric arguments for the QP solver (in the nomenclatrue of QPSolverRelaxed)

	VectorSlice 			qp_g		= qp_grad();
	const MatrixWithOp& 	qp_G 		= s.rHL().get_k(0);
	const value_type		qp_etaL 	= 0.0;
	SpVectorSlice			qp_dL(NULL,0,0,n-r);	// If nz() == 0 then no simple bounds
	SpVectorSlice			qp_dU(NULL,0,0,n-r);
	const MatrixWithOp		*qp_E		= NULL;
	BLAS_Cpp::Transp		qp_trans_E	= BLAS_Cpp::no_trans;
	VectorSlice             qp_b;
	SpVectorSlice			qp_eL(NULL,0,0,n);
	SpVectorSlice			qp_eU(NULL,0,0,n);
	const MatrixWithOp		*qp_F		= NULL;
	BLAS_Cpp::Transp		qp_trans_F	= BLAS_Cpp::no_trans;
	VectorSlice				qp_f;
	value_type				qp_eta      = 0;
	VectorSlice				qp_d		= pz_k();
	SpVector				*qp_nu		= NULL;
	SpVector				*qp_mu		= NULL;
	VectorSlice				qp_Ed;
	VectorSlice				qp_lambda;

	SpVector _nu_pz, _nu_Dpz; 	// Possible storage for multiplers for separate inequality
								// constriants for pz.
	Vector _Dpz;				// Possible storage for D*pz computed by QP solver?

	//
	// Determine if we can use simple bounds on pz.
	// 
	// If we have a variable reduction null space matrix
	// (with any choice for Y) then:
	// 
	// d = Z*pz + (1-eta) * Y*py
	// 
	// [ d(dep)   ]  = [ D ] * pz  + (1-eta) * [ Ypy(dep)   ]
	// [ d(indep) ]    [ I ]                   [ Ypy(indep) ]
	// 
	// For a cooridinate decomposition (Y = [ I ; 0 ]) then Ypy(indep) = 0 and
	// in this case the bounds on d(indep) become simple bounds on pz even
	// with the relaxation.
	// 
	// Otherwise, we can not use simple variable bounds and implement the
	// relaxation properly.
	// 

	const ZVarReductMatrix
		*Zvr = dynamic_cast<const ZVarReductMatrix*>( &Z_k );
	Range1D
		indep 	= Zvr ? Zvr->indep() : Range1D(),
		dep		= Zvr ? Zvr->dep()   : Range1D();

	const bool
		use_simple_pz_bounds = ( Zvr!=NULL && norm_inf(s.Ypy().get_k(0).v()(indep))==0.0 );

	if( use_simple_pz_bounds ) {

		// Set simple bound constraints on pz
		qp_dL.bind( bl(indep) );
		qp_dU.bind( bu(indep) );
		qp_nu = &( _nu_pz = s.nu().get_k(0)(indep) );	// warm start?
		
		// Set general inequality constraints for D*pz
		qp_E = &Zvr->D();
		qp_b.bind( Ypy_k(dep) );
		qp_eL.bind( bl(dep) );
		qp_eU.bind( bu(dep) );
		qp_mu = &( _nu_Dpz = s.nu().get_k(0)(dep) );	// warm start?
		_Dpz.resize(r);
		qp_Ed.bind(_Dpz());	// Part of Zpz_k will be updated directly!

	}
	else {

		// There are no simple bounds!
		qp_nu = &_nu_pz; // Should be left with zero nonzeros

		// Set general inequality constraints for Z*pz
		qp_E = &Z_k;
		qp_b.bind( Ypy_k() );
		qp_eL.bind( bl() );
		qp_eU.bind( bu() );
		qp_mu = &s.nu().set_k(0);	// warm start??
		Vector &Zpz_k = s.Zpz().set_k(0).v();
		Zpz_k.resize(n);
		qp_Ed.bind(Zpz_k());	// Zpz_k will be updated directly!
	}

	// Set the general equality constriants (if they exist)
	Vector q(m-r);
	Range1D undecomp = s.con_undecomp();
	if( m > r ) {
		// q = U_k * py_k + c_k(undecomp)
		V_MtV( &q, s.U().get_k(0), BLAS_Cpp::no_trans, s.py().get_k(0)() );
		Vp_V( &q(), s.c().get_k(0).v()(undecomp) );
		// Must resize for the undecomposed constriants if it has not already been
		if( !s.lambda().updated_k(0) )
			s.lambda().set_k(0).v().resize(n,0.0);
		qp_F = &s.V().get_k(0);
		qp_f.bind( q() );
		qp_lambda.bind( s.lambda().set_k(0).v()(undecomp) );
	}

	// Setup the rest of the arguments

	QPSolverRelaxed::EOutputLevel
		qp_olevel;
	switch( olevel ) {
		case PRINT_NOTHING:
			qp_olevel = QPSolverRelaxed::PRINT_NONE;
			break;
		case PRINT_BASIC_ALGORITHM_INFO:
			qp_olevel = QPSolverRelaxed::PRINT_NONE;
			break;
		case PRINT_ALGORITHM_STEPS:
			qp_olevel = QPSolverRelaxed::PRINT_BASIC_INFO;
			break;
		case PRINT_ACTIVE_SET:
			qp_olevel = QPSolverRelaxed::PRINT_ITER_SUMMARY;
			break;
		case PRINT_VECTORS:
			qp_olevel = QPSolverRelaxed::PRINT_ITER_VECTORS;
			break;
		case PRINT_ITERATION_QUANTITIES:
			qp_olevel = QPSolverRelaxed::PRINT_EVERY_THING;
			break;
		default:
			assert(0);
	}

	//
	// Solve the QP
	// 
	const QPSolverStats::ESolutionType
		solution_type =
		qp_solver().solve_qp(
			 int(olevel) == int(PRINT_NOTHING) ? NULL : &out
			, qp_olevel
			, algo.algo_cntr().check_results()
				? QPSolverRelaxed::RUN_TESTS :  QPSolverRelaxed::NO_TESTS
			, qp_g, qp_G, qp_etaL, qp_dL, qp_dU
			, qp_E, qp_trans_E, qp_E ? &qp_b : NULL
				, qp_E ? &qp_eL : NULL, qp_E ? &qp_eU : NULL
			, qp_F, qp_trans_F, qp_F ? &qp_f : NULL
			, NULL
			, &qp_eta, &qp_d
			, qp_nu
			, qp_mu, qp_E ? &qp_Ed : NULL
			, qp_F ? &qp_lambda : NULL, NULL
			);

	//
	// Check the optimality conditions for the QP
	//
	std::ostringstream omsg;
	bool throw_qp_failure = false;
	if(		qp_testing() == QP_TEST
		|| ( qp_testing() == QP_TEST_DEFAULT && algo.algo_cntr().check_results() )  )
	{
		if( int(olevel) >= int(PRINT_ALGORITHM_STEPS) ) {
			out	<< "\nChecking the optimality conditions of the reduced QP subproblem ...\n";
		}
		if(!qp_tester().check_optimality_conditions(
			  solution_type
			, int(olevel) == int(PRINT_NOTHING) ? NULL : &out
			, int(olevel) >= int(PRINT_VECTORS) ? true : false
			, int(olevel) >= int(PRINT_ITERATION_QUANTITIES) ? true : false
			, qp_g, qp_G, qp_etaL, qp_dL, qp_dU
			, qp_E, qp_trans_E, qp_E ? &qp_b : NULL
				, qp_E ? &qp_eL : NULL, qp_E ? &qp_eU : NULL
			, qp_F, qp_trans_F, qp_F ? &qp_f : NULL
			, NULL
			, &qp_eta, &qp_d
			, qp_nu
			, qp_mu, qp_E ? &qp_Ed : NULL
			, qp_F ? &qp_lambda : NULL, NULL
			))
		{
			omsg << "\n*** Alert! at least one of the QP optimality conditions did not check out.\n";
			if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
				out << omsg.str();
			}
			throw_qp_failure = true;
		}
	}

	//
	// Set the solution
	//
	if( use_simple_pz_bounds ) {
		// Set nu_k from nu_indep and nu_dep
		SpVector &nu_k = s.nu().set_k(0);
		nu_k.resize(n,n-r+1); // Room for all active + 1 (relaxation free)
		typedef SpVector::element_type ele_t;
		if(_nu_Dpz.nz()) {
			// nu_dep first from _nu_Dpz
			const SpVector::difference_type off = _nu_Dpz.offset();
			for( SpVectorSlice::const_iterator itr = _nu_Dpz.begin(); itr != _nu_Dpz.end(); ++itr ) {
				nu_k.add_element( ele_t( itr->indice()+off,itr->value() ) );
			}
		}
		if(_nu_pz.nz() ) {
			// nu_indep last from _nu_pz (offset by r dependent variables)
			const SpVector::difference_type off = _nu_pz.offset();
			for( SpVectorSlice::const_iterator itr = _nu_pz.begin(); itr != _nu_pz.end(); ++itr ) {
				nu_k.add_element( ele_t( itr->indice() + off + r, itr->value() ) );
			}
		}
		nu_k.assume_sorted(true);
		if(check_results)
			nu_k.assert_valid_and_sorted();

		// Set Zpz_k
		Vector &Zpz_k = s.Zpz().set_k(0).v();
		Zpz_k.resize(n);
		Zpz_k(dep) = _Dpz;
		Zpz_k(indep) = pz_k;
	}
	else {
		// Everything should already be set!
	}

	// Set the solution statistics
	qp_solver_stats_(s).set_k(0) = qp_solver().get_qp_stats();

	// Cut back Ypy_k = (1-eta) * Ypy_k
	const value_type eps = std::numeric_limits<value_type>::epsilon();
	if( fabs(qp_eta - 0.0) > eps ) {
		if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out
				<< "\n*** Alert! the QP was infeasible (eta = "<<qp_eta<<").  Cutting back Ypy_k = (1.0 - eta)*Ypy_k  ...\n";
		}
		Vt_S( &s.Ypy().get_k(0).v()() , 1.0 - qp_eta );
	}

	// eta_k
	s.eta().set_k(0) = qp_eta;

	// *** End code to isolate out to be reused

	//
	// Modify the solution if we have to!
	// 
	switch(solution_type) {
		case QPSolverStats::OPTIMAL_SOLUTION:
			break;	// we are good!
		case QPSolverStats::PRIMAL_FEASIBLE_POINT:
		{
			omsg
				<< "\n*** Alert! the returned QP solution is PRIMAL_FEASIBLE_POINT but not optimal!\n";
			if( primal_feasible_point_error() )
				omsg
					<< "\n*** primal_feasible_point_error == true, this is an error!\n";
			if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
				out << omsg.str();
			}
			throw_qp_failure = primal_feasible_point_error();
			break;
		}	
		case QPSolverStats::DUAL_FEASIBLE_POINT:
		{
			omsg
				<< "\n*** Alert! the returned QP solution is DUAL_FEASIBLE_POINT"
				<< "\n*** but not optimal so we cut back the step ...\n";
			if( dual_feasible_point_error() )
				omsg
					<< "\n*** dual_feasible_point_error == true, this is an error!\n";
			if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
				out << omsg.str();
			}
			// Cut back the step to fit in the bounds
			// 
			// dl <= u*(Ypy_k+Zpz_k) <= du
			//
			Vector zero(n), d_tmp;
			zero = 0.0;
			V_VpV( &d_tmp, s.Ypy().get_k(0)(), s.Zpz().get_k(0)() );
			const std::pair<value_type,value_type>
				u_steps = max_near_feas_step( zero(), d_tmp(), dl, du, 0.0 );
			const value_type
				u = std::_MIN( u_steps.first, 1.0 ); // largest positive step size
			s.alpha().set_k(0) = u;
			if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
				out	<< "\nFinding u s.t. dl <= u*(Ypy_k+Zpz_k) <= du\n"
					<< "max step length u = " << u << std::endl
					<< "alpha_k = u = " << s.alpha().get_k(0) << std::endl;
			}
			throw_qp_failure = dual_feasible_point_error();
			break;
		}	
		case QPSolverStats::SUBOPTIMAL_POINT:
		{
			omsg
				<< "\n*** Alert!, the returned QP solution is SUBOPTIMAL_POINT!\n";
			if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
				out << omsg.str();
			}
			throw_qp_failure = true;
			break;
		}
		default:
			assert(0);	// should not happen!
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\n||pz||inf    = " << s.pz().get_k(0).norm_inf()
			<< "\nnu.nz()      = " << s.nu().get_k(0).nz()
			<< "\nmax(|nu(i)|) = " << norm_inf( s.nu().get_k(0)() )
			<< "\nmin(|nu(i)|) = " << min_abs( s.nu().get_k(0)() )
			;
		if( m > r ) out << "\n||lambda_k(undecomp)||inf = " << s.lambda().get_k(0).norm_inf();
		out	<< "\n||Zpz||2     = " << s.Zpz().get_k(0).norm_2()
			;
		if(qp_eta > 0.0) out << "\n||Ypy||2 = " << s.Ypy().get_k(0).norm_2();
		out << std::endl;
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\npz_k = \n" << s.pz().get_k(0)();
		out << "\nnu_k = \n" << s.nu().get_k(0)();
		if( m > r ) out << "\nlambda_k(undecomp) = \n" << s.lambda().get_k(0).v()(undecomp);
		out << "\nZpz_k = \n" << s.Zpz().get_k(0)();
		if(qp_eta > 0.0) out << "\nYpy = \n" << s.Ypy().get_k(0)();
	}

	if( qp_eta == 1.0 ) {
		omsg
			<< "IndepDirecWithBoundsStd_Step::do_step(...) : Error, a QP relaxation parameter\n"
			<< "of eta = " << qp_eta << " was calculated and therefore it must be assumed\n"
			<< "that the NLP's constraints are infeasible\n"
			<< "Throwing an InfeasibleConstraints exception!\n";
		if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out << omsg.str();
		}
		throw InfeasibleConstraints(omsg.str());
	}

	if( throw_qp_failure )
		throw QPFailure( omsg.str(), qp_solver().get_qp_stats() );

	return true;
}

void IndepDirecWithBoundsStd_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Calculate the null space step (independent direction) by solving a constrainted QP\n"
		<< L << "qp_grad_k = rGf_k\n"
		<< L << "if w_k is updated then set qp_grad_k = qp_grad_k + zeta_k * w_k\n"
		<< L << "bl = d_bounds_k.l - Ypy_k\n"
		<< L << "bu = d_bounds_k.u - Ypy_k\n"
		<< L << "if check_results == true then\n"
		<< L << "    assert that bl and bu are valid and sorted\n"
		<< L << "end\n"
		<< L << "etaL = 0.0\n"
		<< L << "*** Determine if we can use simple bounds on pz or not\n"
		<< L << "if Z_k is a variable reduction null space matrix and norm(Ypy_k(indep),0) == 0 then\n"
		<< L << "    use_simple_pz_bounds = true\n"
		<< L << "else\n"
		<< L << "    use_simple_pz_bounds = false\n"
		<< L << "end\n"
		<< L << "*** Setup QP arguments\n"
		<< L << "qp_g = qp_grad_k\n"
		<< L << "qp_G = rHL_k\n"
		<< L << "if use_simple_pz_bounds == true then\n"
		<< L << "    qp_dL = bl(indep),  qp_dU = bu(indep)\n"
		<< L << "    qp_E  = Z_k.D,      qp_b  = Ypy_k(dep)\n"
		<< L << "    qp_eL = bl(dep),    qp_eU = bu(dep)\n"
		<< L << "else\n"
		<< L << "    qp_dL = -inf,       qp_dU = +inf\n"
		<< L << "    qp_E  = Z_k,        qp_b  = Ypy_k\n"
		<< L << "    qp_eL = bl,         qp_eU = bu\n"
		<< L << "end\n"
		<< L << "if m > r then\n"
		<< L << "    qp_F  = V_k,        qp_f  = c_k(undecomp) + U_k * py_k\n"
		<< L << "else\n"
		<< L << "    qp_F  = empty,      qp_f  = empty\n"
		<< L << "end\n"
		<< L << "Given active set statistics (act_set_stats_km1)\n"
		<< L << "    frac_same = max(num_active-num_adds-num_drops,0)/(num_active)\n"
		<< L << "Use a warm start when frac_same >= warm_start_frac\n"
		<< L << "Solve the following QP to compute qp_d, qp_eta, qp_Ed = qp_E * qp_d\n"
		<< L << ",qp_nu, qp_mu and qp_lambda (" << typeid(qp_solver()).name() << "):\n"
		<< L << "    min      qp_g' * qp_d + 1/2 * qp_d' * qp_G * qp_d + M(eta)\n"
		<< L << "    qp_d <: R^(n-r)\n"
		<< L << "    s.t.\n"
		<< L << "             etaL  <=  eta\n"
		<< L << "             qp_dL <= qp_d <= qp_dU                          [qp_nu]\n"
		<< L << "             qp_eL <= qp_E * qp_d + (1-eta)*qp_b  <= qp_eU   [qp_mu]\n"
		<< L << "             qp_F * d_qp + (1-eta) * qp_f = 0                [qp_lambda]\n"
		<< L << "if (qp_testing==QP_TEST) or (fd_deriv_testing==QP_TEST_DEFAULT\n"
		<< L << "and check_results==true) then\n"
		<< L << "    Check the optimality conditions of the above QP\n"
		<< L << "    if the optimality conditions do not check out then\n"
		<< L << "        set throw_qp_failure = true\n"
		<< L << "    end\n"
		<< L << "end\n"
		<< L << "*** Set the solution to the QP subproblem\n"
		<< L << "pz_k  = qp_d\n"
		<< L << "eta_k = qp_eta\n"
		<< L << "if use_simple_pz_bounds == true then\n"
		<< L << "    nu_k(dep)    = qp_mu,  nu_k(indep)  = qp_nu\n"
		<< L << "    Zpz_k(dep)   = qp_Ed,  Zpz_k(indep) = pz_k\n"
		<< L << "else\n"
		<< L << "    nu_k  = qp_mu\n"
		<< L << "    Zpz_k = qp_Ed\n"
		<< L << "end\n"
		<< L << "if m > r then\n"
		<< L << "    lambda_k(undecomp) = qp_lambda\n"
		<< L << "end\n"
		<< L << "if eta_k > 0 then set Ypy_k = (1-eta_k) * Ypy_k\n"
		<< L << "if QP solution is suboptimal then\n"
		<< L << "    throw_qp_failure = true\n"
		<< L << "elseif QP solution is primal feasible (not optimal) then\n"
		<< L << "    throw_qp_failure = primal_feasible_point_error\n"
		<< L << "elseif QP solution is dual feasible (not optimal) then\n"
		<< L << "    find max u s.t.\n"
		<< L << "        d_bounds_k.l <= u*(Ypy_k+Zpz_k) <= d_bounds_k.u\n"
		<< L << "    alpha_k = u\n"
		<< L << "    throw_qp_failure = dual_feasible_point_error\n"
		<< L << "end\n"
		<< L << "if eta_k == 1.0 then\n"
		<< L << "   The constraints are infeasible!\n"
		<< L << "   throw InfeasibleConstraints(...)\n"
		<< L << "end\n"
		<< L << "if throw_qp_failure == true then\n"
		<< L << "    throw QPFailure(...)\n"
		<< L << "end\n"
		;
}

}	// end namespace ReducedSpaceSQPPack
