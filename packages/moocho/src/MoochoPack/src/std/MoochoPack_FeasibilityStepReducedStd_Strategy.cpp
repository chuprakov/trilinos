// ////////////////////////////////////////////////////////////////////////////////
// FeasibilityStepReducedStd_Strategy.cpp

#include "ReducedSpaceSQPPack/include/std/FeasibilityStepReducedStd_Strategy.h"
#include "ReducedSpaceSQPPack/include/rSQPAlgo.h"
#include "ReducedSpaceSQPPack/include/rSQPState.h"
#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackExceptions.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "ConstrainedOptimizationPack/include/ZVarReductMatrix.h"
#include "ConstrainedOptimizationPack/include/QPSolverStats.h"
#include "ConstrainedOptimizationPack/include/MatrixSymIdentity.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "SparseLinAlgPack/include/MatrixFactorized.h"
#include "SparseLinAlgPack/include/sparse_bounds.h"
#include "LinAlgPack/include/VectorOut.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "Misc/include/WorkspacePack.h"
#include "Misc/include/dynamic_cast_verbose.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

namespace ReducedSpaceSQPPack {

FeasibilityStepReducedStd_Strategy::FeasibilityStepReducedStd_Strategy(
	const quasi_range_space_step_ptr_t   &quasi_range_space_step
	,const qp_solver_ptr_t               &qp_solver
	,const qp_tester_ptr_t               &qp_tester
	,EQPObjective                        qp_objective
	,EQPTesting                          qp_testing
	)
	:quasi_range_space_step_(quasi_range_space_step)
	,qp_solver_(qp_solver)
	,qp_tester_(qp_tester)
	,qp_objective_(qp_objective)
	,qp_testing_(qp_testing)
	,current_k_(-1)
{}

bool FeasibilityStepReducedStd_Strategy::compute_feasibility_step(
	std::ostream& out, EJournalOutputLevel olevel, rSQPAlgo *algo, rSQPState *s
	,const VectorSlice& xo, const VectorSlice& c_xo, VectorSlice* w
  	)
{
	using SparseLinAlgPack::sparse_bounds_itr             sparse_bounds_itr;
	using ConstrainedOptimizationPack::MatrixSymIdentity  MatrixSymIdentity;
	using ConstrainedOptimizationPack::QPSolverStats      QPSolverStats;
	using DynamicCastHelperPack::dyn_cast;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	// problem dimensions
	const size_type
		n = algo->nlp().n(),
		m = algo->nlp().m(),
		r = algo->nlp().r();

	// Compute the quasi-range space step Ywy
	wsp::Workspace<value_type> Ywy_ws(wss,xo.size());
	VectorSlice                Ywy(&Ywy_ws[0],Ywy_ws.size());
	if(!quasi_range_space_step().solve_quasi_range_space_step(
		out,olevel,algo,s,xo,c_xo,&Ywy ))
		return false;

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\n||Ywy||2     = " << LinAlgPack::norm_2(Ywy);
		out << std::endl;
	}
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nYwy = \n" << Ywy;
	}

	//
	// Set the bounds on the null space QP subproblem:
	//
	// d_bounds_k.l <= (xo - x_k) + (1-eta) * Ywy + Z*wz <= d_bounds_k.u
	// =>
	// bl <= Z*wz - eta * Ywy <= bu
	//
	// where:   bl = d_bounds_k.l - (xo - x_k) - Ywy
	//          bu = d_bounds_k.u - (xo - x_k) - Ywy
	//
	// Above, fix the variables that are at an active bound as equalities
	// to maintain the same active-set.
	//
	const SparseBounds
		&d_bounds = d_bounds_(*s).get_k(0);
	const SpVectorSlice
		&dl = d_bounds.l,
		&du = d_bounds.u;
	const Vector
		&x_k = s->x().get_k(0).v();
	const SpVector
		&nu_k = s->nu().get_k(0);
	assert(nu_k.is_sorted());
	SpVector bl(n,n), bu(n,n);
	sparse_bounds_itr
		d_bnds_itr(dl.begin(),dl.end(),dl.offset(),du.begin(),du.end(),du.offset());
	SpVector::const_iterator
		nu_itr     = nu_k.begin(),
		nu_end     = nu_k.end();
	for( ; !d_bnds_itr.at_end(); ++d_bnds_itr ) {
		typedef SpVectorSlice::element_type ele_t;
		const size_type i = d_bnds_itr.indice();
		while( nu_itr != nu_end && nu_itr->indice() + nu_k.offset() < i )
			++nu_itr;
		if( nu_itr != nu_end && nu_itr->indice() + nu_k.offset() == i ) {
			const value_type
				act_bnd = nu_itr->value() > 0.0 ? d_bnds_itr.ubound() : d_bnds_itr.lbound();
			bl.add_element(ele_t( i, act_bnd - xo(i) + x_k(i) - Ywy(i) ));
			bu.add_element(ele_t( i, act_bnd - xo(i) + x_k(i) - Ywy(i) ));
		}
		else {
			if( d_bnds_itr.lbound() != -d_bnds_itr.big_bnd() )
				bl.add_element(ele_t(i,d_bnds_itr.lbound()  - xo(i) + x_k(i) - Ywy(i) ));
			if( d_bnds_itr.ubound() != +d_bnds_itr.big_bnd() )
				bu.add_element(ele_t(i, d_bnds_itr.ubound() - xo(i) + x_k(i) - Ywy(i) ));
		}
	}
	bl.assume_sorted(true);
	bu.assume_sorted(true);
	//
	// Setup the objective function for the null space QP subproblem
	//
	// 
	// OBJ_MIN_FULL_STEP
	//    min 1/2 * (Y*wy + Z*wz)'*(Y*wy + Z*wz) = 1/2*wy'*Y'*Y*wy + (Z'*Y*wy)'*wz + 1/2*wz'*(Z'*Z)*wz
	//    => grad = (Z'*Y*wy), Hess = Z'*Z
	//
	// OBJ_MIN_WZ
	//    min 1/2 * wz'*wz => grad = 0, Hess = I
	//
	// OBJ_RSQP
	//    min qp_grad_k'*wz + 1/2 * wz'*rHL_k*wz
	//    => grad = qp_grad, Hess = rHL_k
	//
	const MatrixWithOp
		&Z_k = s->Z().get_k(0);
	if( current_k_ != s->k() ) {
		if( qp_objective() != OBJ_RSQP )
			grad_store_.resize(n-r);
		if( qp_objective() == OBJ_MIN_FULL_STEP )
			Hess_store_.resize(n-r,n-r);
	}
	VectorSlice grad;
	switch(qp_objective())
	{
	    case OBJ_MIN_FULL_STEP: // grad = (Z'*Ywy), Hess = Z'*Z
		{
			grad.bind( grad_store_() );
			LinAlgOpPack::V_MtV( &grad, Z_k, BLAS_Cpp::trans, Ywy );
			assert(0); // ToDo: Need to implement initialization of Hess_store = Z'*Z (add to MatrixWithOp?)
			break;
		}
	    case OBJ_MIN_NULL_SPACE_STEP: // grad = 0, Hess = I
		{
			grad.bind( grad_store_() );
			MatrixSymIdentity
				*H_ptr = NULL;
			if( Hess_ptr_.get() == NULL || dynamic_cast<const MatrixSymIdentity*>(Hess_ptr_.get()) == NULL )
				Hess_ptr_ = new MatrixSymIdentity;
			if( current_k_ != s->k() ) {
				H_ptr = const_cast<MatrixSymIdentity*>(dynamic_cast<const MatrixSymIdentity*>(Hess_ptr_.get()));
				assert(H_ptr); // Should not be null!
				H_ptr->init_setup(n-r,1.0);
				grad = 0.0;
			}
			break;
		}
	    case OBJ_RSQP: // grad = qp_grad, Hess = rHL_k
		{
			grad.bind( s->qp_grad().get_k(0)() );
			Hess_ptr_ = Hess_ptr_t( &s->rHL().get_k(0), false );
			break;
		}
	    defaut:
			assert(0); // Not a valid option
	}

	//
	// Solve the null space subproblem
	//

	wsp::Workspace<value_type>  wz_ws(wss,n-r),Zwz_ws(wss,n);
	VectorSlice                 wz(&wz_ws[0],wz_ws.size());
	VectorSlice                 Zwz(&Zwz_ws[0],Zwz_ws.size());
	value_type                  qp_eta      = 0;

	bool throw_qp_failure = false;

	if( bl.nz() == 0 && bu.nz() == 0 && m-r == 0 ) {
		//
		// Just solve the unconstrainted QP
		//
		// wz = - inv(Hess)*grad
#ifdef _WINDOWS
		const MatrixFactorized &Hess = dynamic_cast<const MatrixFactorized&>(*Hess_ptr_);
#else
		const MatrixFactorized &Hess = dyn_cast<const MatrixFactorized>(*Hess_ptr_);
#endif
		SparseLinAlgPack::V_InvMtV( &wz, Hess, BLAS_Cpp::no_trans, grad );
		LinAlgPack::Vt_S(&wz,-1.0);
		// Zwz = Z*wz
		LinAlgOpPack::V_MtV( &Zwz, Z_k, BLAS_Cpp::no_trans, wz );
	}
	else {

		//
		// Set the arguments to the QP subproblem
		//

		VectorSlice 			qp_g		= grad;
		const MatrixWithOp& 	qp_G 		= *Hess_ptr_;
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
		VectorSlice				qp_d		= wz;
		SpVector				*qp_nu		= NULL;
		SpVector				*qp_mu		= NULL;
		VectorSlice				qp_Ed;
		VectorSlice				qp_lambda;

		SpVector _nu_wz, _nu_Dwz, 	// Possible storage for multiplers for separate inequality
			_nu;				// constriants for wz.
		Vector _Dwz;				// Possible storage for D*wz computed by QP solver?

		//
		// Determine if we can use simple bounds on wz.
		// 
		// If we have a variable reduction null space matrix
		// (with any choice for Y) then:
		// 
		// w = Z*wz + (1-eta) * Y*wy
		// 
		// [ w(dep)   ]  = [ D ] * wz  + (1-eta) * [ Ywy(dep)   ]
		// [ w(indep) ]    [ I ]                   [ Ywy(indep) ]
		// 
		// For a cooridinate decomposition (Y = [ I ; 0 ]) then Ywy(indep) = 0 and
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
			use_simple_pz_bounds = ( Zvr!=NULL && norm_inf(Ywy(indep))==0.0 );

		if( use_simple_pz_bounds ) {

			// Set simple bound constraints on pz
			qp_dL.bind( bl(indep) );
			qp_dU.bind( bu(indep) );
			qp_nu = &( _nu_wz = s->nu().get_k(0)(indep) );	// warm start?
		
			// Set general inequality constraints for D*pz
			qp_E = &Zvr->D();
			qp_b.bind( Ywy(dep) );
			qp_eL.bind( bl(dep) );
			qp_eU.bind( bu(dep) );
			qp_mu = &( _nu_Dwz = s->nu().get_k(0)(dep) );	// warm start?
			_Dwz.resize(r);
			qp_Ed.bind(_Dwz());	// Part of Zwz will be updated directly!

		}
		else {

			// Set general inequality constraints for Z*pz
			qp_E = &Z_k;
			qp_b.bind( Ywy() );
			qp_eL.bind( bl() );
			qp_eU.bind( bu() );
			qp_mu = &(_nu = s->nu().get_k(0));	// warm start??
			qp_Ed.bind(Zwz);	// Zwz
		}

		// Set the general equality constriants (if they exist)
		Vector q(m-r);
		Range1D undecomp = s->con_undecomp();
		if( m > r ) {
			assert(0); // ToDo: Implement when needed!
		}

		// Setup the rest of the arguments

		QPSolverRelaxed::EOutputLevel
			qp_olevel;
		switch( olevel ) {
		    case PRINT_NOTHING:
				qp_olevel = QPSolverRelaxed::PRINT_NONE;
				break;
		    case PRINT_BASIC_ALGORITHM_INFO:
				qp_olevel = QPSolverRelaxed::PRINT_BASIC_INFO;
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
				, algo->algo_cntr().check_results()
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
		if(		qp_testing() == QP_TEST
				|| ( qp_testing() == QP_TEST_DEFAULT && algo->algo_cntr().check_results() )  )
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

		if( solution_type !=  QPSolverStats::OPTIMAL_SOLUTION ) {
			if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
				out << "\nCould not solve the QP!\n";
			}
			return false;
		}

		//
		// Set the solution
		//
		if( use_simple_pz_bounds ) {
			// Set Zwz
			Zwz(dep)   = _Dwz;
			Zwz(indep) = wz;
		}
		else {
			// Everything should already be set!
		}

		// Cut back Ywy = (1-eta) * Ywy
		const value_type eps = std::numeric_limits<value_type>::epsilon();
		if( fabs(qp_eta - 0.0) > eps ) {
			if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
				out
					<< "\n*** Alert! the QP was infeasible (eta = "<<qp_eta<<").  Cutting back Ywy_k = (1.0 - eta)*Ywy  ...\n";
			}
			Vt_S( &Ywy , 1.0 - qp_eta );
		}
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\n||wz||inf    = " << LinAlgPack::norm_inf(wz);
		out	<< "\n||Zwz||2     = " << LinAlgPack::norm_2(Zwz);
		if(qp_eta > 0.0) out << "\n||Ypy||2 = " << LinAlgPack::norm_2(Ywy);
		out << std::endl;
	}
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nwz = \n" << wz;
		out << "\nZwz = \n" << Zwz;
		if(qp_eta > 0.0) out << "\nYwy = \n" << Ywy;
	}
	if( qp_eta == 1.0 ) {
		std::ostringstream omsg;
		omsg
			<< "FeasibilityStepReducedStd_Strategy::compute_feasibility_step(...) : "
			<< "Error, a QP relaxation parameter\n"
			<< "of eta = " << qp_eta << " was calculated and therefore it must be assumed\n"
			<< "that the NLP's constraints are infeasible\n"
			<< "Throwing an InfeasibleConstraints exception!\n";
		if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out << omsg.str();
		}
		throw InfeasibleConstraints(omsg.str());
	}

	//
	// Set the full step
	//
	// w = Ywy + Zwz
	//
	LinAlgPack::V_VpV( w, Ywy, Zwz );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\n||w||inf    = " << LinAlgPack::norm_inf(*w);
		out << std::endl;
	}
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nw = \n" << *w;
	}

	if( throw_qp_failure )
		return false;
	return true;
}

void FeasibilityStepReducedStd_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
	out << L << "*** Computes the step by solving range and null space problems\n"
		<< L << "*** ToDo: Finish documentation!\n";
}

} // end namespace ReducedSpaceSQPPack
