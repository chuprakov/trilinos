// ///////////////////////////////////////////////////////
// ReducedHessianSecantUpdateBFGSProjected_Strategy.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	
// disable VC 5.0 warnings about truncated identifier names (templates).
#pragma warning(disable : 4503)	

#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateBFGSProjected_Strategy.h"
#include "ReducedSpaceSQPPack/include/rSQPAlgo.h"
#include "ReducedSpaceSQPPack/include/rSQPState.h"
#include "ConstrainedOptimizationPack/include/MatrixSymAddDelUpdateable.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/include/MatrixWithOpOut.h"
#include "SparseLinAlgPack/include/GenPermMatrixSlice.h"
#include "SparseLinAlgPack/include/GenPermMatrixSliceOp.h"
#include "SparseLinAlgPack/include/MatrixSymInitDiagonal.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "Misc/include/dynamic_cast_verbose.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

namespace ReducedSpaceSQPPack {

ReducedHessianSecantUpdateBFGSProjected_Strategy::ReducedHessianSecantUpdateBFGSProjected_Strategy(
	const bfgs_update_ptr_t&      bfgs_update
	,value_type                   proj_start_act_set_frac
	,value_type                   super_basic_mult_drop_tol
	)
	: bfgs_update_(bfgs_update)
	, proj_start_act_set_frac_(proj_start_act_set_frac)
	, super_basic_mult_drop_tol_(super_basic_mult_drop_tol)
{}

bool ReducedHessianSecantUpdateBFGSProjected_Strategy::perform_update(
	VectorSlice* s_bfgs, VectorSlice* y_bfgs, bool first_update
	,std::ostream& out, EJournalOutputLevel olevel, rSQPAlgo *algo, rSQPState *s
	,MatrixWithOp *rHL_k
	)
{
	using std::setw;
	using std::endl;
	using std::right;
	using DynamicCastHelperPack::dyn_cast;
	using LinAlgOpPack::V_MtV;
	using SparseLinAlgPack::norm_inf;
	typedef ConstrainedOptimizationPack::MatrixHessianSuperBasic MHSB_t;

#ifdef _WINDOWS
	MHSB_t &rHL_super = dynamic_cast<MHSB_t&>(*rHL_k);
#else
	MHSB_t &rHL_super = dyn_cast<MHSB_t>(*rHL_k);
#endif

	const GenPermMatrixSlice
		&Q_R = rHL_super.Q_R(),
		&Q_X = rHL_super.Q_X();
	const size_type
		r    = algo->nlp().r(),
		n_pz = Q_R.rows();

	bool do_projected_rHL_RR = false;

	if( Q_R.is_identity() ) {
		// Determine when to start adding and removing rows/cols form rHL_RR
		if( act_set_stats_(*s).updated_k(-1) ) {
			if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
				out	<< "\nDetermining if projected BFGS updating of superbasics should begin ...\n";
			}
			ActSetStats &stats = act_set_stats_(*s).get_k(-1);
			const SpVector
				&nu_km1 = s->nu().get_k(-1);
			const SpVectorSlice
				nu_indep = nu_km1(s->var_indep());
			const size_type
				num_active_indep = nu_indep.nz(),
				num_adds         = stats.num_adds(),
				num_drops        = stats.num_drops();
			const value_type
				frac_same
				= ( num_adds == ActSetStats::NOT_KNOWN || num_active_indep == 0
					? 0.0
					: std::_MAX(((double)(num_active_indep)-num_adds-num_drops) / num_active_indep, 0.0 ) );
			do_projected_rHL_RR = ( num_active_indep > 0 && frac_same >= proj_start_act_set_frac() );
			if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
				out << "\nnum_active_indep = " << num_active_indep;
				if( num_active_indep ) {
					out	<< "\nmax(num_active_indep-num_adds-num_drops,0)/(num_active_indep) = "
						<< "max("<<num_active_indep<<"-"<<num_adds<<"-"<<num_drops<<",0)/("<<num_active_indep<<") = "
						<< frac_same;
					if( do_projected_rHL_RR )
						out << " >= ";
					else
						out << " < ";
					out << "proj_start_act_set_frac = " << proj_start_act_set_frac();
				}
				if( do_projected_rHL_RR )
					out << "\nStart performing projected BFGS updating of superbasic variables only!\n";
				else
					out << "\nJust keep BFGS updating the full reduced Hessian rHL for now!\n";
			}
			if( do_projected_rHL_RR ) {
				//
				// Eliminate those rows/cols from rHL_RR for fixed variables and reinitialize rHL
				//
				if( i_x_free_.size() < n_pz ) { // Only need to resize these once
					i_x_free_.resize(n_pz); 
					i_x_fixed_.resize(n_pz);
					bnd_fixed_.resize(n_pz);
				}
				// Loop through and set i_x_free and i_x_fixed
				if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
					out << "\nDetermining which fixed variables to remove from rHL_RR (can remove all but one)...\n";
				}
				const value_type
					max_nu_indep = norm_inf(nu_indep);
				const bool
					all_fixed = n_pz == nu_indep.nz();
				if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
					out << "\nmax{|nu_k(indep)|,i=r+1...n} = " << max_nu_indep << std::endl;
				}
				const int prec = out.precision();
				if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ACTIVE_SET) ) {
					out << endl
						<< right << setw(10)      << "i"
						<< right << setw(prec+12) << "nu(i)"
						<< right << setw(8)       << "status"
						<< endl
						<< right << setw(10)      << "--------"
						<< right << setw(prec+12) << "--------"
						<< right << setw(8)       << "------"
						<< endl;
				}
				SpVector::const_iterator
					nu_itr = nu_indep.begin(),
					nu_end = nu_indep.end();
				SpVector::difference_type
					nu_o = nu_indep.offset();
				i_x_t::iterator
					i_x_free_itr  = i_x_free_.begin(),
					i_x_fixed_itr = i_x_fixed_.begin();
				bnd_fixed_t::iterator
					bnd_fixed_itr = bnd_fixed_.begin();
				size_type
					n_pz_X = 0, // We will count these
					n_pz_R = 0;
				bool kept_one = false;
                {for( size_type i_indep = 1; i_indep <= n_pz; ++i_indep ) {
					if( nu_itr != nu_end && (nu_itr->indice() + nu_o) == i_indep ) {
						const value_type
							abs_val = ::fabs(nu_itr->value());
						const bool
							keep = ( (all_fixed && abs_val == max_nu_indep && !kept_one)
									 || abs_val > super_basic_mult_drop_tol() );
						if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ACTIVE_SET) ) {
							out << right << setw(10)      << i_indep + r
								<< right << setw(prec+12) << nu_itr->value()
								<< right << setw(8)       << (keep ? "keep" : "drop")
								<< endl;
						}
						if(!keep) {
							*i_x_fixed_itr++ = i_indep;
							*bnd_fixed_itr++
								= ( nu_itr->value() > 0.0 ? MHSB_t::UPPER : MHSB_t::LOWER );
							// ToDo: Consider fixed variable bounds
							++n_pz_X;
						}
						else {
							kept_one = true;
						}
						++nu_itr;
						if(!keep) continue;
					}
					*i_x_free_itr++ = i_indep;
					++n_pz_R;
				}}
				assert( i_x_free_itr  - i_x_free_.begin()  == n_pz_R );
				assert( i_x_fixed_itr - i_x_fixed_.begin() == n_pz_X );
				assert( bnd_fixed_itr - bnd_fixed_.begin() == n_pz_X );
				// Delete rows/cols from rHL_RR for fixed variables
#ifdef _WINDOWS
				MatrixSymAddDelUpdateable
					&rHL_RR = dynamic_cast<MatrixSymAddDelUpdateable&>(
						const_cast<MatrixSymWithOpFactorized&>(*rHL_super.B_RR_ptr())
						);
#else
				MatrixSymAddDelUpdateable
					&rHL_RR = dyn_cast<MatrixSymAddDelUpdateable>(
						const_cast<MatrixSymWithOpFactorized&>(*rHL_super.B_RR_ptr())
						);
#endif
				if( n_pz_R < n_pz ) {
					if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
						out << "\nDeleting n_pz_X = " << n_pz_X << " rows/columns from rHL_RR for fixed independent variables...\n";
					}
					{for( size_type k = 0; k < n_pz_X; ++k ) {
						rHL_RR.delete_update( i_x_fixed_[k]-k, false );
					}}
					assert( rHL_RR.rows() == n_pz_R );
					if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
						out << "\nrHL_RR after variables where removed =\n" << *rHL_super.B_RR_ptr();
					}
					// Initialize B_XX = I for now
#ifdef _WINDOWS
					MatrixSymInitDiagonal
						&B_XX = dynamic_cast<MatrixSymInitDiagonal&>(
							const_cast<MatrixSymWithOp&>(*rHL_super.B_XX_ptr())
							);
#else
					MatrixSymInitDiagonal
						&B_XX = dyn_cast<MatrixSymInitDiagonal>(
							const_cast<MatrixSymWithOp&>(*rHL_super.B_XX_ptr())
							);
#endif
					B_XX.init_identity( n_pz_X ); // ToDo: We may want scaled element?
					// Reinitialize rHL for new active set
					rHL_super.initialize(
						n_pz, n_pz_R, &i_x_free_[0], &i_x_fixed_[0], &bnd_fixed_[0]
						,rHL_super.B_RR_ptr(),NULL,BLAS_Cpp::no_trans,rHL_super.B_XX_ptr()
						);
					if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
						out << "\nFull rHL after reinitialization but before BFGS update =\n" << *rHL_k;
					}
				}
				else {
					do_projected_rHL_RR = false;
					if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
						out << "\nThere where no variable bounds with |nu(i)| large enough to consider dropping\n";
					}
				}
			}
		}
	}
	else {
		// Modify rHL_RR by adding and dropping rows/cols for freeded and fixed variables
		do_projected_rHL_RR = true;

		// ToDo: Implement this!
	}

	// Perform the BFGS update
	if( do_projected_rHL_RR ) {
		// Perform BFGS update on smaller rHL_RR.
		// By the time we get here rHL_RR should be resize and ready to update
		const GenPermMatrixSlice
			&Q_R = rHL_super.Q_R(),
			&Q_X = rHL_super.Q_X();
		const size_type
			n_pz_R = Q_R.cols(),
			n_pz_X = Q_X.cols();
		assert( n_pz_R + n_pz_X == n_pz );
		// Get projected BFGS update vectors y_bfgs_R, s_bfgs_R
		Vector y_bfgs_R; // y_bfgs_R = Q_R'*y_bfgs
		V_MtV( &y_bfgs_R, Q_R, BLAS_Cpp::trans, *y_bfgs );
		Vector s_bfgs_R; // s_bfgs_R = Q_R'*s_bfgs
		V_MtV( &s_bfgs_R, Q_R, BLAS_Cpp::trans, *s_bfgs );
		// Update rHL_RR
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out << "\nPerform BFGS update on " << n_pz_R << " x " << n_pz_R << " projected reduced Hessian for the superbasic variables where B = rHL_RR...\n";
		}
		bfgs_update().perform_update(
			&s_bfgs_R(),&y_bfgs_R(),first_update,out,olevel,algo->algo_cntr().check_results()
			,const_cast<MatrixWithOp*>(static_cast<const MatrixWithOp*>(rHL_super.B_RR_ptr().get()))
			,&quasi_newton_stats_(*s).set_k(0)
			);
	}
	else {
		// Update the full reduced Hessain matrix (rHL = rHL_RR)
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out << "\nPerform BFGS update on the full reduced Hessian where B = rHL...\n";
		}
		bfgs_update().perform_update(
			s_bfgs,y_bfgs,first_update,out,olevel,algo->algo_cntr().check_results()
			,const_cast<MatrixWithOp*>(static_cast<const MatrixWithOp*>(rHL_super.B_RR_ptr().get()))
			,&quasi_newton_stats_(*s).set_k(0)
			);
	}

	return true;
}

void ReducedHessianSecantUpdateBFGSProjected_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Perform BFGS update on only free independent (super basic) variables.\n"
		<< L << "ToDo: Finish implementation!\n";
}

}  // end namespace ReducedSpaceSQPPack
