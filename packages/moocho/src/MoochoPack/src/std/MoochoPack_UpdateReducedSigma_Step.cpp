// ////////////////////////////////////////////////////////////////////////////
// UpdateReducedSigma_Step.cpp
//
// Copyright (C) 2001
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
#include <typeinfo>
#include <iostream>

#include "AbstractLinAlgPack/include/MatrixSymDiagonalStd.h"
#include "AbstractLinAlgPack/include/MatrixWithOp.h"
#include "AbstractLinAlgPack/include/MatrixSymWithOpNonsingular.h"
#include "AbstractLinAlgPack/include/MultiVectorMutable.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/VectorWithOpOut.h"
#include "AbstractLinAlgPack/include/VectorAuxiliaryOps.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "AbstractLinAlgPack/include/assert_print_nan_inf.h"
#include "ConstrainedOptimizationPack/include/MatrixIdentConcat.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "NLPInterfacePack/include/NLPFirstOrderInfo.h"
#include "ReducedSpaceSQPPack/include/ipState.h"
#include "ReducedSpaceSQPPack/include/std/UpdateReducedSigma_Step.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"

#include "StringToIntMap.h"

#include "dynamic_cast_verbose.h"

namespace ReducedSpaceSQPPack {

		/** Constructor.
		 */
UpdateReducedSigma_Step::UpdateReducedSigma_Step(
  const e_update_methods update_method
  )
	:
	update_method_(update_method)
	{}

bool UpdateReducedSigma_Step::do_step(
  Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
	{
	using DynamicCastHelperPack::dyn_cast;
	using GeneralIterationPack::print_algorithm_step;

	rSQPAlgo            &algo   = dyn_cast<rSQPAlgo>(_algo);
	ipState             &s      = dyn_cast<ipState>(_algo.state());
	NLPFirstOrderInfo   &nlp    = dyn_cast<NLPFirstOrderInfo>(algo.nlp());
	
	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();
	
	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
		{
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
		}

	switch (update_method_)
		{
		case ALWAYS_EXPLICIT:
			{
			FormReducedSigmaExplicitly(_algo);
			if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
				{
				out << "update_method is always_explicit\n"
					<< " forming Reduced Sigma Explicitly\n";
				}
			}
			break;
		case BFGS_PRIMAL:
		case BFGS_DUAL_NO_CORRECTION:
		case BFGS_DUAL_EXPLICIT_CORRECTION:
		case BFGS_DUAL_SCALING_CORRECTION:
			{
			THROW_EXCEPTION(true, std::logic_error, "Option BFGS_primal not handled yet in UpdateReducedSigma\n");
			}
			break;
		default:
			assert(0); // local error ?
		};

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) 
		{
		out << "rHB=\n";
		s.rHB().get_k(0).output(out);
		}

	return true;
	}


void UpdateReducedSigma_Step::print_step(
  const Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
	{
	out << L << "\n#Update Z^T*Sigma*Z\n"
		<< L << "if (update_method == always_explicit) then\n"
		<< L << "   sigma = invXl*Vl-invXu*Vu\n"
		<< L << "   sigma_I = sigma.sub_view(Z.I_rng)\n"
		<< L << "   sigma_D_sqrt = (sigma.sub_view(Z.D_rng))^1/2\n"
		<< L << "   J = sigma_D_sqrt*Z\n"
		<< L << "   rHB = J^T*J + sigma_I\n"
		<< L << "elsif (update_method == BFGS_???) then\n"
		<< L << "   # NOT IMPLEMENTED YET\n"
		<< L << "end\n";
	}


void UpdateReducedSigma_Step::FormReducedSigmaExplicitly(Algorithm& _algo)
	{
	using DynamicCastHelperPack::dyn_cast;
 	using LinAlgOpPack::Mp_MtM;
 	using LinAlgOpPack::M_MtM;
 	using LinAlgOpPack::M_StM;

	// Calculate Reduced Sigma directly from
	// Sigma = invXl*Vl + invXu*Vu
	// Z_kT*Sigma*Z_k

	rSQPAlgo            &algo   = dyn_cast<rSQPAlgo>(_algo);
	ipState             &s      = dyn_cast<ipState>(_algo.state());

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// Get the iteration quantities
	const MatrixIdentConcat& Z = dyn_cast<MatrixIdentConcat>(s.Z().get_k(0));
	const MatrixSymDiagonalStd& invXl = s.invXl().get_k(0);
	const MatrixSymDiagonalStd& invXu = s.invXu().get_k(0);
	const MatrixSymDiagonalStd& Vl = s.Vl().get_k(0);
	const MatrixSymDiagonalStd& Vu = s.Vu().get_k(0);
	
	MatrixSymDiagonalStd& sigma = s.Sigma().set_k(0);

	MatrixSymWithOpNonsingular& rHB = dyn_cast<MatrixSymWithOpNonsingular>(s.rHB().set_k(0));
	if (rHB.cols() != Z.cols())
		{
		// Initialize space in rHB
		dyn_cast<MatrixSymInitDiagonal>(rHB).init_identity(Z.space_rows(), 0.0);
		}
	
	// Calculate Sigma
	sigma.diag()=0;
	ele_wise_prod(1.0, invXl.diag(), Vl.diag(), &sigma.diag());
	
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
		{
		out << "sigma_l = \n"
			<< sigma.diag();
			
		}
	
    MemMngPack::ref_count_ptr<VectorWithOpMutable> temp = invXl.diag().space().create_member(0);
	ele_wise_prod(1.0, invXu.diag(), Vu.diag(), temp.get());
	sigma.diag().axpy(1.0, *temp);

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
		{
		out << "sigma_l + sigma_u = \n"
			<< sigma.diag();
			
		}

	// Calculate the cross term first...
	VectorWithOpMutable& w_sigma = s.w_sigma().set_k(0);
	*temp = 0;
	ele_wise_prod(1.0, s.Ypy().get_k(0), sigma.diag(), temp.get());

	using LinAlgOpPack::V_MtV;
	V_MtV(&w_sigma, Z, BLAS_Cpp::trans, *temp);

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
		{
		out << "w = \n";
		w_sigma.output(out);
		}
	
	// Calculate Reduced Sigma
	// Try sigma^1/2 making use of dependent and independent variables
	MemMngPack::ref_count_ptr<MatrixSymDiagonalStd> sigma_I_sqrt = 
		MemMngPack::rcp(
		  new MatrixSymDiagonalStd(
			sigma.diag().sub_view(Z.I_rng())->clone()
			)
		  );


	MemMngPack::ref_count_ptr<MatrixSymDiagonalStd> sigma_D_sqrt = 
		MemMngPack::rcp(
		  new MatrixSymDiagonalStd(
			sigma.diag().sub_view(Z.D_rng())->clone()
			)
		  );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
		{
		out << "sigma_I = \n";
		sigma_I_sqrt->output(out);
		}
	
	ele_wise_sqrt(sigma_I_sqrt->diag());
	ele_wise_sqrt(sigma_D_sqrt->diag());

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
		{
		out << "sigma_D_sqrt = \n";
		sigma_D_sqrt->output(out);
		}

	MemMngPack::ref_count_ptr<MultiVectorMutable> J = sigma_D_sqrt->space_cols().create_members(Z.cols());
	M_MtM((MatrixWithOp*)J.get(), *sigma_D_sqrt, BLAS_Cpp::no_trans, Z.D(), BLAS_Cpp::no_trans);

	J->syrk(BLAS_Cpp::trans, 1.0, 0.0, &rHB);

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
		{
		out << "J'*J = \n";
		rHB.output(out);
		}	

	// This is really bad, but I cannot add a diagonal matrix to a PosDefCholFactor one
	MemMngPack::ref_count_ptr<MultiVectorMutable> Ji = sigma_I_sqrt->space_cols().create_members(sigma_I_sqrt->cols());
	using LinAlgOpPack::assign;
	assign(Ji.get(), *sigma_I_sqrt, BLAS_Cpp::no_trans);
	Ji->syrk(BLAS_Cpp::trans, 1.0, 1.0, &rHB);
  
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
		{
		out << "rHB = \n";
		rHB.output(out);
		}	

	//std::cout << "rHB=\n";
 	//rHB.output(std::cout);
	
	/*
 	// Try the full unspecialised calculation, this is expensive, but will
	// serve as a debug for the more efficient calculations.

	VectorSpace::multi_vec_mut_ptr_t sigma_Z = Z.space_cols().create_members(Z.cols());
	M_MtM((MatrixWithOp*)sigma_Z.get(), sigma, BLAS_Cpp::no_trans, Z, BLAS_Cpp::no_trans);

	//std::cout << "sigma_Z\n";
	//sigma_Z->output(std::cout);

	VectorSpace::multi_vec_mut_ptr_t ZT_sigma_Z = Z.space_rows().create_members(Z.cols());
	M_MtM((MatrixWithOp*)ZT_sigma_Z.get(), (MatrixWithOp&)Z, BLAS_Cpp::trans, (MatrixWithOp&)*sigma_Z, BLAS_Cpp::no_trans);

	std::cout << "ZT_sigma_Z=\n";
	ZT_sigma_Z->output(std::cout);
	*/
	}


namespace {

const int local_num_options = 1;

enum local_EOptions 
	{
	UPDATE_METHOD
	};

const char* local_SOptions[local_num_options] = 
	{
	"update_method",
	};

const int num_update_methods = 5;

const char* s_update_methods[num_update_methods] = 
	{
	"always_explicit",
	"BFGS_primal",
	"BFGS_dual_no_correction",
	"BFGS_dual_explicit_correction",
	"BFGS_dual_scaling_correction"
	};

}

 
UpdateReducedSigma_StepSetOptions::UpdateReducedSigma_StepSetOptions(
  UpdateReducedSigma_Step* target
  , const char opt_grp_name[] )
	:
	OptionsFromStreamPack::SetOptionsFromStreamNode(
	  opt_grp_name, local_num_options, local_SOptions ),
	OptionsFromStreamPack::SetOptionsToTargetBase< UpdateReducedSigma_Step >( target )
	{
	}

void UpdateReducedSigma_StepSetOptions::set_option( 
  int option_num, const std::string& option_value )
	{
	using OptionsFromStreamPack::StringToIntMap;

	typedef UpdateReducedSigma_Step target_t;
	switch( (local_EOptions)option_num ) 
		{
		case UPDATE_METHOD:
			{
			StringToIntMap	config_map( UpdateReducedSigma_opt_grp_name, num_update_methods, s_update_methods );
			target().update_method( (UpdateReducedSigma_Step::e_update_methods) config_map( option_value.c_str() ) );
			}
			break;
		default:
			assert(0);	// Local error only?
		}
	}

}; // end namespace ReducedSpaceSQPPack
