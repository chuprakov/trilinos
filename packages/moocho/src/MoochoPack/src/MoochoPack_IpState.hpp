// ////////////////////////////////////////////////////////////////////////////////////
// MoochoPack_IpState.hpp
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
//

#if !defined IP_STATE_H
#define IP_STATE_H

#include "MoochoPack_NLPAlgoState.hpp"

namespace MoochoPack {

// Iteration Quantity Strings
extern const std::string barrier_parameter_name;
extern const std::string barrier_obj_name;
extern const std::string grad_barrier_obj_name;
extern const std::string e_tol_name;
extern const std::string comp_err_mu_name;
extern const std::string Vu_name;
extern const std::string Vl_name;
extern const std::string invXu_name;
extern const std::string invXl_name;
extern const std::string rHB_name;
extern const std::string B_name;
extern const std::string Sigma_name;
extern const std::string w_sigma_name;
extern const std::string dvl_name;
extern const std::string dvu_name;
extern const std::string alpha_vl_name;
extern const std::string alpha_vu_name;


class IpState 
	: public MoochoPack::NLPAlgoState
	{

	public:
		///********** Iteration Quantities **************

		/// mu: barrier parameter
		STATE_SCALAR_IQ_DECL(barrier_parameter)

		/// barrier_obj: objective value with 
		//   barrier term included
		STATE_SCALAR_IQ_DECL(barrier_obj)

		/// grad_barrier_obj: gradient of the objective
		//   with barrier term included
		STATE_VECTOR_IQ_DECL(grad_barrier_obj)

		/// e_tol: current error tolerance for inner loop
		STATE_SCALAR_IQ_DECL(e_tol)

		/// comp_err_mu: perturbed complementarity error for barrier sub problem
		STATE_SCALAR_IQ_DECL(comp_err_mu)

		/// Vu - diagonal matrix of upper bound multipliers
		STATE_IQ_DECL(MatrixSymDiagStd, Vu)

		/// Vl - diagonal matrix of lower bound multipliers
		STATE_IQ_DECL(MatrixSymDiagStd, Vl)

		/// invXu - (Xu)^-1 - matrix of 1/(xu-x) diagonal
		STATE_IQ_DECL(MatrixSymDiagStd, invXu)

		/// invXl - (Xl)^-1 - matrix of 1/(x-xl) diagonal
		STATE_IQ_DECL(MatrixSymDiagStd, invXl)

		/// rHB - reduced Hessian of the barrier term (Z_Sigma_Z)
		STATE_IQ_DECL(MatrixSymOp, rHB)

		/// B - overall reduced 'Hessian' (Z_W_Z+Z_Sigma_Z)
		STATE_IQ_DECL(MatrixSymOp, B)

		/// Full space Sigma (invXl*Vl-invXu*Vu)
		STATE_IQ_DECL(MatrixSymDiagStd, Sigma)

		/// w_sigma:  crossterm correction for sigma (Z' * Sigma * Y * py)
		STATE_VECTOR_IQ_DECL(w_sigma) 

		/// dvl:  Search direction for lower bound multipliers ( n x 1 )
		STATE_VECTOR_IQ_DECL(dvl)

		/// dvu:  Search direction for upper bound multipliers ( n x 1 )
		STATE_VECTOR_IQ_DECL(dvu)

		/// alpha_vl: step size for vl
		STATE_SCALAR_IQ_DECL(alpha_vl)

		/// alpha_vl: step size for vu
		STATE_SCALAR_IQ_DECL(alpha_vu)

		///
		/** Construct
		 *
		 * 
		 */
		IpState(
		  const decomp_sys_ptr_t& decomp_sys   = Teuchos::null
		  ,const vec_space_ptr_t& space_x      = Teuchos::null
		  ,const vec_space_ptr_t& space_c      = Teuchos::null
		  ,const vec_space_ptr_t& space_range  = Teuchos::null
		  ,const vec_space_ptr_t& space_null   = Teuchos::null
		  );

		virtual ~IpState();

	}; // end class IpState

} // end namespace MoochoPack



#endif // if !defined IP_STATE_H
