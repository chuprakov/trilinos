// ////////////////////////////////////////////////////////////////////
// IpState.cpp
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
#include "MoochoPack_IpState.hpp"
#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"

// Iteration Quantity Strings
extern const std::string MoochoPack::barrier_parameter_name = "barrier_parameter";
extern const std::string MoochoPack::barrier_obj_name = "barrier_obj";
extern const std::string MoochoPack::grad_barrier_obj_name = "grad_barrier_obj";
extern const std::string MoochoPack::e_tol_name = "e_tol";
extern const std::string MoochoPack::comp_err_mu_name = "comp_err_mu";
extern const std::string MoochoPack::Vu_name = "Vu";
extern const std::string MoochoPack::Vl_name = "Vl";
extern const std::string MoochoPack::invXu_name = "invXu";
extern const std::string MoochoPack::invXl_name = "invXl";
extern const std::string MoochoPack::rHB_name = "rHB";
extern const std::string MoochoPack::B_name = "B";
extern const std::string MoochoPack::Sigma_name = "Sigma";
extern const std::string MoochoPack::w_sigma_name = "w_sigma";
extern const std::string MoochoPack::dvl_name = "dvl";
extern const std::string MoochoPack::dvu_name = "dvu";
extern const std::string MoochoPack::alpha_vl_name = "alpha_vl";
extern const std::string MoochoPack::alpha_vu_name = "alpha_vu";

namespace MoochoPack {

IpState::IpState(
  const decomp_sys_ptr_t& decomp_sys
  ,const vec_space_ptr_t& space_x
  ,const vec_space_ptr_t& space_c
  ,const vec_space_ptr_t& space_range
  ,const vec_space_ptr_t& space_null
  )
	:
	NLPAlgoState(decomp_sys, space_x, space_c, space_range, space_null)
	{
	}

IpState::~IpState()
	{
	}

///********** Iteration Quantities **************

STATE_SCALAR_IQ_DEF(IpState, barrier_parameter, barrier_parameter_name)

STATE_SCALAR_IQ_DEF(IpState, barrier_obj, barrier_obj_name)

STATE_VECTOR_IQ_DEF(IpState, grad_barrier_obj, grad_barrier_obj_name, get_space_x(), VST_SPACE_X)

STATE_SCALAR_IQ_DEF(IpState, e_tol, e_tol_name)

STATE_SCALAR_IQ_DEF(IpState, comp_err_mu, comp_err_mu_name)

STATE_IQ_DEF(IpState, MatrixSymDiagStd, Vu, Vu_name)

STATE_IQ_DEF(IpState, MatrixSymDiagStd, Vl, Vl_name)

STATE_IQ_DEF(IpState, MatrixSymDiagStd, invXu, invXu_name)

STATE_IQ_DEF(IpState, MatrixSymDiagStd, invXl, invXl_name)

STATE_IQ_DEF(IpState, MatrixSymOp, rHB, rHB_name)

STATE_IQ_DEF(IpState, MatrixSymOp, B, B_name)

STATE_IQ_DEF(IpState, MatrixSymDiagStd, Sigma, Sigma_name)

STATE_VECTOR_IQ_DEF(IpState, w_sigma, w_sigma_name, get_space_null(), VST_SPACE_NULL )  
STATE_VECTOR_IQ_DEF(IpState, dvl, dvl_name, get_space_x(), VST_SPACE_X)
STATE_VECTOR_IQ_DEF(IpState, dvu, dvu_name, get_space_x(), VST_SPACE_X)

STATE_SCALAR_IQ_DEF(IpState, alpha_vl, alpha_vl_name)
STATE_SCALAR_IQ_DEF(IpState, alpha_vu, alpha_vu_name)
} // end namespace MoochoPack
