// ////////////////////////////////////////////////////////////////////
// rSQPState.cpp
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

#pragma warning(disable : 4786)	// too long class name for debugger warning

#include <sstream>
#include <typeinfo>

#include "ReducedSpaceSQPPack/include/rSQPState.h"
#include "ConstrainedOptimizationPack/include/MeritFuncNLP.h"
#include "AbstractLinAlgPack/include/MatrixSymWithOp.h"
#include "AbstractLinAlgPack/include/MatrixWithOpNonsingular.h"

// rSQPState iteration quantities names

// Iteration Info
const std::string ReducedSpaceSQPPack::num_basis_name		= "num_basis";
// NLP Problem Info 
const std::string ReducedSpaceSQPPack::x_name				= "x";
const std::string ReducedSpaceSQPPack::f_name				= "f";
const std::string ReducedSpaceSQPPack::Gf_name				= "Gf";
const std::string ReducedSpaceSQPPack::HL_name				= "HL";
const std::string ReducedSpaceSQPPack::c_name				= "c";
const std::string ReducedSpaceSQPPack::h_name				= "h";
const std::string ReducedSpaceSQPPack::Gc_name				= "Gc";
const std::string ReducedSpaceSQPPack::Gh_name				= "Gh";
// Constraint Gradient Null Space / Range Space Decomposition Info
const std::string ReducedSpaceSQPPack::Y_name				= "Y";
const std::string ReducedSpaceSQPPack::Z_name				= "Z";
const std::string ReducedSpaceSQPPack::R_name				= "R";
const std::string ReducedSpaceSQPPack::Uy_name				= "Uy";
const std::string ReducedSpaceSQPPack::Uz_name				= "Uz";
const std::string ReducedSpaceSQPPack::Vy_name				= "Vy";
const std::string ReducedSpaceSQPPack::Vz_name				= "Vz";
// Search Direction Info
const std::string ReducedSpaceSQPPack::py_name				= "py";
const std::string ReducedSpaceSQPPack::Ypy_name				= "Ypy";
const std::string ReducedSpaceSQPPack::pz_name				= "pz";
const std::string ReducedSpaceSQPPack::Zpz_name				= "Zpz";
const std::string ReducedSpaceSQPPack::d_name				= "d";
// Reduced QP Subproblem Info
const std::string ReducedSpaceSQPPack::rGf_name				= "rGf";
const std::string ReducedSpaceSQPPack::rHL_name				= "rHL";
const std::string ReducedSpaceSQPPack::w_name				= "w";
const std::string ReducedSpaceSQPPack::zeta_name			= "zeta";
const std::string ReducedSpaceSQPPack::qp_grad_name			= "qp_grad";
const std::string ReducedSpaceSQPPack::eta_name				= "eta";
// Global Convergence Info
const std::string ReducedSpaceSQPPack::alpha_name			= "alpha";
const std::string ReducedSpaceSQPPack::merit_func_nlp_name	= "merit_func_nlp";
const std::string ReducedSpaceSQPPack::mu_name				= "mu";
const std::string ReducedSpaceSQPPack::phi_name				= "phi";
// KKT Info
const std::string ReducedSpaceSQPPack::opt_kkt_err_name		= "opt_kkt_err";
const std::string ReducedSpaceSQPPack::feas_kkt_err_name	= "feas_kkt_err";
const std::string ReducedSpaceSQPPack::GL_name				= "GL";
const std::string ReducedSpaceSQPPack::rGL_name				= "rGL";
const std::string ReducedSpaceSQPPack::lambda_name			= "lambda";
const std::string ReducedSpaceSQPPack::lambdaI_name			= "lambdaI";
const std::string ReducedSpaceSQPPack::nu_name				= "nu";

namespace ReducedSpaceSQPPack {

// rSQPState

rSQPState::rSQPState(
	const decomp_sys_ptr_t& decomp_sys
	,const space_x_ptr_t&   space_x
	,const space_x_ptr_t&   space_c
	,const space_x_ptr_t&   space_h
	,const space_x_ptr_t&   space_range
	,const space_x_ptr_t&   space_null
	)
	:decomp_sys_(decomp_sys)
	,space_x_(space_x)
	,space_c_(space_c)
	,space_h_(space_h)
	,space_range_(space_range)
	,space_null_(space_null)
{}

// Iteration Info

RSQP_STATE_INDEX_IQ_DEF(  rSQPState,                  num_basis, num_basis_name          )

// NLP Problem Info

RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  x,         x_name,  get_space_x()  )
RSQP_STATE_SCALAR_IQ_DEF( rSQPState,                  f,         f_name                  )
RSQP_STATE_IQ_DEF(        rSQPState, MatrixSymWithOp, HL,        HL_name                 )
RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  Gf,        Gf_name, get_space_x()  )
RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  c,         c_name,  get_space_c()  )
RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  h,         h_name,  get_space_h()  )
RSQP_STATE_IQ_DEF(        rSQPState, MatrixWithOp,    Gc,        Gc_name                 )
RSQP_STATE_IQ_DEF(        rSQPState, MatrixWithOp,    Gh,        Gh_name                 )

// Constraint Gradient Null Space / Range Space Decomposition Info

RSQP_STATE_IQ_DEF(        rSQPState, MatrixWithOp,            Y,  Y_name                  )
RSQP_STATE_IQ_DEF(        rSQPState, MatrixWithOp,            Z,  Z_name                  )
RSQP_STATE_IQ_DEF(        rSQPState, MatrixWithOpNonsingular, R,  R_name                  )
RSQP_STATE_IQ_DEF(        rSQPState, MatrixWithOp,            Uy, Uy_name                 )
RSQP_STATE_IQ_DEF(        rSQPState, MatrixWithOp,            Uz, Uz_name                 )
RSQP_STATE_IQ_DEF(        rSQPState, MatrixWithOp,            Vy, Vy_name                 )
RSQP_STATE_IQ_DEF(        rSQPState, MatrixWithOp,            Vz, Vz_name                 )

// Search Direction Info

RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  py,  py_name,   get_space_range() )
RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  Ypy, Ypy_name,  get_space_x()     )
RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  pz,  pz_name,   get_space_null()  )
RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  Zpz, Zpz_name,  get_space_x()     )
RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  d,   d_name,    get_space_x()     )

// QP Subproblem Info

RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  rGf,     rGf_name,      get_space_null() )
RSQP_STATE_IQ_DEF(        rSQPState, MatrixSymWithOp, rHL,     rHL_name                        )
RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  w,       w_name,        get_space_null() ) 
RSQP_STATE_SCALAR_IQ_DEF( rSQPState,                  zeta,    zeta_name                       )
RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  qp_grad, qp_grad_name,  get_space_null() )
RSQP_STATE_SCALAR_IQ_DEF( rSQPState,                  eta,     eta_name                        )

// Global Convergence Info

RSQP_STATE_SCALAR_IQ_DEF( rSQPState,                  alpha,          alpha_name          )
RSQP_STATE_IQ_DEF(        rSQPState, MeritFuncNLP,    merit_func_nlp, merit_func_nlp_name )
RSQP_STATE_SCALAR_IQ_DEF( rSQPState,                  mu,             mu_name             )
RSQP_STATE_SCALAR_IQ_DEF( rSQPState,                  phi,            phi_name            )

// KKT Info

RSQP_STATE_SCALAR_IQ_DEF( rSQPState,                  opt_kkt_err,    opt_kkt_err_name                     )
RSQP_STATE_SCALAR_IQ_DEF( rSQPState,                  feas_kkt_err,   feas_kkt_err_name                    )
RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  GL,             GL_name,           get_space_x()     )
RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  rGL,            rGL_name,          get_space_null()  )
RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  lambda,         lambda_name,       get_space_c()     )
RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  lambdaI,        lambdaI_name,      get_space_h()     )
RSQP_STATE_VECTOR_IQ_DEF( rSQPState,                  nu,             nu_name,           get_space_x()     )

// protected

void rSQPState::update_iq_id(
	const std::string&                iq_name
	,iq_id_encap*                     iq_id
	) const
{
	namespace rcp = ReferenceCountingPack;
	if(iq_id->iq_id == DOES_NOT_EXIST)
		iq_id->iq_id = this->get_iter_quant_id(iq_name);
	THROW_EXCEPTION(
		iq_id->iq_id == DOES_NOT_EXIST, DoesNotExist
		,"rSQPState::update_iq_id(iq_name,iq_id) : Error, "
		" The iteration quantity with name \'" << iq_name <<
		"\' does not exist!" );
}

void rSQPState::update_index_type_iq_id(
	const std::string&                iq_name
	,iq_id_encap*                     iq_id
	)
{
	namespace rcp = ReferenceCountingPack;
	if(iq_id->iq_id == DOES_NOT_EXIST) {
		iq_id_type
			_iq_id = this->get_iter_quant_id(iq_name);
		if(_iq_id == DOES_NOT_EXIST) {
			iq_id->iq_id = this->set_iter_quant(
				iq_name
				,rcp::rcp(
					new IterQuantityAccessContiguous<index_type>(
						1
						,iq_name
						)
					)
				);
		}
		else {
			iq_id->iq_id = _iq_id;
		}
	}
}

void rSQPState::update_value_type_iq_id(
	const std::string&                iq_name
	,iq_id_encap*                     iq_id
	)
{
	namespace rcp = ReferenceCountingPack;
	if(iq_id->iq_id == DOES_NOT_EXIST) {
		iq_id_type
			_iq_id = this->get_iter_quant_id(iq_name);
		if(_iq_id == DOES_NOT_EXIST) {
			iq_id->iq_id = this->set_iter_quant(
				iq_name
				,rcp::rcp(
					new IterQuantityAccessContiguous<value_type>(
						1
						,iq_name
						)
					)
				);
		}
		else {
			iq_id->iq_id = _iq_id;
		}
	}
}

void rSQPState::update_vector_iq_id(
	const std::string&                iq_name
	,const VectorSpace::space_ptr_t&  vec_space
	,iq_id_encap*                     iq_id
	)
{
	namespace rcp = ReferenceCountingPack;
	if(iq_id->iq_id == DOES_NOT_EXIST) {
		iq_id_type
			_iq_id = this->get_iter_quant_id(iq_name);
		if(_iq_id == DOES_NOT_EXIST) {
			iq_id->iq_id = this->set_iter_quant(
				iq_name
				,rcp::rcp(
					new IterQuantityAccessContiguous<VectorWithOpMutable>(
						1
						,iq_name
						,vec_space
						)
					)
				);
		}
		else {
			iq_id->iq_id = _iq_id;
		}
	}
}

}	// end namespace ReducedSpaceSQPPack
