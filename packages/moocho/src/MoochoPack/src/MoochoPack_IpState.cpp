// ////////////////////////////////////////////////////////////////////
// ipState.cpp
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
#include "ReducedSpaceSQPPack/include/ipState.h"
#include "AbstractLinAlgPack/include/MatrixSymWithOp.h"

// Iteration Quantity Strings
extern const std::string ReducedSpaceSQPPack::Vu_name = "Vu";
extern const std::string ReducedSpaceSQPPack::Vl_name = "Vl";
extern const std::string ReducedSpaceSQPPack::invXu_name = "invXu";
extern const std::string ReducedSpaceSQPPack::invXl_name = "invXl";
extern const std::string ReducedSpaceSQPPack::rHB_name = "rHB";
extern const std::string ReducedSpaceSQPPack::B_name = "B";

namespace ReducedSpaceSQPPack {

ipState::ipState(
  const decomp_sys_ptr_t& decomp_sys   = MemMngPack::null
  ,const vec_space_ptr_t& space_x      = MemMngPack::null
  ,const vec_space_ptr_t& space_c      = MemMngPack::null
  ,const vec_space_ptr_t& space_h      = MemMngPack::null
  ,const vec_space_ptr_t& space_range  = MemMngPack::null
  ,const vec_space_ptr_t& space_null   = MemMngPack::null
  )
	:
	rSQPState(decomp_sys, space_x, space_c, space_h, space_range, space_null)
	{
	}

ipState::~ipState()
	{
	}

///********** Iteration Quantities **************

STATE_IQ_DEF(ipState, MatrixSymWithOp, Vu, Vu_name);

STATE_IQ_DEF(ipState, MatrixSymWithOp, Vl, Vl_name);

STATE_IQ_DEF(ipState, MatrixSymWithOp, invXu, invXu_name);

STATE_IQ_DEF(ipState, MatrixSymWithOp, invXl, invXl_name);

STATE_IQ_DEF(ipState, MatrixSymWithOp, rHB, rHB_name);

STATE_IQ_DEF(ipState, MatrixSymWithOp, B, B_name);

} // end namespace ReducedSpaceSQPPack
