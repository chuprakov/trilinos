// ////////////////////////////////////////////////////////////////////////////////////
// ipState.h
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

#include "ReducedSpaceSQPPack/include/rSQPState.h"

namespace ReducedSpaceSQPPack {

// Iteration Quantity Strings
extern const std::string Vu_name;
extern const std::string Vl_name;
extern const std::string invXu_name;
extern const std::string invXl_name;
extern const std::string rHB_name;
extern const std::string B_name;


class ipState 
	: public ReducedSpaceSQPPack::rSQPState
	{

		///********** Iteration Quantities **************

		/// Vu - diagonal matrix of upper bound multipliers
		STATE_IQ_DECL(MatrixSymWithOp, Vu);

		/// Vl - diagonal matrix of lower bound multipliers
		STATE_IQ_DECL(MatrixSymWithOp, Vl);

		/// invXu - (Xu)^-1 - matrix of 1/(xu-x) diagonal
		STATE_IQ_DECL(MatrixSymWithOp, invXu);

		/// invXl - (Xl)^-1 - matrix of 1/(x-xl) diagonal
		STATE_IQ_DECL(MatrixSymWithOp, invXl);

		/// rHB - reduced Hessian of the barrier term (Z_Sigma_Z)
		STATE_IQ_DECL(MatrixSymWithOp, rHB);

		/// B - overall reduced 'Hessian' (Z_W_Z+Z_Sigma_Z)
		STATE_IQ_DECL(MatrixSymWithOp, B);


		///
		/** Construct
		 *
		 * 
		 */
		ipState(
		  const decomp_sys_ptr_t& decomp_sys   = MemMngPack::null
		  ,const vec_space_ptr_t& space_x      = MemMngPack::null
		  ,const vec_space_ptr_t& space_c      = MemMngPack::null
		  ,const vec_space_ptr_t& space_h      = MemMngPack::null
		  ,const vec_space_ptr_t& space_range  = MemMngPack::null
		  ,const vec_space_ptr_t& space_null   = MemMngPack::null
		  );

		virtual ~ipState();

	}; // end class ipState

} // end namespace ReducedSpaceSQPPack



#endif // if !defined IP_STATE_H
