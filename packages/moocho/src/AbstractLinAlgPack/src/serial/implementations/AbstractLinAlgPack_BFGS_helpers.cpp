// ///////////////////////////////////////////////////////////
// BFGS_helpers.cpp
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

#include <math.h>
#include <limits>
#include <ostream>

#include "AbstractLinAlgPack_BFGS_helpers.hpp"
#include "AbstractLinAlgPack_Vector.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"

bool AbstractLinAlgPack::BFGS_sTy_suff_p_d(
	const Vector    &s
	,const Vector   &y
	,const value_type     *sTy_in
	,std::ostream         *out
	,const char           func_name[]
	)
{
	const value_type
		sTy          = sTy_in ? *sTy_in : AbstractLinAlgPack::dot(s,y),
		nrm_s        = s.norm_2(),
		nrm_y        = y.norm_2(),
		sqrt_macheps = ::sqrt(std::numeric_limits<value_type>::epsilon()),
		min_sTy      = sqrt_macheps * nrm_s * nrm_y;
	// Skip update if: s'*y < sqrt(macheps)*||s||2*||y||2 (Dennis and Schnabel, A9.4.2)
	const bool
		sufficiently_p_d = sTy > min_sTy;
	if( !sufficiently_p_d && out ) {
		if(func_name)
			*out << func_name << " : ";
		*out
			<< "Error, s'*y = " << sTy << " < sqrt(mach_eps) * ||s||2 * ||y||2 = "
			<< sqrt_macheps << " * " << nrm_s << " * " << nrm_y << " = " << min_sTy
			<< "\nTherefore the BFGS update is illdefined!\n";
	}
	return sufficiently_p_d;
}
