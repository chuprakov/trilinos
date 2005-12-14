// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Moocho_ConfigDefs.hpp"
#include "check_nan_inf.h"

// Computed values to compare to
RTOp_value_type
    RTOp_pos_inf = +1.0/sin(0.0),
    RTOp_neg_inf = -1.0/sin(0.0),
    RTOp_pos_nan = +0.0/sin(0.0),
    RTOp_neg_nan = -0.0/sin(0.0);

int RTOp_is_nan( RTOp_value_type val )
{
#if defined(_INTEL_CXX)
  return _isnan(val) != 0;
#else
  return val == RTOp_pos_nan || val == RTOp_neg_nan || val != val;
#endif
}

int RTOp_is_inf( RTOp_value_type val )
{
    return val == RTOp_pos_inf || val == RTOp_neg_inf; // IEEE math
}

int RTOp_is_nan_inf( RTOp_value_type val )
{
  return
#if defined(_INTEL_CXX)
    _isnan(val) != 0
#else
	  val == RTOp_pos_nan || val == RTOp_neg_nan || val != val
#endif
	  || val == RTOp_pos_inf || val == RTOp_neg_inf;
}
