/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009
 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TRel3DSizeBarrier.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TRel3DSizeBarrier.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

TRel3DSizeBarrier::~TRel3DSizeBarrier() {}

std::string TRel3DSizeBarrier::get_name() const
  { return "SizeBarrier"; }

bool TRel3DSizeBarrier::evaluate( const MsqMatrix<3,3>& T, 
                                    double& result, 
                                    MsqError& )
{
  double d = det(T);
  if (invalid_determinant(d)) {
    result = 0.0;
    return false;
  }
  result = d + 1.0/d - 2.0;
  return true;  
}

bool TRel3DSizeBarrier::evaluate_with_grad( const MsqMatrix<3,3>& T,
                                              double& result,
                                              MsqMatrix<3,3>& deriv_wrt_T,
                                              MsqError& err )
{
  double d = det(T);
  if (invalid_determinant(d)) {
    result = 0.0;
    return false;
  }
  result = d + 1.0/d - 2.0;
  deriv_wrt_T = (1 - 1/(d*d)) * transpose_adj(T);
  return true;  
}
  
bool TRel3DSizeBarrier::evaluate_with_hess( const MsqMatrix<3,3>& T,
                                              double& result,
                                              MsqMatrix<3,3>& deriv_wrt_T,
                                              MsqMatrix<3,3> second_wrt_T[6],
                                              MsqError& err )
{
  double d = det(T);
  if (invalid_determinant(d)) {
    result = 0.0;
    return false;
  }
  result = d + 1.0/d - 2.0;
  MsqMatrix<3,3> adjt = transpose_adj(T);
  const double f = 1 - 1/(d*d);
  deriv_wrt_T = f * adjt;
  
  set_scaled_outer_product( second_wrt_T, 2/(d*d*d), adjt );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_T, f, T );
  return true;  
}

} // namespace Mesquite
