// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

template <typename ScalarT>
FEApp::ConstantNodeBCStrategy<ScalarT>::
ConstantNodeBCStrategy(unsigned int solution_index, 
		       unsigned int residual_index,
		       const ScalarT& value) :
  sol_index(solution_index),
  res_index(residual_index),
  val(value),
  offsets(1)
{
  offsets[0] = residual_index;
}

template <typename ScalarT>
FEApp::ConstantNodeBCStrategy<ScalarT>::
~ConstantNodeBCStrategy()
{
}

template <typename ScalarT>
const std::vector<unsigned int>&
FEApp::ConstantNodeBCStrategy<ScalarT>::
getOffsets() const
{
  return offsets;
}

template <typename ScalarT>
void
FEApp::ConstantNodeBCStrategy<ScalarT>::
evaluateResidual(const std::vector<ScalarT>* dot,
		 const std::vector<ScalarT>& solution,
		 std::vector<ScalarT>& residual) const
{
  residual[res_index] = solution[sol_index] - val;
}

template <typename ScalarT>
void
FEApp::ConstantNodeBCStrategy<ScalarT>::
setValue(const ScalarT& value, bool mark_constant)
{
  val = value;
  if (mark_constant)
    Sacado::MarkConstant<ScalarT>::eval(val);
}

template <typename ScalarT>
const ScalarT&
FEApp::ConstantNodeBCStrategy<ScalarT>::
getValue() const
{
  return val;
}
