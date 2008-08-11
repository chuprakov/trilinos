/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
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
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "SundanceBasisFamily.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceDiscreteFuncElement.hpp"

using namespace TSFExtended;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;





int BasisFamily::order() const 
{
  return ptr()->order();
}

int BasisFamily::dim() const 
{
  return ptr()->dim();
}

bool BasisFamily::operator==(const BasisFamily& other) const
{
  return !(*this < other || other < *this);
}

unsigned int BasisFamily::size(const Array<BasisFamily>& b)
{
  unsigned int rtn = 0;
  for (unsigned int i=0; i<b.size(); i++) rtn += b[i].dim();
  return rtn;
}

int BasisFamily::nReferenceDOFs(const CellType& maximalCellType,
  const CellType& cellType) const 
{
  return ptr()->nReferenceDOFs(maximalCellType, cellType);
}

BasisFamily BasisFamily::getBasis(const Expr& expr)
{
  TEST_FOR_EXCEPTION(expr.size() > 1, RuntimeError, "non-scalar expression in BasisFamily::getBasis()");

  /* First try to process assuming the input is an unknown function */
  const UnknownFuncElement* u 
    = dynamic_cast<const UnknownFuncElement*>(expr[0].ptr().get());
  if (u != 0)
    {
      const UnknownFunctionData* data = UnknownFunctionData::getData(u);
      return data->basis()[u->myIndex()];
    }

  /* Next try to process assuming the input is a test function */
  const TestFuncElement* t 
    = dynamic_cast<const TestFuncElement*>(expr[0].ptr().get());
  if (t != 0)
    {
      const TestFunctionData* data = TestFunctionData::getData(t);
      return data->basis()[t->myIndex()];
    }

  /* Next try to process assuming the input is a discrete function */
  const DiscreteFuncElement* d
    = dynamic_cast<const DiscreteFuncElement*>(expr[0].ptr().get());
  if (d != 0)
    {
      const DiscreteFunctionData* data = DiscreteFunctionData::getData(d);
      return data->discreteSpace().basis()[d->myIndex()];
    }

  TEST_FOR_EXCEPTION(u==0 && t==0 && d==0, RuntimeError,
                     "BasisFamily::getBasis() argument is not a function");
  return BasisFamily();
  
}




namespace SundanceStdFwk
{

Array<int> vectorDimStructure(const Array<BasisFamily>& basis)
{
  Array<int> rtn;
  for (unsigned int i=0; i<basis.size(); i++) rtn.append(basis[i].dim());
  return rtn;
}


Array<int> vectorDimStructure(const BasisFamily& basis)
{
  return vectorDimStructure(tuple(basis));
}

}
