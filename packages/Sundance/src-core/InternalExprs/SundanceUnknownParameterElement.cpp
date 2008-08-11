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

#include "SundanceUnknownParameterElement.hpp"

#include "SundanceFunctionalDeriv.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceCoordDeriv.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

UnknownParameterElement
::UnknownParameterElement(const string& name,
  const string& suffix,
  int commonFuncID,
  int myIndex)
	: UnknownFuncElement(rcp(new UnknownFuncDataStub()), name, 
    suffix, commonFuncID, myIndex),
    SpatiallyConstantExpr()
{}



Set<MultipleDeriv> 
UnknownParameterElement::internalFindW(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;

  if (order==0) 
    {
      if (!evalPtIsZero()) rtn.put(MultipleDeriv());
    }
  else if (order==1)
    {
      Deriv d = new FunctionalDeriv(this, MultiIndex());
      MultipleDeriv md;
      md.put(d);
      rtn.put(md);
    }

  return rtn;
}



Set<MultipleDeriv> 
UnknownParameterElement::internalFindV(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;

  return rtn;
}


Set<MultipleDeriv> 
UnknownParameterElement::internalFindC(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;

  if (order==1)
    {
      Deriv d = new FunctionalDeriv(this, MultiIndex());
      MultipleDeriv md;
      md.put(d);
      rtn.put(md);
    }
  return rtn.intersection(UnknownFuncElement::findR(order, context));
}



Evaluator* UnknownParameterElement
::createEvaluator(const EvaluatableExpr* expr,
                  const EvalContext& context) const 
{
  return SymbolicFuncElement::createEvaluator(expr, context);
}


const Parameter* UnknownParameterElement::parameterValue() const 
{
  const Parameter* p = dynamic_cast<const Parameter*>(evalPt());
  TEST_FOR_EXCEPTION(p==0, InternalError, 
                     "UnknownParameter evalPt() is not a Parameter");
  return p;
}

Parameter* UnknownParameterElement::parameterValue()  
{
  Parameter* p = dynamic_cast<Parameter*>(evalPt());
  TEST_FOR_EXCEPTION(p==0, InternalError, 
                     "UnknownParameter evalPt() is not a Parameter");
  return p;
}






XMLObject UnknownParameterElement::toXML() const 
{
	XMLObject rtn("UnknownParameterElement");
	rtn.addAttribute("name", name());
	return rtn;
}

