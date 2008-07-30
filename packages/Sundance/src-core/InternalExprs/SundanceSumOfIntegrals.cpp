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

#include "SundanceSumOfIntegrals.hpp"
#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceSpectralExpr.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

SumOfIntegrals::SumOfIntegrals(const RefCountPtr<CellFilterStub>& region,
                               const Expr& expr,
                               const RefCountPtr<QuadratureFamilyStub>& quad)
  : ScalarExpr(), regions_(),
    quad_(),
    expr_(),
    cellSetToIndexMap_(),
    quadToIndexMap_()
{
  addTerm(region, expr, quad, 1);
}


Expr SumOfIntegrals::filterSpectral(const Expr& expr) const 
{
  const SpectralExpr* se = dynamic_cast<const SpectralExpr*>(expr.ptr().get());
  if (se != 0) return se->getCoeff(0);
  return expr;
}



void SumOfIntegrals::addTerm(const RefCountPtr<CellFilterStub>& regionPtr,
                             const Expr& expr,
                             const RefCountPtr<QuadratureFamilyStub>& quadPtr, 
                             int sign)
{
  int d = regions_.length();

  Expr ex = filterSpectral(expr);

  OrderedHandle<CellFilterStub> region(regionPtr);
  OrderedHandle<QuadratureFamilyStub> quad(quadPtr);

  if (cellSetToIndexMap_.containsKey(region))
    {
      d = cellSetToIndexMap_.get(region);
    }
  else
    {
      regions_.append(region);
      quadToIndexMap_.resize(d+1);
      quad_.resize(d+1);
      expr_.resize(d+1);
      cellSetToIndexMap_.put(region, d);
    }
  
  int q = quad_[d].length();
  if (quadToIndexMap_[d].containsKey(quad))
    {
      q = quadToIndexMap_[d].get(quad);
      if (sign > 0)
        {
          expr_[d][q] = expr_[d][q] + ex;
        }
      else
        {
          expr_[d][q] = expr_[d][q] - ex;
        }
    }
  else
    {
      quad_[d].append(quad);
      if (sign > 0)
        {
          expr_[d].append(ex);
        }
      else
        {
          expr_[d].append(-ex);
        }
      quadToIndexMap_[d].put(quad, q);
    }
}

void SumOfIntegrals::merge(const SumOfIntegrals* other, int sign) 
{
  for (unsigned int d=0; d<other->regions_.size(); d++)
    {
      for (int q=0; q<other->numTerms(d); q++)
        {
          addTerm(other->regions_[d].ptr(), other->expr_[d][q],
                  other->quad_[d][q].ptr(), sign);
        }
    }
}

void SumOfIntegrals::multiplyByConstant(const SpatiallyConstantExpr* expr) 
{
  double a = expr->value();

  for (unsigned int d=0; d<regions_.size(); d++)
    {
      for (int q=0; q<numTerms(d); q++)
        {
          expr_[d][q] = a*expr_[d][q];
        }
    }
}

void SumOfIntegrals::changeSign()
{
  for (unsigned int d=0; d<regions_.size(); d++)
    {
      for (int q=0; q<numTerms(d); q++)
        {
          expr_[d][q] = -expr_[d][q];
        }
    }
}

Set<int> SumOfIntegrals::funcsOnRegion(int d, const Set<int>& funcSet) const 
{
  Set<int> rtn;
  for (unsigned int t=0; t<expr_[d].size(); t++)
    { 
      expr_[d][t].ptr()->accumulateFuncSet(rtn, funcSet);
    }
  return rtn;
}


bool SumOfIntegrals::integralHasTestFunctions(int d) const 
{
  for (unsigned int t=0; t<expr_[d].size(); t++)
    { 
      if (expr_[d][t].ptr()->hasTestFunctions()) return true;
    }
  return false;
}



RefCountPtr<CellFilterStub> SumOfIntegrals::nullRegion() const
{
  for (unsigned int d=0; d<regions_.size(); d++)
    {
      if (!regions_[d].ptr()->isNullRegion())
        {
          return regions_[d].ptr()->makeNullRegion();
        }
    }
  TEST_FOR_EXCEPTION(true, RuntimeError,
                     "SumOfIntegrals::nullRegion() called on a sum "
                     "of integrals with no non-null regions");

  return RefCountPtr<CellFilterStub>();
}

bool SumOfIntegrals::isIndependentOf(const Expr& u) const
{
  for (unsigned int d=0; d<expr_.size(); d++)
  {
    for (unsigned int t=0; t<expr_[d].size(); t++)
    {
      if (!expr_[d][t].isIndependentOf(u)) return false;
    }
  }
  return true;
}

bool SumOfIntegrals::isLinearForm(const Expr& u) const
{
  for (unsigned int d=0; d<expr_.size(); d++)
  {
    for (unsigned int t=0; t<expr_[d].size(); t++)
    {
      if (!expr_[d][t].isLinearForm(u)) return false;
    }
  }
  return true;
}

bool SumOfIntegrals::isQuadraticForm(const Expr& u) const
{
  for (unsigned int d=0; d<expr_.size(); d++)
  {
    for (unsigned int t=0; t<expr_[d].size(); t++)
    {
      if (!expr_[d][t].isQuadraticForm(u)) return false;
    }
  }
  return true;
}


bool SumOfIntegrals::everyTermHasTestFunctions() const
{
  for (unsigned int d=0; d<expr_.size(); d++)
  {
    for (unsigned int t=0; t<expr_[d].size(); t++)
    {
      if (!expr_[d][t].everyTermHasTestFunctions()) return false;
    }
  }
  return true;
}


bool SumOfIntegrals::isLinearInTests() const
{
  for (unsigned int d=0; d<expr_.size(); d++)
  {
    for (unsigned int t=0; t<expr_[d].size(); t++)
    {
      if (!expr_[d][t].isLinearInTests()) return false;
    }
  }
  return true;
}

bool SumOfIntegrals::hasTestFunctions() const
{
  for (unsigned int d=0; d<expr_.size(); d++)
  {
    for (unsigned int t=0; t<expr_[d].size(); t++)
    {
      if (expr_[d][t].hasTestFunctions()) return true;
    }
  }
  return false;
}



ostream& SumOfIntegrals::toText(ostream& os, bool paren) const
{
  os << "Sum of Integrals[" << endl;
  for (unsigned int d=0; d<regions_.size(); d++)
    {
      for (unsigned int t=0; t<quad_[d].size(); t++)
        { 
          os << "Integral[" << endl;
          os << regions_[d].ptr()->toXML() << endl;
          os << "quad rule: " << quad_[d][t].ptr()->toXML() << endl;
          os << "expr: " << expr_[d][t].toString() << endl;
          os << "]" << endl;
        }

    }
  os << "]" << endl;

  return os;
}

ostream& SumOfIntegrals::toLatex(ostream& os, bool paren) const
{
  TEST_FOR_EXCEPTION(true, InternalError, 
                     "SumOfIntegrals::toLatex is undefined");
  return os;
}

XMLObject SumOfIntegrals::toXML() const 
{
  XMLObject rtn("SumOfIntegrals");
  for (unsigned int d=0; d<regions_.size(); d++)
    {
      XMLObject child("Integral");
      rtn.addChild(child);
      for (unsigned int t=0; t<quad_[d].size(); t++)
        {
          child.addChild(quad_[d][t].ptr()->toXML());
          child.addChild(expr_[d][t].toXML());
        }
    }

  return rtn;
}


bool SumOfIntegrals::lessThan(const ScalarExpr* other) const
{
  const SumOfIntegrals* f = dynamic_cast<const SumOfIntegrals*>(other);
  TEST_FOR_EXCEPTION(f==0, InternalError, "cast should never fail at this point");
  
  if (regions_.size() < f->regions_.size()) return true;
  if (regions_.size() > f->regions_.size()) return false;

  for (unsigned int i=0; i<regions_.size(); i++)
    {
      if (regions_[i] < f->regions_[i]) return true;
      if (f->regions_[i] < regions_[i]) return false;

      if (quad_[i] < f->quad_[i]) return true;
      if (f->quad_[i] < quad_[i]) return false;

      if (expr_[i].size() < f->expr_[i].size()) return true;
      if (expr_[i].size() > f->expr_[i].size()) return false;

      for (unsigned int j=0; j<expr_[i].size(); j++)
        {
          if (expr_[i][j].lessThan(f->expr_[i][j])) return true;
          if (f->expr_[i][j].lessThan(expr_[i][j])) return false;
        }
    }
  
  return false;
}
