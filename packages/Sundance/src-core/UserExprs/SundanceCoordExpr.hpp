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

#ifndef SUNDANCE_COORDEXPR_H
#define SUNDANCE_COORDEXPR_H

#include "SundanceFuncElementBase.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "SundanceCoordExprEvaluator.hpp"

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace SundanceCore::Internal;

  /** 
   * CoordExpr is an 
   * expression that returns one of the cartesian coordinates for
   * each point at which it evaluated. Which coordinate (i.e., \f$x\f$,
   * \f$y\f$, or \f$z\f$) to be represented is specified by the first
   * argement to the constructor. 
   */
  class CoordExpr : public FuncElementBase,
                    public EvaluatableExpr,
                    public GenericEvaluatorFactory<CoordExpr, CoordExprEvaluator>
    {
    public:
      /** */
      CoordExpr(int dir, const string& name="");

      /** */
      virtual ~CoordExpr() {;}

      /** */
      virtual XMLObject toXML() const ;

#ifndef DOXYGEN_DEVELOPER_ONLY
      /** */
      int dir() const {return dir_;}


    
      /** */
      virtual Set<MultipleDeriv> 
      internalFindW(int order, const EvalContext& context) const ;
      
      /** */
      virtual Set<MultipleDeriv> 
      internalFindC(int order, const EvalContext& context) const ;
      
      /** */
      virtual Set<MultipleDeriv> 
      internalFindV(int order, const EvalContext& context) const ;
      
      /** */
      virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}

      /** Ordering operator for use in transforming exprs to standard form */
      virtual bool lessThan(const ScalarExpr* other) const ;

      static string coordName(int dir, const string& name);

    private:
      int dir_;
#endif  /* DOXYGEN_DEVELOPER_ONLY */

    };
}

#endif
