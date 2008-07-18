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

#ifndef SUNDANCE_GROUPER_H
#define SUNDANCE_GROUPER_H

#include "SundanceDefs.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceIntegralGroup.hpp"
#include "SundanceEquationSet.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;

namespace Internal
{
using namespace Teuchos;
using TSFExtended::ObjectWithVerbosity;  

/** 
 * Grouper
 */
class GrouperBase
  : public TSFExtended::ParameterControlledObjectWithVerbosity<GrouperBase>
{
public:
  /** */
  GrouperBase() {}

  /** */
  GrouperBase(const ParameterList& verbParams)
    : ParameterControlledObjectWithVerbosity<GrouperBase>("Grouper", verbParams)
    {}
             /** */
    virtual ~GrouperBase(){;}

  /** */
  virtual void findGroups(const EquationSet& eqn,
    const CellType& maxCellType,
    int spatialDim,
    const CellType& cellType,
    int cellDim,
    const QuadratureFamily& quad,
    const RefCountPtr<SparsitySuperset>& sparsity,
    Array<IntegralGroup>& groups) const = 0 ;

      /** */
      static RefCountPtr<ParameterList> defaultVerbParams()
        {
          static RefCountPtr<ParameterList> rtn = rcp(new ParameterList("Grouper"));
          static int first = true;
          if (first)
          {
            rtn->set<int>("global", 0);
            first = false;
          }
          return rtn;
        }

protected:
  void extractWeakForm(const EquationSet& eqn,
    const MultipleDeriv& functionalDeriv,
    BasisFamily& testBasis, 
    BasisFamily& unkBasis,
    MultiIndex& miTest, MultiIndex& miUnk,
    int& rawVarID, int& rawUnkID,  
    int& reducedTestID, int& reducedUnkID, 
    int& testBlock, int& unkBlock, 
    bool& isOneForm) const ;
                              
};

}
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
