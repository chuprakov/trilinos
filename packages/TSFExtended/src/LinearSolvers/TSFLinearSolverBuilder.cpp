/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/
/* @HEADER@ */

#include "TSFLinearSolverBuilder.hpp"
#include "TSFAztecSolver.hpp"
#include "TSFBICGSTABSolver.hpp"
#include "TSFGMRESSolver.hpp"
#include "Teuchos_XMLParameterListReader.hpp"

using namespace TSFExtended;
using namespace TSFExtendedOps;
using namespace Teuchos;

LinearSolver<double> LinearSolverBuilder::createSolver(const XMLObject& xml)
{
  string name = xml.getRequired("name");
  TEST_FOR_EXCEPTION(name != "Linear Solver", runtime_error,
                     "solver builder expected name [Linear Solver], got "
                     << name);

  string solverType;
  for (int i=0; i<xml.numChildren(); i++)
    {
      XMLObject child = xml.getChild(i);
      if (child.hasAttribute("Solver"))
        {
          solverType = child.getRequired("Solver");
          break;
        }
    }
  
  if (solverType=="Aztec")
    {
      XMLParameterListReader reader;
      ParameterList params = reader.toParameterList(xml);
      return new AztecSolver(params);
    }
  else if (solverType=="TSF")
    {
      XMLParameterListReader reader;
      ParameterList params = reader.toParameterList(xml);
      string method = params.template get<string>("Method");
      if (method=="BICGSTAB") 
        {
          return new BICGSTABSolver<double>(params);
        }
       else if (method=="GMRES")
         {
           return new GMRESSolver<double>(params);
         }
    }

  TEST_FOR_EXCEPTION(true, runtime_error, 
                     "Could not create a solver from XML object " 
                     << xml);
  return LinearSolver<double>();
    
}

LinearSolver<double> LinearSolverBuilder::createSolver(const ParameterList& params)
{
  TEST_FOR_EXCEPTION(!params.isSublist("Linear Solver"), runtime_error,
                     "did not find Linear Solver sublist in " << params);


  ParameterList solverSublist = params.sublist("Linear Solver");

  const string& solverType = getParameter<string>(solverSublist, "Type");

  if (solverType=="Aztec")
    {
      return new AztecSolver(solverSublist);
    }
  else if (solverType=="TSF")
    {
      const string& solverMethod = getParameter<string>(solverSublist, "Method");
      if (solverMethod=="BICGSTAB") 
        {
          return new BICGSTABSolver<double>(solverSublist);
        }
      else if (solverMethod=="GMRES")
        {
          return new GMRESSolver<double>(solverSublist);
        }
    }

  TEST_FOR_EXCEPTION(true, runtime_error, 
                     "Could not create a solver from parameter list " 
                     << params);
  return LinearSolver<double>();
    
}

