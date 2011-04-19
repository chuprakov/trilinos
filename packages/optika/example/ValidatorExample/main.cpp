// @HEADER
// ***********************************************************************
// 
//         Optika: A Tool For Developing Parameter Obtaining GUIs
//                Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, with Sandia Corporation, the 
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Kurtis Nusbaum (klnusbaum@gmail.com) 
// 
// ***********************************************************************
// @HEADER
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Optika_GUI.hpp"
int main(){
 /*
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  * !!!!!!!!!!!!!!!!!!              ATTENTION              !!!!!!!!!!!!!!!!!!!!
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  * !!!!  PLEASE VIEW THE BASIC XML EXAMPLE FIRST BEFORE READING THIS    !!!!!!
  * !!!!  EXAMPLE. IT PROVIDES FUNDAMENTAL  KNOWLEDGE THAT WILL BE VERY  !!!!!!
  * !!!!  HELPFUL IN UNDERSTANDING THIS EXAMPLE.                         !!!!!!
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  */ 
  
  /**
   * Since validators can get unruly, we're going to do this example in XML.
   * so let's head over the the inputs.xml file for the rest of this example.
   */
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::rcp;
  RCP<ParameterList> myValidParameters = rcp(new ParameterList);
  Optika::getInput("inputs.xml", myValidParameters);
  return 0;
}

