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
#include "Optika_GUI.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_VerboseObject.hpp"

/*
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!!!!!!!!!!!!!!!!!!     ATTENTION      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!!!!!!!!!!! You might want to look at the BasicExample first. !!!!!!!!!!
 * !!!!!!!!!!!!! It's not ciritical, but it has some important     !!!!!!!!!!
 * !!!!!!!!!!!!! information that you might find helpful.          !!!!!!!!!!
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 */

int main(int argc, char* argv[])
{
  /**
   * The end result of this example is exactly the same thing as we created
   * in the BasicExample. The difference is this time, instead of building
   * our GUI in source code, we're going to build it via XML. This approach
   * has a number of advantages.
   *    1. It's a lot cleaner. The source code approach can get pretty ugly
   *       especially when we get into some of the more complicated examples.
   *    2. Going along with the same vein of thought above, cleaner code is
   *       code that's easier to maintain.
   *    3. There's no need to recompile your GUI every time you change
   *       something. Since everything is loaded dynamically from the XML,
   *       there's no source code to recompile.
   */

  /**
   * The first thing we've gotta do is read in the ParameterList from the XML
   * so we can pass it to Optika. It's pretty easy to do.
   */
  Teuchos::RCP<Teuchos::ParameterList> myParameters;

  /**
   * Then we just call getInput! There's a little more to it, so let's
   * head on over to the inputs.xml file to see what's going on there.
   */
  Optika::getInput("inputs.xml", myParameters);

  return 0;
}

