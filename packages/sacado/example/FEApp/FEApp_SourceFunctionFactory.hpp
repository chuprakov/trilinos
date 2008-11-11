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

#ifndef FEAPP_SOURCEFUNCTIONFACTORY_HPP
#define FEAPP_SOURCEFUNCTIONFACTORY_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "FEApp_AbstractSourceFunction.hpp"

#include "Sacado_ScalarParameterLibrary.hpp"

namespace FEApp {

  /*!
   * \brief A factory class to instantiate AbstractSourceFunction objects
   */
  template <typename EvalT>
  class SourceFunctionFactory {
  public:

    //! Default constructor
    SourceFunctionFactory(
	       const Teuchos::RCP<Teuchos::ParameterList>& funcParams,
	       const Teuchos::RCP<ParamLib>& paramLib);

    //! Destructor
    virtual ~SourceFunctionFactory() {}

    virtual Teuchos::RCP< FEApp::AbstractSourceFunction<EvalT> >
    create();

  private:

    //! Private to prohibit copying
    SourceFunctionFactory(const SourceFunctionFactory&);

    //! Private to prohibit copying
    SourceFunctionFactory& operator=(const SourceFunctionFactory&);

  protected:

    //! Parameter list specifying what strategy to create
    Teuchos::RCP<Teuchos::ParameterList> funcParams;

    //! Parameter library
    Teuchos::RCP<ParamLib> paramLib;

  };

}

// Include implementation
#ifndef SACADO_ETI
#include "FEApp_SourceFunctionFactoryImpl.hpp"
#endif 

#endif // FEAPP_SOURCEFUNCTIONFACTORY_HPP
