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

#ifndef FEAPP_BRUSSELATORNODEBCSTRATEGY_HPP
#define FEAPP_BRUSSELATORNODEBCSTRATEGY_HPP

#include "Teuchos_RCP.hpp"

#include "FEApp_AbstractNodeBCStrategy.hpp"

#include "Sacado_ScalarParameterLibrary.hpp"

namespace FEApp {

  template <typename EvalT>
  class BrusselatorNodeBCStrategy : 
    public FEApp::AbstractNodeBCStrategy<EvalT> {
  public:

    //! Scalar type
    typedef typename FEApp::AbstractNodeBCStrategy<EvalT>::ScalarT ScalarT;

    //! Constructor
    BrusselatorNodeBCStrategy(
                const ScalarT& alpha, const ScalarT& beta,
                const Teuchos::RCP<ParamLib>& paramLib);

    //! Destructor
    virtual ~BrusselatorNodeBCStrategy();

    //! Get residual offsets
    const std::vector<unsigned int>& getOffsets() const;

    //! Evaluate BC residual
    virtual void evaluateResidual(const std::vector<ScalarT>* dot,
                                  const std::vector<ScalarT>& solution,
                                  std::vector<ScalarT>& residual) const;

  private:
    
    //! Private to prohibit copying
    BrusselatorNodeBCStrategy(const BrusselatorNodeBCStrategy&);

    //! Private to prohibit copying
    BrusselatorNodeBCStrategy& operator=(const BrusselatorNodeBCStrategy&);

  protected:

    //! Value of alpha
    mutable ScalarT alpha;

    //! Value of beta
    mutable ScalarT beta;

    //! Residual offsets
    std::vector<unsigned int> offsets;

    //! Parameter library
    Teuchos::RCP<ParamLib> pl;

  };

  class BrusselatorNodeBCStrategy_TemplateBuilder {
  public:
    BrusselatorNodeBCStrategy_TemplateBuilder(
        double alpha, double beta,
        const Teuchos::RCP<ParamLib>& paramLib) :
      a(alpha), b(beta), pl(paramLib) {}
    template <typename T>
    Teuchos::RCP<FEApp::AbstractNodeBCStrategy_NTBase> build() const {
      return Teuchos::rcp( new BrusselatorNodeBCStrategy<T>(a,b,pl));
    }
  protected:
    double a, b;
    Teuchos::RCP<ParamLib> pl;
  };

}

// Include implementation
#ifndef SACADO_ETI
#include "FEApp_BrusselatorNodeBCStrategyImpl.hpp"
#endif 

#endif // FEAPP_BRUSSELATORNODEBCSTRATEGY_HPP
