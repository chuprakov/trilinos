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

#ifndef FEAPP_ABSTRACTNODEBCSTRATEGY_HPP
#define FEAPP_ABSTRACTNODEBCSTRATEGY_HPP

#include <vector>

#include "FEApp_AbstractNodeBCStrategy_NTBase.hpp"
#include "FEApp_TemplateTypes.hpp"

#include "Sacado_TemplateManager.hpp"
#include "Sacado_mpl_placeholders.hpp"

using namespace Sacado::mpl::placeholders;

namespace FEApp {

  /*!
   * \brief Abstract interface for a nodal boundary condition
   */
  template <typename EvalT>
  class AbstractNodeBCStrategy : public FEApp::AbstractNodeBCStrategy_NTBase {
  public:

    //! Scalar type
    typedef typename FEApp::EvaluationTraits::apply<EvalT>::type ScalarT;

    //! Constructor
    AbstractNodeBCStrategy() {};

    //! Destructor
    virtual ~AbstractNodeBCStrategy() {};

    //! Evaluate BC residual
    virtual void evaluateResidual(const std::vector<ScalarT>* dot,
                                  const std::vector<ScalarT>& solution,
                                  std::vector<ScalarT>& residual) const = 0;

  private:

    //! Private to prohibit copying
    AbstractNodeBCStrategy(const AbstractNodeBCStrategy&);

    //! Private to prohibit copying
    AbstractNodeBCStrategy& operator=(const AbstractNodeBCStrategy&);

  };

  template <typename TypeSeq>
  class AbstractNodeBCStrategy_TemplateManager : 
    public Sacado::TemplateManager<TypeSeq, 
                                   FEApp::AbstractNodeBCStrategy_NTBase,
                                   AbstractNodeBCStrategy<_> > {
  public:
    AbstractNodeBCStrategy_TemplateManager() {}
    ~AbstractNodeBCStrategy_TemplateManager() {}
  };

}

#endif // FEAPP_ABSTRACTNODEBC_HPP
