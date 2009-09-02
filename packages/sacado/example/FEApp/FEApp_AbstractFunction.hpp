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

#ifndef FEAPP_ABSTRACTFUNCTION_HPP
#define FEAPP_ABSTRACTFUNCTION_HPP

#include <vector>
#include "FEApp_TemplateTypes.hpp"

namespace FEApp {

  /*!
   * \brief Abstract interface for representing a PDE material function
   */
  template <typename EvalT>
  class AbstractFunction {
  public:

    //! Scalar type
    typedef typename FEApp::EvaluationTraits::apply<EvalT>::type ScalarT;
  
    //! Default constructor
    AbstractFunction() {};

    //! Destructor
    virtual ~AbstractFunction() {};

    //! Evaluate function at quadrature points
    virtual void
    evaluate(const std::vector<double>& quad_points,
             std::vector<ScalarT>& value) const = 0;

  private:
    
    //! Private to prohibit copying
    AbstractFunction(const AbstractFunction&);

    //! Private to prohibit copying
    AbstractFunction& operator=(const AbstractFunction&);

  };

}

#endif // FEAPP_ABSTRACTFUNCTION_HPP
