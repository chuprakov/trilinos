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

#ifndef FEAPP_CONSTANTNODEBCSTRATEGY_HPP
#define FEAPP_CONSTANTNODEBCSTRATEGY_HPP

#include "FEApp_AbstractNodeBCStrategy.hpp"

namespace FEApp {

  template <typename ScalarT>
  class ConstantNodeBCStrategy : 
    public FEApp::AbstractNodeBCStrategy<ScalarT> {
  public:

    //! Constructor
    ConstantNodeBCStrategy(unsigned int solution_index, 
			   unsigned int residual_index,
			   const ScalarT& value);

    //! Destructor
    virtual ~ConstantNodeBCStrategy();

    //! Get residual offsets
    const std::vector<unsigned int>& getOffsets() const;

    //! Evaluate BC residual
    virtual void evaluateResidual(const std::vector<ScalarT>* dot,
				  const std::vector<ScalarT>& solution,
				  std::vector<ScalarT>& residual) const;

    //! Set value of BC
    void setValue(const ScalarT& value, bool mark_constant);

    //! Get value of BC
    const ScalarT& getValue() const;

  private:
    
    //! Private to prohibit copying
    ConstantNodeBCStrategy(const ConstantNodeBCStrategy&);

    //! Private to prohibit copying
    ConstantNodeBCStrategy& operator=(const ConstantNodeBCStrategy&);

  protected:
    
    //! Index of solution variable
    unsigned int sol_index;

    //! Index of residual
    unsigned int res_index;

    //! Value of BC
    ScalarT val;

    //! Residual offsets
    std::vector<unsigned int> offsets;

  };

  class ConstantNodeBCStrategy_TemplateBuilder {
  public:
    ConstantNodeBCStrategy_TemplateBuilder(unsigned int solution_index, 
					   unsigned int residual_index,
					   double value) :
      sol_index(solution_index), res_index(residual_index), val(value) {}
    template <typename T>
    Teuchos::RefCountPtr<FEApp::AbstractNodeBCStrategy_NTBase> build() const {
      return Teuchos::rcp( new ConstantNodeBCStrategy<T>(sol_index, 
							 res_index, 
							 val));
    }
  protected:
    unsigned int sol_index;
    unsigned int res_index;
    double val;
  };

}

// Include implementation
#ifndef SACADO_ETI
#include "FEApp_ConstantNodeBCStrategyImpl.hpp"
#endif 

#endif // FEAPP_CONSTANTNODEBCSTRATEGY_HPP
