#ifndef __Panzer_ResponseEvaluatorFactory_Functional_impl_hpp__
#define __Panzer_ResponseEvaluatorFactory_Functional_impl_hpp__

#include <string>

#include "Panzer_config.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_Integrator_Scalar.hpp"
#include "Panzer_ResponseScatterEvaluator_Functional.hpp"
#include "Panzer_Response_Functional.hpp"

namespace panzer {

template <typename EvalT>
Teuchos::RCP<ResponseBase> ResponseEvaluatorFactory_Functional<EvalT>::
buildResponseObject(const std::string & responseName) const
{ return Teuchos::rcp(new Response_Functional<EvalT>(responseName,comm_,linearObjFactory_)); }

template <typename EvalT>
void ResponseEvaluatorFactory_Functional<EvalT>::
buildAndRegisterEvaluators(const std::string & responseName,
                           PHX::FieldManager<panzer::Traits> & fm,
                           const panzer::PhysicsBlock & physicsBlock,
                           const Teuchos::ParameterList & user_data) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;


   // build integration evaluator (integrate over element)
   if(requiresCellIntegral_) {
     std::string field = (quadPointField_=="" ? responseName : quadPointField_);

     // build integration rule to use in cell integral
     RCP<IntegrationRule> ir = rcp(new IntegrationRule(cubatureDegree_,physicsBlock.cellData()));

     Teuchos::ParameterList pl;
     pl.set("Integral Name",field);
     pl.set("Integrand Name",field);
     pl.set("IR",ir);

     Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval 
         = Teuchos::rcp(new Integrator_Scalar<EvalT,panzer::Traits>(pl));
 
     fm.template registerEvaluator<EvalT>(eval);
   }

   // build scatter evaluator
   {
     std::string field = (quadPointField_=="" ? responseName : quadPointField_);

     // build useful evaluator
     Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval 
         = Teuchos::rcp(new ResponseScatterEvaluator_Functional<EvalT,panzer::Traits>(field,responseName,physicsBlock.cellData()));

     fm.template registerEvaluator<EvalT>(eval);

     // require last field
     fm.template requireField<EvalT>(*eval->evaluatedFields()[0]);
   }
}

template <typename EvalT>
bool ResponseEvaluatorFactory_Functional<EvalT>::
typeSupported() const
{
   return false;
}

template < >
bool ResponseEvaluatorFactory_Functional<panzer::Traits::Residual>::
typeSupported() const
{
  return true;
}

template < >
bool ResponseEvaluatorFactory_Functional<panzer::Traits::Jacobian>::
typeSupported() const
{
  return linearObjFactory_!=Teuchos::null;
}

#ifdef HAVE_STOKHOS
template < >
bool ResponseEvaluatorFactory_Functional<panzer::Traits::SGResidual>::
typeSupported() const
{
  return true;
}
#endif

}

#endif
