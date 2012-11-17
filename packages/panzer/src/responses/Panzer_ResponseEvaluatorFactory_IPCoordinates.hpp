#ifndef __Panzer_ResponseEvaluatorFactory_IPCoordinates_hpp__
#define __Panzer_ResponseEvaluatorFactory_IPCoordinates_hpp__

#include <string>

#include "Panzer_config.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ResponseEvaluatorFactory.hpp"

#include <mpi.h>

namespace panzer {

/** This class defines a response based on a functional.
  */
template <typename EvalT> 
class ResponseEvaluatorFactory_IPCoordinates : public ResponseEvaluatorFactory<EvalT> {
public:

   ResponseEvaluatorFactory_IPCoordinates(int cubatureDegree)
     : cubatureDegree_(cubatureDegree) {}

   virtual ~ResponseEvaluatorFactory_IPCoordinates() {}
 
   /** Build the response object used by this factory. This object
     * assumes the role of the scatter target and will be accessible
     * by all the evaluators in the field managers. 
     *
     * \param[in] responseName Name of response to be built. This
     *                         name will be used for looking up
     *                         the response in the <code>GlobalEvaluationDataContainer</code>
     *                         object.
     */
   virtual Teuchos::RCP<ResponseBase> buildResponseObject(const std::string & responseName,
                                                          const std::vector<std::string> & eBlocks) const;
   virtual Teuchos::RCP<ResponseBase> buildResponseObject(const std::string & responseName,
                                                          const std::vector<std::pair<std::string,std::string> > & sideset_blocks) const
   { TEUCHOS_ASSERT(false); return Teuchos::null; }

   /** Build and register evaluators for a response on a particular physics
     * block. 
     *
     * \param[in] responseName The name of the response to be constructed
     *                         by these evaluators.
     * \param[in,out] fm Field manager to be fuild with the evaluators.
     * \param[in] physicsBlock What physics block is being used for constructing
     *                         the evaluators
     * \param[in] user_data The user data parameter list, this stores things
     *                      that the user may find useful.
     */
   virtual void buildAndRegisterEvaluators(const std::string & responseName,
                                           PHX::FieldManager<panzer::Traits> & fm,
                                           const panzer::PhysicsBlock & physicsBlock,
                                           const Teuchos::ParameterList & user_data) const;

protected:
   //! Accessor method for Cubature degree (can be used by sub classes)
   int getCubatureDegree() const { return cubatureDegree_; }

private:
   int cubatureDegree_;
};

}

#include "Panzer_ResponseEvaluatorFactory_IPCoordinates_impl.hpp"

#endif
