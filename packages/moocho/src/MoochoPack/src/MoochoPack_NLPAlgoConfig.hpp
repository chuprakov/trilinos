// ////////////////////////////////////////////////////////////////////////////
// rSQPAlgo_Config.h

#ifndef RSQP_ALGO_CONFIG_H
#define RSQP_ALGO_CONFIG_H

#include "ReducedSpaceSQPPackTypes.h"

namespace ReducedSpaceSQPPack {

class rSQPAlgoContainer;
class rSQPAlgoInterface;

///
/** Responsible for configuring the containter with an rSQP algorithm object, configuring
  * the algo with step, state etc. objects and initailizing the algorithm before the
  * interations start.
  */
class rSQPAlgo_Config {
public:

	/** @name exceptions */
	//@{

	/// Thrown if NLP type is incompatible with this config
	class InvalidNLPType : public std::logic_error
	{public: InvalidNLPType(const std::string& what_arg) : std::logic_error(what_arg) {}};

	//@}

	// Virtual destructor
	virtual ~rSQPAlgo_Config() {}

	/// Configure the rSQP algorithm container with an rSQP algorithm object (non-templated).
	virtual void config_algo_cntr(rSQPAlgoContainer& algo_cntr, std::ostream* trase_out = 0) = 0;

	/// Initialize the rSQP algorithm object for the start of the rSQP iterations
	virtual void init_algo(rSQPAlgoInterface& algo) = 0;

};	// end class rSQPAlgo_Config

}	// end namespace ReducedSpaceSQPPack

#endif	// RSQP_ALGO_CONFIG_H