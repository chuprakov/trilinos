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
	virtual void config_algo_cntr(rSQPAlgoContainer& algo_cntr) = 0;

	/// Deconfigure the rSQP algorithm container object.
	virtual void deconfig_algo_cntr(rSQPAlgoContainer& algo_cntr) = 0;

	/// Initialize the rSQP algorithm object for the start of the rSQP iterations
	virtual void init_algo(rSQPAlgoInterface& algo) = 0;

	///
	/** Prints a description of the algorithm.
	  *
	  */
	virtual void print_algorithm(std::ostream& out) const = 0;

	///
	virtual void print_state() const = 0;

};	// end class rSQPAlgo_Config

}	// end namespace ReducedSpaceSQPPack

#endif	// RSQP_ALGO_CONFIG_H