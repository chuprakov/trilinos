// ////////////////////////////////////////////////////////////////////////////
// rSQPAlgoClientInterface.h
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#ifndef RSQP_ALGO_CLIENT_INTERFACE_H
#define RSQP_ALGO_CLIENT_INTERFACE_H

#include "rSQPSolverClientInterface.h"

namespace ReducedSpaceSQPPack {

///
/** Interface that smart clients use to set the configuration
  * object that defines the rSQP algorithm to be used to solve
  * the NLP.
  */
class rSQPAlgoClientInterface : public rSQPSolverClientInterface {
public:

	/** @name Public Types */
	//@{

	///
	typedef ReferenceCountingPack::ref_count_ptr<rSQPAlgo_Config>	config_ptr_t;

	//@}

	/** @name «std comp» members for config. */
	//@{

	///
	virtual void set_config(const config_ptr_t& config) = 0;
	///
	virtual config_ptr_t& get_config() = 0;
	///
	virtual const config_ptr_t& get_config() const = 0;
	///
	virtual rSQPAlgo_Config& config() = 0;
	///
	virtual const rSQPAlgo_Config& config() const = 0;

	//@}

	///
	/** Call that causes the algorithm to be configured.
	  *
	  * Causes the config object to configure the algorithm
	  * to be ready to solve an NLP or print the algorithm.
	  *
	  * May be called after the nlp, track and config objects
	  * are set.
	  *
	  * Must be  called before print_algorithm(...) is called.
	  */
	virtual void configure_algorithm(std::ostream* trase_out = 0) = 0;

	/// Print the configured algorithm
	virtual void print_algorithm(std::ostream& out) const = 0;

};	// end class rSQPAlgoClientInterface

}	// end namespace ReducedSpaceSQPPack

#endif	// RSQP_ALGO_CLIENT_INTERFACE_H
