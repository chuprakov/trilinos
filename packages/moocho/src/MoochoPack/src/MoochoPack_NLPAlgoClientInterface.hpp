// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_NLPAlgoClientInterface.hpp
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

#include "MoochoPack_NLPSolverClientInterface.hpp"

namespace MoochoPack {

///
/** Interface that smart clients use to set the algorithm configuration
 * object that defines the rSQP algorithm to be used to solve the NLP.
 *
 * ToDo: Finish documentation!
 */
class NLPAlgoClientInterface : public NLPSolverClientInterface {
public:

	/** @name Public Types */
	//@{

	///
	typedef Teuchos::RefCountPtr<NLPAlgoConfig>	config_ptr_t;

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
	virtual NLPAlgoConfig& config() = 0;
	///
	virtual const NLPAlgoConfig& config() const = 0;

	//@}
	
	///
	/** Causes the algorithm to be configured.
	 *
	 * Causes the \c config object to configure the algorithm
	 * to be ready to solve an NLP or print the algorithm.
	 *
	 * May be called after the \c nlp, \c track and \c config objects
	 * are set.
	 *
	 * Must be  called before \c print_algorithm() or \c find_min() are called.
	  */
	virtual void configure_algorithm(std::ostream* trase_out = 0) = 0;

	/// Print the configured algorithm
	virtual void print_algorithm(std::ostream& out) const = 0;

private:

#ifdef DOXYGEN_COMPILE // Strictly for doxygen diagrams
	///
	NLPAlgoConfig    *config;
#endif

};	// end class NLPAlgoClientInterface

}	// end namespace MoochoPack

#endif	// RSQP_ALGO_CLIENT_INTERFACE_H
