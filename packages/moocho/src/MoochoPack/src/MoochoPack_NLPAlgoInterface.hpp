// ////////////////////////////////////////////////////////////////////////////
// rSQPAlgoInterface.h

#ifndef RSQP_ALGO_INTERFACE_H
#define RSQP_ALGO_INTERFACE_H

#include "ReducedSpaceSQPPackTypes.h"
#include "rSQPSolverClientInterface.h"

namespace ReducedSpaceSQPPack {

///
/** Interface rSQPAlgoContainer uses to access rSQPAlgo.
  */
class rSQPAlgoInterface {
public:

	// Virtual destructor
	virtual ~rSQPAlgoInterface() {}

	/// Return the state object
	virtual const rSQPState& retrieve_state() const = 0;

	///
	/** Start the iterations.
	  *
	  * This function returns true if the solution was found and false
	  * if the maximum number of iterations was reached before the
	  * solution was found.
	  */
	virtual rSQPSolverClientInterface::EFindMinReturn dispatch() = 0;

	/** @name Set options on Algorithm */
	//@{

	///
	virtual void interface_print_algorithm(std::ostream& out) const = 0;
	///
	virtual void interface_set_algo_timing( bool algo_timing ) = 0;
	///
	virtual bool interface_algo_timing() const = 0;
	///
	virtual void interface_print_algorithm_times( std::ostream& out ) const = 0;

	//@}

};	// end class rSQPAlgoInterface

}	// end namespace ReducedSpaceSQPPack

#endif	// RSQP_ALGO_INTERFACE_H