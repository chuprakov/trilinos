// ////////////////////////////////////////////////////////////////////////////
// rSQPAlgo.h

#ifndef RSQP_ALGO_H
#define RSQP_ALGO_H

#include "rSQPAlgoInterface.h"
#include "rSQPAlgoContainer.h"
#include "rSQPState.h"
#include "GeneralIterationPack/include/Algorithm.h"
#include "NLPInterfacePack/include/NLPReduced.h"

namespace ReducedSpaceSQPPack {

///
/** rSQP Algorithm control class.
  */
class rSQPAlgo : public rSQPAlgoInterface, public GeneralIterationPack::Algorithm {
public:

	/** @name Public Types */
	//@{

	//@}

	/// Constructs with no step, added_step, pre_step, post_step, state, or decomp_sys objects added.
	rSQPAlgo()
		: nlp_(0), algo_cntr_(0), trust_region_calc_(false), first_step_poss_(1)
	{}

	/** @name <<std aggr>> stereotype members for algo_cntr.
	  */
	//@{

	///
	void set_algo_cntr(rSQPAlgoContainer* algo_cntr);
	///
	rSQPAlgoContainer* get_algo_cntr()
	{	return algo_cntr_; }
	///
	rSQPAlgoContainer& algo_cntr()
	{	return StandardCompositionRelationshipsPack::role_name(algo_cntr_, false, "algo_cntr"); }
	///
	const rSQPAlgoContainer& algo_cntr() const
	{	return StandardCompositionRelationshipsPack::role_name(algo_cntr_, false, "algo_cntr"); }

	//@}

	/** @name <<std aggr>> stereotype members for nlp.
	  */
	//@{

	///
	void set_nlp(NLPReduced* nlp)
	{	nlp_ = nlp; }
	///
	NLPReduced* get_nlp()
	{	return nlp_; }
	///
	const NLPReduced* get_nlp() const
	{	return nlp_; }
	///
	NLPReduced& nlp()
	{	return StandardCompositionRelationshipsPack::role_name(nlp_, false, "nlp"); }
	///
	const NLPReduced& nlp() const
	{	return StandardCompositionRelationshipsPack::role_name(nlp_, false, "nlp"); }

	//@}

	// ////////////////////////////////////////////////////////////////////////////
	/** @name Trust region specialty.
	  *
	  * These operations apply to rSQP algorithms that use a trust region
	  * and have to peform trust region reduction minor loops from time to time.
	  * They are primarily an optimization that allow Stategy objects to be
	  * skipped if thier output variables or internal state will not change
	  * for a trust region reduction minor loop.
	  */
	//@{

	/// Set a flag for an active trust region reduction minor loop calculation.
	virtual void set_trust_region_calc(bool trust_region_calc)
	{	trust_region_calc_ = trust_region_calc; }
	
	/// See if a trust region reduction minor loop calculation is being performed.
	virtual bool is_trust_region_calc()
	{	return trust_region_calc_; }

	//@}

	///
	rSQPState& rsqp_state()
	{	return dynamic_cast<rSQPState&>(state()); }

	///
	const rSQPState& rsqp_state() const
	{	return dynamic_cast<const rSQPState&>(state()); }

	///
	void do_step_first(Algorithm::poss_type first_step_poss)
	{	first_step_poss_ = first_step_poss; }

	/** @name Overridden form rSQPAlgoInteface */
	//@{	
	
	///
	const rSQPState& retrieve_state() const;

	///
	void set_max_run_time(double max_run_time);
	///
	double return_max_run_time() const;

	///
	/** This is the main control function for the rSQP algorithm.  It is responsible
	  * for executing the main loop of the algorithm and does so sequentially until
	  * one of the Strategy objects that it calls to perform steps in the algorithm
	  * returns false indicating a request for a discontinous jump to another step.
	  * This is usually used to initiate a minor loop and #this->disbatch()# starts again
	  * from the step set by #this->do_step_next(stepid)#.
	  *
	  * This function returns true if the solution has been found (#this->found_min()# was called)
	  * or false if the maximum number of iterations has been reached
	  * (#this->state().k() > #this->algo_cntr().max_iter()#).
	  *
	  * This function starts executing with the step set by the last call
	  * to #this->do_step_next(stepid)#.  This allows the config object
	  * to start from the middle of a rSQP algorithm if it would like.  To start
	  * from the very beginning call #do_step_next(rSQPAlgo_Pack::FirstStep)#
	  * before you call #disbatch()#
	  *
	  */
	rSQPSolverClientInterface::EFindMinReturn dispatch();

	///
	void interface_print_algorithm(std::ostream& out) const;
	///
	void interface_set_algo_timing( bool algo_timing );
	///
	bool interface_algo_timing() const;
	///
	void interface_print_algorithm_times( std::ostream& out ) const;

	//@}

	/// overridden from Algorihth.

	///
	void print_algorithm(std::ostream& out) const;

protected:
	NLPReduced*				nlp_;			// reduced NLP being minimized
	rSQPAlgoContainer*		algo_cntr_;		// algorithm container object

	// Trust region reduction minor loop calculation flag
	bool trust_region_calc_;

	// First step to execute
	Algorithm::poss_type first_step_poss_;

};	// end class rSQPAlgo

}	// end namespace ReducedSpaceSQPPack

#endif	// RSQP_ALGO_H