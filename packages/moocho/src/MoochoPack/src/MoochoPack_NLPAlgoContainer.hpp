// ////////////////////////////////////////////////////////////////////////////
// rSQPAlgoContainer.h

#ifndef RSQP_ALGO_CONTAINER_H
#define RSQP_ALGO_CONTAINER_H

#include "rSQPAlgoClientInterface.h"
#include "rSQPAlgoInterface.h"
#include "rSQPAlgo_Config.h"

namespace ReducedSpaceSQPPack {

///
/** Implementation for rSQPAlgo solver.
  *
  * Acts as a container for rSQPAlgo.  This class is hidden from clients
  * by not exposing it to them in header files.
  */
class rSQPAlgoContainer : public rSQPAlgoClientInterface {
public:

	/// Construct a container with no configuration object set.
	rSQPAlgoContainer()
	{}

	/** @name <<std aggr>> stereotype members for algo.
	  */
	//@{

	///
	void set_algo(rSQPAlgoInterface* algo)
	{	algo_ = algo; }
	///
	rSQPAlgoInterface* get_algo()
	{	return algo_; }
	///
	const rSQPAlgoInterface* get_algo() const
	{	return algo_; }
	///
	rSQPAlgoInterface& algo()
	{	return StandardCompositionRelationshipsPack::role_name(algo_, false, "algo"); }
	///
	const rSQPAlgoInterface& algo() const
	{	return StandardCompositionRelationshipsPack::role_name(algo_, false, "algo"); }

	//@}

	/** @name Overridden.*/
	//@{

	///
	void set_max_run_time(double max_run_time);
	///
	double max_run_time() const;

	/** @name «std comp» members for config. */
	//@{

	///
	void set_config(const config_ptr_t& config);
	///
	config_ptr_t& get_config()
	{	return config_ ; }
	///
	const config_ptr_t& get_config() const
	{	return config_; }
	///
	rSQPAlgo_Config& config()
	{	return *config_; }
	///
	const rSQPAlgo_Config& config() const
	{	return *config_; }

	//@}

	///
	EFindMinReturn find_min();
	///	
	const rSQPState& state() const;
	///
	void configure_algorithm();
	///
	void print_algorithm(std::ostream& out) const;
	///
	void set_algo_timing( bool algo_timing );
	///
	bool algo_timing() const;
	///
	void print_algorithm_times( std::ostream& out ) const;

	//@}

private:
	config_ptr_t			config_;
	rSQPAlgoInterface		*algo_;

	// Assert that the object has been set up properly and throw exception if it has not
	void assert_valid_setup() const;

	// Not defined and not to be called
	rSQPAlgoContainer(const rSQPAlgoContainer&);
	rSQPAlgoContainer& operator=(const rSQPAlgoContainer&);

};	// end class rSQPAlgoContainer

}	// end namespace ReducedSpaceSQPPack

#endif	// RSQP_ALGO_CONTAINER_H