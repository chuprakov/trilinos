// ///////////////////////////////////////////////////////////////
// Algorithm.cpp
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

#include <assert.h>
#include <signal.h>

#include <iterator>
#include <numeric>
#include <typeinfo>

#include "IterationPack/src/Algorithm.hpp"
#include "stpwatch.hpp"
#include "Teuchos_TestForException.hpp"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

// Define to see MPI/interrupt deugging output
//#define ITERATION_PACK_ALGORITHM_SHOW_MPI_DEBUG_INFO

// Define of the MPI implementation receives signals on all processes
//#define ITERATION_PACK_ALGORITHM_SIGNALS_ON_ALL_PROCESSES;

extern "C" {

void sig_handler_interrupt_algorithm( int signum )
{
	IterationPack::Algorithm::interrupt();
}

} // extern "C"

namespace {

template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }

} // end namespace

// ToDo: change step_itr and assoc_step_itr to just return iterators without
// asserting if the names exist.  This will be more useful.

namespace IterationPack {

// constructors / destructor

Algorithm::Algorithm()
	:running_state_(NOT_RUNNING), max_iter_(100)
	,max_run_time_(std::numeric_limits<double>::max())
	,next_step_name_(0), do_step_next_called_(false), reconfigured_(false)
	,time_stats_computed_(false)
{
	// Register the signal handler for the SIGINT
	signal( SIGINT, &sig_handler_interrupt_algorithm );
	interrupt_called_ = false;
	processed_user_interrupt_ = false;
	// Set MPI info
	num_proc_ = 1;
	proc_rank_ = 0;
#ifdef HAVE_MPI
	// ToDo: Allow the specification of another communicator if needed!
	MPI_Comm_size( MPI_COMM_WORLD, &num_proc_ );
	MPI_Comm_rank( MPI_COMM_WORLD, &proc_rank_ );
#ifdef ITERATION_PACK_ALGORITHM_SHOW_MPI_DEBUG_INFO
	std::cerr << "\np=" << proc_rank_ << ": Algorithm::Algorithm() being called (num_proc = "<<num_proc_<<") ... \n";
#endif
#endif // HAVE_MPI
}

Algorithm::~Algorithm()
{}

// maximum iterations

void Algorithm::max_iter(size_t max_iter)
{	max_iter_ = max_iter; }

size_t Algorithm::max_iter() const
{	return max_iter_; }

// maximum run tine

void Algorithm::max_run_time(double max_run_time)
{	max_run_time_ = max_run_time; }

double Algorithm::max_run_time() const
{	return max_run_time_; }


// step information / access

int Algorithm::num_steps() const
{	return steps_.size(); }

Algorithm::poss_type Algorithm::get_step_poss(const std::string& step_name) const
{	
	steps_t::const_iterator itr = step_itr(step_name);
	return itr == steps_.end() ? DOES_NOT_EXIST : std::distance( steps_.begin(), itr ) + 1;
}

const std::string& Algorithm::get_step_name(poss_type step_poss) const
{	return steps_[validate(step_poss) - 1].name; }

Algorithm::step_ptr_t& Algorithm::get_step(poss_type step_poss)
{	return steps_[validate(step_poss) - 1].step_ptr; }

const Algorithm::step_ptr_t& Algorithm::get_step(poss_type step_poss) const
{	return steps_[validate(step_poss) - 1].step_ptr; }

// pre/post step information / access

int Algorithm::num_assoc_steps(poss_type step_poss, EAssocStepType type) const
{	return assoc_steps_[validate(step_poss) - 1][type].size(); }

Algorithm::poss_type Algorithm::get_assoc_step_poss(poss_type step_poss, EAssocStepType type
	,const std::string& assoc_step_name) const
{	
	// ToDo: change to return DOES_NOT_EXIST if it does not exist.
	const assoc_steps_ele_list_t &assoc_list = assoc_steps_[validate(step_poss) - 1][type];
	assoc_steps_ele_list_t::const_iterator itr = assoc_step_itr(assoc_list,assoc_step_name);
	return itr == assoc_list.end() ? DOES_NOT_EXIST : std::distance( assoc_list.begin() , itr ) + 1;
}

const std::string& Algorithm::get_assoc_step_name(poss_type step_poss, EAssocStepType type
	, poss_type assoc_step_poss) const
{
	const assoc_steps_ele_list_t &assoc_list= assoc_steps_[validate(step_poss) - 1][type];
	validate(assoc_list,assoc_step_poss);
	assoc_steps_ele_list_t::const_iterator itr = assoc_list.begin();
	std::advance( itr, assoc_step_poss - 1 );
	return (*itr).name;
}

Algorithm::step_ptr_t& Algorithm::get_assoc_step(poss_type step_poss, EAssocStepType type
	, poss_type assoc_step_poss)
{
	assoc_steps_ele_list_t &assoc_list= assoc_steps_[validate(step_poss) - 1][type];
	validate(assoc_list,assoc_step_poss);
	assoc_steps_ele_list_t::iterator itr = assoc_list.begin();
	std::advance( itr, assoc_step_poss - 1 );
	return (*itr).step_ptr;
}

const Algorithm::step_ptr_t& Algorithm::get_assoc_step(poss_type step_poss, EAssocStepType type
	, poss_type assoc_step_poss) const
{
	const assoc_steps_ele_list_t &assoc_list= assoc_steps_[validate(step_poss) - 1][type];
	validate(assoc_list,assoc_step_poss);
	assoc_steps_ele_list_t::const_iterator itr = assoc_list.begin();
	std::advance( itr, assoc_step_poss - 1 );
	return (*itr).step_ptr;
}

// step manipulation

void Algorithm::insert_step(poss_type step_poss, const std::string& step_name, const step_ptr_t& step)
{
	validate_not_in_state(RUNNING);
	TEST_FOR_EXCEPTION(
		step.get() == NULL, std::invalid_argument
		,"Algorithm::insert_step(...) : A step with the name = \'" << step_name
		<< "\' being inserted into the position  = " << step_poss
		<< " has step.get() == NULL!" );
	// Make sure a step with this name does not already exist.
	steps_t::iterator itr;
	if( steps_.end() != ( itr = step_itr(step_name) ) )
		TEST_FOR_EXCEPTION(
			true, AlreadyExists
			,"Algorithm::insert_step(...) : A step with the name = " << step_name
			<< " already exists at step_poss = " << std::distance(steps_.begin(),itr) + 1 );
	// insert the step in such a way that any container can be used for steps_
	itr = steps_.begin();
	std::advance ( itr , validate(step_poss,+1) - 1 );
	steps_.insert( itr , steps_ele_t(step,step_name) );
	// insert the assoc_step element in such a way that any container can be used for assoc_steps_
	assoc_steps_t::iterator a_itr = assoc_steps_.begin();
	std::advance ( a_itr , step_poss - 1 );
	assoc_steps_.insert( a_itr , assoc_steps_ele_t() );
}

void Algorithm::change_step_name(poss_type step_poss, const std::string& new_name)
{
	validate_not_in_state(RUNNING);
	if(running_state() == RUNNING_BEING_CONFIGURED) {
		validate_not_curr_step(validate(step_poss));
		validate_not_next_step(steps_[step_poss - 1].name);
	}
	steps_[step_poss - 1].name = new_name;
}

void Algorithm::replace_step(poss_type step_poss, const step_ptr_t& step)
{
	validate_not_in_state(RUNNING);
	if(running_state() == RUNNING_BEING_CONFIGURED)	validate_not_curr_step(validate(step_poss));
	steps_[step_poss - 1].step_ptr = step;
}

void Algorithm::remove_step(poss_type step_poss)
{
	validate_not_in_state(RUNNING);
	if(running_state() == RUNNING_BEING_CONFIGURED) {
		validate_not_curr_step(validate(step_poss));
		validate_not_next_step(steps_[step_poss - 1].name);
	}
	// remove the step in such a way that any container can be used for steps_
	steps_t::iterator itr = steps_.begin();
	std::advance ( itr , validate(step_poss) - 1 );
	steps_.erase( itr );
	// remove the assoc_step element in such a way that any container can be used for assoc_steps_
	assoc_steps_t::iterator a_itr = assoc_steps_.begin();
	std::advance ( a_itr , step_poss - 1 );
	assoc_steps_.erase( a_itr );
}

// pre/post step manipulation

void Algorithm::insert_assoc_step(poss_type step_poss, EAssocStepType type, poss_type assoc_step_poss
	, const std::string& assoc_step_name, const step_ptr_t& assoc_step)
{
	validate_not_in_state(RUNNING);
	TEST_FOR_EXCEPTION(
		assoc_step.get() == NULL, std::invalid_argument
		,"Algorithm::insert_assoc_step(...) : A step with the name = \'" << assoc_step_name
		<< "\' being inserted into the position  = " << step_poss
		<< "." << ( type == PRE_STEP
					? (int)assoc_step_poss - num_assoc_steps(step_poss,type) - 1
					: assoc_step_poss )
		<< " has assoc_step.get() == NULL!" );
	if(running_state() == RUNNING_BEING_CONFIGURED) validate_not_curr_step(validate(step_poss));
	// Make sure an associated step with this name does not already exist.
	assoc_steps_ele_list_t &assoc_list = assoc_steps_[step_poss - 1][type];
	validate(assoc_list,assoc_step_poss,+1);
	assoc_steps_ele_list_t::iterator itr = assoc_list.begin();
	char assoc_type_name[2][10] = { "PRE_STEP" , "POST_STEP" };
	if( assoc_list.end() != ( itr = assoc_step_itr(assoc_list,assoc_step_name) ) )
		TEST_FOR_EXCEPTION(
			true, AlreadyExists
			,"Algorithm::insert_assoc_step(...) : An associated step of type = "
			<<	assoc_type_name[type]
			<< " with the name = " << assoc_step_name
			<< " already exists at step_poss = " << step_poss
			<< " and assoc_step_poss = " <<  std::distance(assoc_list.begin(),itr) + 1 );
	// insert an associated step in such a way that any container could be used.
	itr = assoc_list.begin();
	std::advance( itr, assoc_step_poss - 1 );
	assoc_list.insert( itr , assoc_steps_ele_list_ele_t(assoc_step,assoc_step_name) );
}

void Algorithm::remove_assoc_step(poss_type step_poss, EAssocStepType type, poss_type assoc_step_poss)
{
	validate_not_in_state(RUNNING);
	if(running_state() == RUNNING_BEING_CONFIGURED) validate_not_curr_step(validate(step_poss));
	validate(step_poss);
	assoc_steps_ele_list_t &assos_list = assoc_steps_[step_poss - 1][type];
	validate(assos_list,assoc_step_poss);
	assoc_steps_ele_list_t::iterator itr = assos_list.begin();
	std::advance( itr, assoc_step_poss - 1 );
	assos_list.erase( itr );
}

//  runtime configuration updating control

void Algorithm::begin_config_update()
{
	validate_in_state(RUNNING);
	saved_next_step_name_ = *next_step_name_;
	saved_curr_step_name_ = steps_[curr_step_poss_ - 1].name;
	running_state_ = RUNNING_BEING_CONFIGURED;
}

void Algorithm::end_config_update()
{
	validate_in_state(RUNNING_BEING_CONFIGURED);

	// update next_step_poss_ and next_step_name_.
	steps_t::iterator itr = step_itr(saved_next_step_name_);
	assert( itr != steps_.end() );	// the step with this name should not have been deleted or changed.
	next_step_poss_ = std::distance( steps_.begin() , itr ) + 1;
	next_step_name_ = &(*itr).name;

	// update curr_step_poss_
	itr = step_itr(saved_curr_step_name_);
	assert( itr != steps_.end() );	// the step with this name should not have been deleted or changed.
	curr_step_poss_ = std::distance( steps_.begin() , itr ) + 1;

	// inform the step objects that *this has changes.
	steps_t::iterator			s_itr = steps_.begin();
	assoc_steps_t::iterator		a_itr = assoc_steps_.begin();
	for(; s_itr != steps_.end(); ++s_itr, ++a_itr) {
		(*(*s_itr).step_ptr).inform_updated(*this);
		for(int assoc_type_i = 0; assoc_type_i < 2; ++assoc_type_i) {	// PRE_STEP, POST_STEP
			assoc_steps_ele_list_t::iterator	as_itr = (*a_itr)[assoc_type_i].begin(),
												as_itr_end = (*a_itr)[assoc_type_i].end();
			for(; as_itr != as_itr_end; ++as_itr)
				(*(*as_itr).step_ptr).inform_updated(*this);
		}
	}

	running_state_ = RUNNING;
	reconfigured_ = true;
}

// algorithmic control

void Algorithm::do_step_next(const std::string& step_name)
{
	validate_in_state(RUNNING);
	steps_t::iterator itr = step_itr_and_assert(step_name);
	next_step_poss_ = std::distance( steps_.begin() , itr ) + 1;
	next_step_name_ = &(*itr).name;
	do_step_next_called_ = true;
}

void Algorithm::do_step_next(poss_type step_poss)
{
	validate_in_state(RUNNING);
	const steps_ele_t &ele = steps_[validate(step_poss) - 1];
	next_step_poss_ = step_poss;
	next_step_name_ = &ele.name;
	do_step_next_called_ = true;
}

const std::string& Algorithm::what_is_next_step_name() const
{
	validate_in_state(RUNNING);
	return *next_step_name_;
}

Algorithm::poss_type Algorithm::what_is_next_step_poss() const
{	
	validate_in_state(RUNNING);
	return next_step_poss_;
}

bool Algorithm::do_step(const std::string& step_name)
{
	validate_in_state(RUNNING);
	return imp_do_step( std::distance( steps_.begin() , step_itr_and_assert(step_name) ) + 1 );
}

bool Algorithm::do_step(poss_type step_poss)
{
	validate_in_state(RUNNING);
	return imp_do_step(step_poss);
}

void Algorithm::terminate(bool success)
{	
	validate_in_state(RUNNING);
	terminate_status_ = success ? STATUS_TERMINATE_TRUE : STATUS_TERMINATE_FALSE;
}

// start iterations

EAlgoReturn Algorithm::do_algorithm(poss_type step_poss)
{
	using StopWatchPack::stopwatch;

	validate_in_state(NOT_RUNNING);

	track().initialize();

	try{
	
	terminate_status_ = STATUS_KEEP_RUNNING;
	running_state_ = RUNNING;

	first_k_ = state().k();
	next_step_poss_ = validate(step_poss);
	next_step_name_ = &steps_[step_poss - 1].name;
	
	// Prepair for timing algorithm
	step_times_.resize( algo_timing_ ? (num_steps()+1) * (max_iter()+1+NUM_STEP_TIME_STATS) : 0 );
	if( algo_timing_ ) {
//		step_times_[ max_iter() ] = 0.0;	// flag for statistics not calc. yet.
//		// set iteration totals to zero
//		if( step_times_[(max_iter() + 1 + 5) * num_steps()] != 0.0 )
//			std::fill_n( step_times_.begin() + (max_iter() + 1 + 5) * num_steps(), max_iter(), 0.0 );
		std::fill_n( step_times_.begin(), step_times_.size(), 0.0 );	// Try setting everything to zero?
		time_stats_computed_ = false;
	}
	stopwatch step_timer;
	stopwatch overall_timer;

	overall_timer.start();
	for(;;) {

		curr_step_poss_ = next_step_poss_;
		// Note that curr_step_poss_ may change if there is a runtime
		// change in the configuration of the steps.

		bool keep_on = true;

		// Execute the steps for this step
		
		if( algo_timing_ ) {
			step_timer.reset();
			step_timer.start();
		}

		keep_on = imp_do_step(curr_step_poss_);

		if( algo_timing_ ) {
			const double time = my_max(step_timer.stop(),-1e-50);	// negative somehow (g++ -O2 ?)
			// time for step k for the iteration
			step_times_[state().k()-first_k_+(curr_step_poss_-1)*(max_iter()+1+NUM_STEP_TIME_STATS)] = time;
			// Add to time for the full iteration
			step_times_[state().k()-first_k_+(num_steps())*(max_iter()+1+NUM_STEP_TIME_STATS)] += time;
		}

		// See if a step object called terminate(...)
		if(terminate_status_ != STATUS_KEEP_RUNNING) {
			EAlgoReturn algo_return;
			if( interrupt_status_ == STOP_END_STEP ) {
				algo_return = ( terminate_status_ == STATUS_TERMINATE_TRUE
												? INTERRUPTED_TERMINATE_TRUE
												: INTERRUPTED_TERMINATE_FALSE );
				interrupt_status_ = NOT_INTERRUPTED;
			}
			else {
				algo_return = ( terminate_status_ == STATUS_TERMINATE_TRUE
												? TERMINATE_TRUE
												: TERMINATE_FALSE );
			}
			track().output_final(*this,algo_return);
			running_state_ = NOT_RUNNING;
			return algo_return;
		}

		if(keep_on) {
			// All the step objects returned true so increment the step and loop around

			if( curr_step_poss_ == num_steps() ) {
				//
				// This is the last step in the algorithm
				//
	
				// Output this iteration
				track().output_iteration(*this);

				// Check if the maximum number of iterations has been exceeded.
				if( state().k() - first_k_ >= max_iter() ) {
					running_state_ = NOT_RUNNING;
					track().output_final(*this,MAX_ITER_EXCEEDED);
					return MAX_ITER_EXCEEDED;
				}

				// Check if the maximum runtime has been exceeded.
				if( overall_timer.read() / 60 >= max_run_time() ) {
					running_state_ = NOT_RUNNING;
					track().output_final(*this,MAX_RUN_TIME_EXCEEDED);
					return MAX_RUN_TIME_EXCEEDED;
				}

				// Set if the algorithm was interrupted
				if( interrupt_status_ == STOP_END_ITER ) {
					running_state_ = NOT_RUNNING;
					interrupt_status_ = NOT_INTERRUPTED;
					const EAlgoReturn algo_return = ( interruptTerminateReturn_
																						? INTERRUPTED_TERMINATE_TRUE
																						: INTERRUPTED_TERMINATE_FALSE );
					track().output_final(*this,algo_return);
					return algo_return;
				}

				// Transition the iteration quantities to k = k + 1
				state().next_iteration();

				// Setup to start the major loop over again
				next_step_poss_ = 1;
				next_step_name_ = &steps_[0].name;
			}
			else {
				// else just increment the step
				++next_step_poss_;
				next_step_name_ = &steps_[next_step_poss_ - 1].name;
			}
			continue;	// loop around
		}
		else {
			// some step object returned false from its do_step(..) operation so it
			// should have called do_step_next(...) to request a jump to
			// a specific operation.
			if(!do_step_next_called_)
				TEST_FOR_EXCEPTION(
					true, InvalidControlProtocal
					,"EAlgoReturn Algorithm::do_algorithm(...) :"
					" A step object returned false from its do_step(...) operation"
					" without calling do_step_next(...) to request jump to a specific"
					" step." );
			do_step_next_called_ = false;
			// just loop around and do the step that the step object requested
			// by changing next_step_poss_ by its call to do_step_next(...).
		}
	}	// end for(;;)

	}	// end try
	catch(...) {
		running_state_ = NOT_RUNNING;
		track().output_final( *this,TERMINATE_FALSE );  // This may also throw an exception?
		throw;
	}
}

// algorithm information output

void Algorithm::print_steps(std::ostream& out) const
{
	out << "\n*** Algorithm Steps ***\n\n";
	imp_print_algorithm(out,false);
	out << std::endl;
}

void Algorithm::print_algorithm(std::ostream& out) const
{
	out << "\n*** Iteration Quantities ***\n\n";
	state().dump_iter_quant(out);
	out << std::endl;
	out << "\n*** Algorithm Description ***\n\n";
	imp_print_algorithm(out,true);
	out << std::endl;
}

// Algorithm Timing.

void Algorithm::set_algo_timing( bool algo_timing ) {
	 validate_not_in_state(RUNNING);
	algo_timing_ = algo_timing;
}

bool Algorithm::algo_timing() const {
	return algo_timing_;
}

void Algorithm::print_algorithm_times( std::ostream& out ) const
{
	using std::setw;
	using std::endl;

	validate_not_in_state(RUNNING);

	if( step_times_.size() == 0 ) {
		out << "No step timing was performed\n";
		return;
	}

	const int w = 10;
	const int prec = 4;
	const int n = num_steps();					// Total steps
	const int m = state().k() - first_k_ + 1;	// Total number of iterations performed
	const int mm = max_iter()+1;				// Total number of possible iterations
	const int mmm = mm + NUM_STEP_TIME_STATS;	// total entries in a step_i row
	 
	// Print the header.
	out	<< "\n\n**************************************\n"
		<< "*** Algorithm step CPU times (sec) ***\n";

	// Print the step names.
	out	<< "\nStep names"
		<< "\n----------\n";
	{for( int i = 1; i <= n; ++i ) {
		out	<< i << ") \"" << get_step_name(i) << "\"\n";	
	}}
	out	<< n+1 << ") Iteration total\n";	
	out << endl;

	out << std::right << std::setprecision(prec);

	// Print table header
	out << setw(w) << "" << "  steps 1..." << n+1 << " ->\n\n";
	
	// print step numbers
	out	<< setw(w)	<< " iter k";
	{for( int i = 1; i <= n+1; ++i ) {
		out	<< setw(w) << i;	
	}}
	out << endl;
	out	<< setw(w)	<< "--------";
	{for( int i = 1; i <= n+1; ++i ) {
		out	<< setw(w) << "--------";	
	}}
	out << endl;
	// Print the step times.
	{for( int k = 0; k < m; ++k ) {
		out	<< setw(w)	<< first_k_ + k;
		{for( int i = 0; i < n+1; ++i ) {
			out	<< setw(w) << step_times_[k+i*mmm];	
		}}
		out << endl;
	}}

	// Compute the (1) totals for each step, the (2) average, (3) min and (4) max times
	// per iteration for each step and the (5) precentages for each step.

	compute_final_time_stats();

	// Ouput time statistics.
	
	out	<< setw(w)	<< "--------";
	{for( int i = 1; i <= n+1; ++i ) {
		out	<< setw(w) << "--------";	
	}}

	// Output the total times for each step.
	out << endl;
	out	<< setw(w)	<< "total(sec)";
	{for( int i = 0; i < n+1; ++i ) {
		const double *step_i_times = &step_times_[i*mmm];
		out	<< setw(w) << step_i_times[ mm + TIME_STAT_TOTALS_OFFSET ];	
	}}
	out << endl;

	// Output the average times per iteration
	out	<< setw(w)	<< "av(sec)/k";
	{for( int i = 0; i < n+1; ++i ) {
		const double *step_i_times = &step_times_[i*mmm];
		out	<< setw(w) << step_i_times[ mm + TIME_STAT_AV_OFFSET ];	
	}}
	out << endl;

	// Output the min times per iteration
	out	<< setw(w)	<< "min(sec)";
	{for( int i = 0; i < n+1; ++i ) {
		const double *step_i_times = &step_times_[i*mmm];
		out	<< setw(w) << step_i_times[ mm + TIME_STAT_MIN_OFFSET ];	
	}}
	out << endl;

	// Output the max times per iteration
	out	<< setw(w)	<< "max(sec)";
	{for( int i = 0; i < n+1; ++i ) {
		const double *step_i_times = &step_times_[i*mmm];
		out	<< setw(w) << step_i_times[ mm + TIME_STAT_MAX_OFFSET ];	
	}}
	out << endl;

	// Output the precentage times for each step.
	out	<< setw(w)	<< "% total";
	{for( int i = 0; i < n+1; ++i ) {
		const double *step_i_times = &step_times_[i*mmm];
		out	<< setw(w) << step_i_times[ mm + TIME_STAT_PERCENT_OFFSET ] * 100.0;	
	}}
	out << endl;


	// Print total time for entire algorithm.
	out << "------------------------------" << endl
		<< "total CPU time = " << total_time_ << " sec\n";;
}


void Algorithm::get_step_times_k( int offset, double step_times[] ) const
{
	TEST_FOR_EXCEPTION(
	  step_times_.size() == 0, std::logic_error
	  ,"Algorithm::get_step_times_k(...) : times requested, but no times calculated!"
	  );
	TEST_FOR_EXCEPTION(
	  offset > 0, std::invalid_argument
	  ,"Algorithm::get_step_times_k(...) : Can\'t get times for an iteratin that has not occured yet!."
	  );

	const int n = num_steps();					// Total steps
	const int m = state().k() - first_k_ + 1;	// Total number of iterations performed
	const int mm = max_iter()+1;				// Total number of possible iterations
	const int mmm = mm + NUM_STEP_TIME_STATS;	// total entries in a step_i row
	
	const int k = state().k() + offset;
	{for (int step = 0; step < n+1; ++step) {
		step_times[step] = step_times_[step*mmm + k];		
	}}

}

void Algorithm::get_final_step_stats( size_t step, double* total, double* average, double* min, double* max, double* percent) const
{
	// Compute the (1) totals for each step, the (2) average, (3) min and (4) max times
	// per iteration for each step and the (5) precentages for each step.
	compute_final_time_stats();

	const int n = num_steps();					// Total steps
	const int m = state().k() - first_k_ + 1;	// Total number of iterations performed
	const int mm = max_iter()+1;				// Total number of possible iterations
	const int mmm = mm + NUM_STEP_TIME_STATS;	// total entries in a step_i row

	double* step_i_times = &const_cast<step_times_t&>(step_times_)[step*mmm];
	if (total) {
		*total   = step_i_times[mm + TIME_STAT_TOTALS_OFFSET];
	}
	if (average) {
		*average = step_i_times[mm + TIME_STAT_AV_OFFSET];
	}
	if (min) {
		*min     = step_i_times[mm + TIME_STAT_MIN_OFFSET];
	}
	if (max) {
		*max     = step_i_times[mm + TIME_STAT_MAX_OFFSET];
	}
	if (percent) {
		*percent = step_i_times[mm + TIME_STAT_PERCENT_OFFSET];
	}
}


void Algorithm::compute_final_time_stats() const
{
	if (!time_stats_computed_) {
		time_stats_computed_ = true;
		
		const int n = num_steps();					// Total steps
		const int m = state().k() - first_k_ + 1;	// Total number of iterations performed
		const int mm = max_iter()+1;				// Total number of possible iterations
		const int mmm = mm + NUM_STEP_TIME_STATS;	// total entries in a step_i row
		
		// compute totals for each step (1...n) and the full iteration (n+1)
		double &_total_time = const_cast<double&>(total_time_);
		_total_time = 0.0;
		
		{for( int i = 0; i < n+1; ++i ) {
			double *step_i_times = &const_cast<step_times_t&>(step_times_)[i*mmm];
			// compute total step times (and total algorithm time)
			const double
				step_time = std::accumulate( step_i_times, step_i_times + m, (double)0.0 );
			if(i < n)
				_total_time += step_time;
			step_i_times[ mm + TIME_STAT_TOTALS_OFFSET ] = step_time;
			// compute average per step.
			step_i_times[ mm + TIME_STAT_AV_OFFSET ] = step_time / m;
			// compute min per step
			step_i_times[ mm + TIME_STAT_MIN_OFFSET ]= *std::min_element( step_i_times, step_i_times + m );
			// compute max per step
			step_i_times[ mm + TIME_STAT_MAX_OFFSET ]= *std::max_element( step_i_times, step_i_times + m );
		}}
		
		{for( int i = 0; i < n+1; ++i ) {
			double *step_i_times = &const_cast<step_times_t&>(step_times_)[i*mmm];
			// compute fractions for each step.
			step_i_times[ mm + TIME_STAT_PERCENT_OFFSET ]
				= step_i_times[ mm + TIME_STAT_TOTALS_OFFSET ] / total_time_;
		}}
	}
}

// private

void Algorithm::validate_in_state(ERunningState running_state) const {
	const char running_state_name[3][25] = { "NOT_RUNNING" , "RUNNING", "RUNNING_BEING_CONFIGURED" };
	if(running_state_ != running_state)
		TEST_FOR_EXCEPTION(
			true, InvalidRunningState
			,"Algorithm::validate_in_state(...) : The condition running_state() == "
			<< running_state_name[running_state_] << " has been violated with "
			<< " running_state() = " << running_state_name[running_state] );
}

void Algorithm::validate_not_in_state(ERunningState running_state) const {
	const char running_state_name[3][25] = { "NOT_RUNNING" , "RUNNING", "RUNNING_BEING_CONFIGURED" };
	if(running_state_ == running_state)
		TEST_FOR_EXCEPTION(
			true, InvalidRunningState
			,"Algorithm::validate_not_in_state(...) : The condition running_state() != "
			<< running_state_name[running_state_] << " has been violated" );
}

void Algorithm::validate_not_curr_step(poss_type step_poss) const {
	if(step_poss == curr_step_poss_)
		TEST_FOR_EXCEPTION(
			true, InvalidConfigChange
			,"Algorithm::validate_not_curr_step(step_poss="<<step_poss<<") : "
			"Error, You can not modify the step being currently executed" );
}

void Algorithm::validate_not_next_step(const std::string& step_name) const {
	if( step_name == saved_next_step_name_ )
		TEST_FOR_EXCEPTION(
			true, InvalidConfigChange,
			"Algorithm::validate_not_next_step(step_name): "
			"Error, You can not modify name or remove the step given by "
			"step_name = what_is_next_name() = " << step_name );
}

Algorithm::steps_t::iterator Algorithm::step_itr_and_assert(const std::string& step_name)
{
	steps_t::iterator itr = step_itr(step_name);
	if(itr == steps_.end())
		TEST_FOR_EXCEPTION(
			true, DoesNotExist
			,"Algorithm::step_itr(...) : A step with the name "
			<< step_name << " does not exist." );
	return itr;	
}

Algorithm::steps_t::const_iterator Algorithm::step_itr_and_assert(const std::string& step_name) const
{
	steps_t::const_iterator itr = step_itr(step_name);
	if(itr == steps_.end())
		TEST_FOR_EXCEPTION(
			true, DoesNotExist
			,"Algorithm::step_itr(...) : A step with the name "
			<< step_name << " does not exist." );
	return itr;	
}

bool Algorithm::imp_do_step(poss_type step_poss) {
	curr_step_poss_ = step_poss;
	// do the pre steps in order
	if( !imp_do_assoc_steps(PRE_STEP) ) return false;
	// do the main step
	if( !steps_[curr_step_poss_-1].step_ptr->do_step(*this, curr_step_poss_, DO_MAIN_STEP, 0) ) return false;
	// do the post steps in order
	if( !imp_do_assoc_steps(POST_STEP) ) return false;
	// if you get here all the pre steps, step, and post steps returned true.
	look_for_interrupt();
	if( interrupt_status_ == STOP_END_STEP ) {
		terminate( interruptTerminateReturn_ );
		return false;
	}
	return true;
}

bool Algorithm::imp_do_assoc_steps(EAssocStepType type) {
	assoc_steps_ele_list_t				*assoc_list	= &assoc_steps_[curr_step_poss_ - 1][type];
	assoc_steps_ele_list_t::iterator	itr			= assoc_list->begin();
	int									n			= assoc_list->size();
	for(int i = 1; i <= n; ++itr, ++i) {
		if(reconfigured_) {
			// The associated step just has reconfigured *this
			// so we must update our pointers and iterators.
			// Since it is not allowed for this step or its associated steps
			// to have been changed, the next associated step to
			// execute will not change.
			assoc_list	= &assoc_steps_[curr_step_poss_ - 1][type];
			itr			= assoc_list->begin();
			std::advance( itr, i - 1 );
			reconfigured_ = false;	// This works as long as no one else needs to know
									// if *this has been reconfigured.
		}
		if( !(*(*itr).step_ptr).do_step(*this, curr_step_poss_, do_step_type(type), i) ) return false;
	}
	return true;	// All the associated steps returned true.
}

void Algorithm::imp_print_algorithm(std::ostream& out, bool print_steps) const
{
	const std::string leading_str = "    ";
	
	steps_t::const_iterator				s_itr = steps_.begin();
	assoc_steps_t::const_iterator		a_itr = assoc_steps_.begin();
	poss_type step_i = 1;
	for(; step_i <= num_steps(); ++step_i, ++s_itr, ++a_itr) {
		// list pre_steps (e.q. 2.-3, 2.-2, 2.-1)
		const assoc_steps_ele_list_t &pre_steps = (*a_itr)[PRE_STEP];
		assoc_steps_ele_list_t::const_iterator pre_step_itr = pre_steps.begin();
		for(int pre_step_i = - pre_steps.size(); pre_step_i < 0; ++pre_step_i, ++pre_step_itr) {
			out		<< step_i << "." << pre_step_i << ". \""
					<< (*pre_step_itr).name << "\"\n"
					<< leading_str << "(" << typeid(*(*pre_step_itr).step_ptr).name() << ")\n";
			if(print_steps) {
				(*(*pre_step_itr).step_ptr).print_step( *this, step_i, DO_PRE_STEP
					, pre_steps.size()+pre_step_i+1, out, leading_str );
				out << std::endl;
			}
		}
		// The main step.
		out		<< step_i << ". \"" << (*s_itr).name
				<< "\"\n"
				<< leading_str << "(" << typeid(*(*s_itr).step_ptr).name() << ")\n";
		if(print_steps) {
			(*(*s_itr).step_ptr).print_step( *this, step_i, DO_MAIN_STEP, 0, out, leading_str );
			out << std::endl;
		}
		// list post_steps (e.q. 2.1, 2.2, 2.3)
		const assoc_steps_ele_list_t &post_steps = (*a_itr)[POST_STEP];
		assoc_steps_ele_list_t::const_iterator post_step_itr = post_steps.begin();
		for(int post_step_i = 1; post_step_i <= post_steps.size(); ++post_step_i, ++post_step_itr) {
			out		<< step_i << "." << post_step_i << ". \""
					<< (*post_step_itr).name << "\"\n"
					<< leading_str << "(" << typeid(*(*post_step_itr).step_ptr).name() << ")\n";
			if(print_steps) {
				(*(*post_step_itr).step_ptr).print_step( *this, step_i, DO_POST_STEP, post_step_i
					, out, leading_str );
				out << std::endl;
			}
		}
	}
	if(print_steps) {
		out
			<< step_i << ". \"Major Loop\" :\n"
			<< "    if k >= max_iter then\n"
			<< "        terminate the algorithm\n"
			<< "    elseif run_time() >= max_run_time then\n"
			<< "        terminate the algorithm\n"
			<< "    else\n"
			<< "        k = k + 1\n"
			<< "        goto 1\n"
			<< "    end\n";
	}
}

// validate poss

Algorithm::poss_type Algorithm::validate(poss_type step_poss, int past_end) const
{
    
	TEST_FOR_EXCEPTION(
		step_poss < 1 || steps_.size() + past_end < step_poss, DoesNotExist
		,"Algorithm::validate(step_poss) : The step_poss = " << step_poss
		<< " is not in range of 1 to " << steps_.size() + past_end );
	return step_poss;
}	

Algorithm::poss_type Algorithm::validate(const assoc_steps_ele_list_t& assoc_list
	, poss_type assoc_step_poss, int past_end) const
{
	TEST_FOR_EXCEPTION(
		assoc_step_poss < 1 || assoc_list.size() + past_end < assoc_step_poss, DoesNotExist
		,"Algorithm::validate(assoc_list,assoc_step_poss) : The assoc_step_poss = "
		<< assoc_step_poss << " is not in range of 1 to " << assoc_list.size() + past_end );
	return assoc_step_poss;
}

void Algorithm::look_for_interrupt()
{
	//
	// Get the mode of aborting from the user!
	//
	if( interrupt_called_ && !processed_user_interrupt_ && proc_rank_ == 0 ) {
		//
		// Get the response from the user
		//
		enum EResponse { R_STOP_END_STEP, R_STOP_END_ITER, R_ABORT_NOW };
		EResponse response = R_ABORT_NOW;
		const int max_tries = 3;
		bool valid_response = false;
		for( int tries = 0; !valid_response && tries < max_tries; ++tries ) {
			std::cerr
				<< "\np="<<proc_rank_<<": IterationPack::Algorithm::interrupt(): Received signal SIGINT: Do you want to:\n"
				<< "  (a) Gracefully terminate the algorithm at the end of this step?\n"
				<< "  (b) Gracefully terminate the algorithm at the end of this iteration?\n"
				<< "  (c) Abort the program immediately?\n"
				<< "Answer a, b or c ? ";
			char c_response;
			std::cin >> c_response;
			if( c_response == 'a' || c_response == 'b' ) {
				if( c_response == 'a')
					response = R_STOP_END_STEP;
				else
					response = R_STOP_END_ITER;
				std::cerr
					<< "\np="<<proc_rank_<<": Terminate the algorithm with true (t) or false (f) ? ";
				std::cin >> c_response;
				if( c_response == 't' ) {
					interruptTerminateReturn_ = true;
					valid_response = true;
				}
				else if( c_response == 'f' ) {
					interruptTerminateReturn_ = false;
					valid_response = true;
				}
				else {
					std::cerr	<< "p="<<proc_rank_<<": Invalid response! Expecting \'t\' or \'f\'\n";
				}
				if(valid_response) {
					// ToDo: Write a file if needed (must syn up processors here)!
				}
			}
			else if( c_response == 'c' ) {
					response = R_ABORT_NOW;
					valid_response = true;
			}
			else {
				std::cerr	<< "\np="<<proc_rank_<<": Invalid response! Expecting \'a\', \'b\' or \'c\'\n";
			}
			std::cerr << std::endl;
		}
		if(!valid_response) {
			std::cerr << "p="<<proc_rank_<<": Three strikes, you are out!\n";
		}
		//
		// Interpret the response
		//
		switch(response) {
			case R_STOP_END_STEP: {
				interrupt_status_ = STOP_END_STEP;
				break;
			}
			case R_STOP_END_ITER: {
				interrupt_status_ = STOP_END_ITER;
				break;
			}
			case R_ABORT_NOW: {
				interrupt_status_ = ABORT_PROGRAM;
				break;
			}
			default: {
				assert(0);
			}
		}
		processed_user_interrupt_ = true;
	}
	//
	// Make sure that all of the processes get the same
	// response
	//
#ifdef HAVE_MPI
	const bool query_for_interrupt = true; // ToDo: Make this an external option!
	if( num_proc_ > 1 && query_for_interrupt ) {
		//
		// Here we will do a global reduction to see of a processor has
		// recieved an interrupt.  Here we will do a sum since only the
		// root processes values should be nonzero.
		//
		int sendbuf[2] = { 0, 0 };
		int recvbuf[2] = { 0, 0 };
		if( proc_rank_ == 0 ) {
			sendbuf[0] = (int)interrupt_status_;
			sendbuf[1] = interruptTerminateReturn_ ? 1 : 0;
		}
		// Note: this global reduction will synchronize all of the processors!
#ifdef ITERATION_PACK_ALGORITHM_SHOW_MPI_DEBUG_INFO
		std::cerr	<< "\np="<<proc_rank_<<": IterationPack::Algorithm::interrupt(): Calling MPI_Allreduce(...) ...\n";
#endif
		MPI_Allreduce(
			sendbuf                  // sendbuf
			,recvbuf                 // recvbuf
			,2                       // count
			,MPI_INT                 // datatype
			,MPI_SUM                 // op
			,MPI_COMM_WORLD          // comm (ToDo: Make more general?)
			);
#ifdef ITERATION_PACK_ALGORITHM_SHOW_MPI_DEBUG_INFO
		std::cerr
			<< "\np="<<proc_rank_<<": IterationPack::Algorithm::interrupt(): After MPI_Allreduce(...)"
			<< "\np="<<proc_rank_<<": recvbuf[0] = " << recvbuf[0] << ", recvbuf[1] = " << recvbuf[1] << std::endl;
#endif
		// Set interrupt_status_
		switch( (EInterruptStatus)recvbuf[0] ) {
			case NOT_INTERRUPTED:
				interrupt_status_ = NOT_INTERRUPTED;
				break;
			case STOP_END_STEP:
				interrupt_status_ = STOP_END_STEP;
				break;
			case STOP_END_ITER:
				interrupt_status_ = STOP_END_ITER;
				break;
			case ABORT_PROGRAM:
				interrupt_status_ = ABORT_PROGRAM;
				break;
			default:
				std::cerr
					<< "p=" << proc_rank_ << ": Algorithm::look_for_interrupt(): Error, the globally reduced value of "
					"recvbuf[0] = " << recvbuf[0] << " is not valid!";
				abort();
		}
		// Set interruptTerminateReturn_
		interruptTerminateReturn_ = ( recvbuf[1] == 0 ? false : true );
	}
	//
	// Abort the program now if the user did not already press Ctrl-C again!
	//
	if( interrupt_status_ == ABORT_PROGRAM ) {
		if( proc_rank_ == 0 ) {
			std::cerr << "\np=" << proc_rank_ << ": Aborting the program now!\n";
		}
		abort();
	}
#endif
}

// static

int Algorithm::num_proc_ = 0; // Flag that an Algorithm object has not be allocated yet!

int Algorithm::proc_rank_ = 0;

bool Algorithm::interrupt_called_ = false;

bool Algorithm::processed_user_interrupt_ = false;

Algorithm::EInterruptStatus Algorithm::interrupt_status_ = Algorithm::NOT_INTERRUPTED;

bool Algorithm::interruptTerminateReturn_ = true;

void Algorithm::interrupt()
{
	processed_user_interrupt_ = false;
#ifdef ITERATION_PACK_ALGORITHM_SHOW_MPI_DEBUG_INFO
	std::cerr	<< "\np="<<proc_rank_<<": IterationPack::Algorithm::interrupt() called!\n";
#endif
	//
	// This function assumes that every process will recieve the same
	// signal which I found to be the case with MPICH.  I am not clear
	// what the MPI standard says about interrupts so I can not
	// guarantee that this is 100% portable.  If other behavior is
	// needed, this will have to be compiled in differently.
	//
	// Note: I have found that on MPICH that you can not guarantee that
	// only a single signal will be sent to a slave process so this
	// function will ignore interupts for slave processes.
	//
	// Note that you have to be very careful what you do inside of a
	// signal handler and in general you should only be setting flags or
	// aborting.
	//
	// See if an algorithm is possibly even running yet!
	//
	if( num_proc_ == 0 ) {
		if( proc_rank_ == 0 )
			std::cerr
				<< "\nIterationPack::Algorithm::interrupt(): Received signal SIGINT but an Algorithm "
				<< "object has not been allocated yet and no algorithm is running.\n"
				<< "\nAborting the program now!\n";
		abort();
		return;  // Should not be called!
	}
	//
	// See if we are going to query for an interrupt when running in MPI mode
	//
	const bool query_for_interrupt = true; // ToDo: Make this an external option!
	if( !query_for_interrupt && num_proc_ > 1 ) {
		if( proc_rank_ == 0 )
			std::cerr
				<< "\nIterationPack::Algorithm::interrupt(): Received signal SIGINT but num_proc = "
				<< num_proc_ << " > 1 and query_for_interrupt = false so:\n"
				<< "\nAborting the program now!\n";
		abort();
		return;  // Should not be called!
	}
	//
	// Remember that this interrupt has been called!
	//
	if( proc_rank_ == 0 )
		std::cerr
			<< "\nIterationPack::Algorithm::interrupt(): Received signal SIGINT.  "
			<< "Wait for the end of the current step and respond to an interactive query,  "
			<< "kill the process by sending another signal (i.e. SIGKILL).\n";
	interrupt_called_ = true;
}

} // end namespace IterationPack
