// ////////////////////////////////////////////////////////////////////////////////////
// Algorithm.h

#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <assert.h>

#include <string>
#include <deque>
#include <list>
#include <vector>
#include <sstream>
#include <algorithm>

#include "AlgorithmState.h"
#include "AlgorithmTrack.h"
#include "AlgorithmStep.h"

#include "Misc/include/ref_count_ptr.h"

namespace GeneralIterationPack {

// ToDo: 7/31/98: Finish documentation.

///
/** Acks as the central hub for an iterative algorithm.
  *
  * This class is the center for a framework for iterative algorithms.
  * These iterative algorithms are of the form:
  *
	\begin{verbatim}

	Step1------>Step2------>Step3------>Step4------>Step5
	 /|\         /|\         /|\          |           |
	  |           |           |_Minor L 1_|           |
	  |           |                       |           |
	  |           |_____Minor Loop 2______|           |
	  |                                               |
	  |_______________Major Loop (k = k+1)____________|

	\end{verbatim}
  *
  * For the typical iteration the steps are executed sequantially from Step1 to Step2
  * and then control loops around the Major Loop to Step1 again.
  * Durring some iterations however
  * Minor Loop 1 may be executed several times before control is continued alone the
  * major loop.  The same may also apply to Minor Loop 2.
  *
  * To allow for greater algorithmic control any step object can take over the role of
  * #Algorithm# and take over complete control the algorithm.  For examle, Step 4 may
  * need to Execute Step1->Step3->Step2->Step5 before returning algorithmic control to
  * #Algorithm#.
  *
  * #Algorithm# executes the steps of the algorithm through step objects of the
  * base type #AlgorithmStep#.  In addition to major step objects as shown above
  * there are also PreStep and PostStep objects.
  * These are steps that are intimatly associated with a major step object and
  * will always (well almost always) be exectuted alone with a major step.
  *
  * ToDo: Finish documentation.
  */
class Algorithm {
public:

	/** @name Public types */
	//@{

	///
	typedef ReferenceCountingPack::ref_count_ptr<AlgorithmState>		state_ptr_t;
	///
	typedef ReferenceCountingPack::ref_count_ptr<AlgorithmTrack>		track_ptr_t;
	///
	typedef ReferenceCountingPack::ref_count_ptr<AlgorithmStep>			step_ptr_t;
	///
	typedef size_t														poss_type;
	///
	enum { DOES_NOT_EXIST = 1000 };	// never be that many steps
	///
	enum ERunningState { NOT_RUNNING = 0, RUNNING = 1, RUNNING_BEING_CONFIGURED = 2 };

	/// Thrown if name or id does not exist
	class DoesNotExist : public std::logic_error
	{public: DoesNotExist(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown if name already exists
	class AlreadyExists : public std::logic_error
	{public: AlreadyExists(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown if an invalid control protocal is used.
	class InvalidControlProtocal : public std::logic_error
	{public: InvalidControlProtocal(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown if a member function is called while #this# is in an invalid running state..
	class InvalidRunningState : public std::logic_error
	{public: InvalidRunningState(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown if a member function is called while #this# is in an invalid running state..
	class InvalidConfigChange : public std::logic_error
	{public: InvalidConfigChange(const std::string& what_arg) : std::logic_error(what_arg) {}};

	//@}

	///
	/** Constructs an algorithm with no steps and a default of max_iter() == 100.
	  *
	  */
	Algorithm();

	///
	virtual ~Algorithm();

	/** @name «std comp» members for state */
	//@{

	///
	void set_state(const state_ptr_t& state);
	///
	state_ptr_t& get_state();
	///
	const state_ptr_t& get_state() const;
	///
	AlgorithmState& state();
	///
	const AlgorithmState& state() const;

	//@}

	/** @name «std comp» members for track */
	//@{

	///
	void set_track(const track_ptr_t& track);
	///
	track_ptr_t& get_track();
	///
	const track_ptr_t& get_track() const;
	///
	AlgorithmTrack& track();
	///
	const AlgorithmTrack& track() const;

	//@}

	///
	ERunningState running_state() const;

	/** @name maximum iterations */
	//@{
	
	///
	virtual void max_iter(size_t max_iter);
	///
	virtual size_t max_iter() const;

	//@}

	/** @name maximum runtime (in min).
	  *
	  * The runtime is checked at the end of each iteration and if it exceeds
	  * this value then the algorithm is terminated.
	  */
	//@{
	
	///
	virtual void max_run_time(double max_iter);
	///
	virtual double max_run_time() const;

	//@}

	/** @name step information / access
	  *
	  * These functions provide information as to the number of major steps
	  * , their possitions given their names and their names given their
	  * possitions.
	  *
	  * In addition, access is given to the step objects themselves through
	  * the ref_count_ptr<...> objects that are used to manage thier memory.
	  * Using this type of direct access allows clients to take over memory
	  * management if needed and to call the step objects in any order
	  * and thereby taking over control of the algorithm.
	  *
	  * These functions can be invoked in any state of the algorithm.
	  */
	//@{

	/// Return the number of main steps
	virtual int num_steps() const;

	///
	/** Return the possition in the major loop of a named step.
	  *
	  * If a step with this name does not exist then the value
	  * DOES_NOT_EXIST will be returned.
	  */
	virtual poss_type get_step_poss(const std::string& step_name) const;

	///
	/** Return the name of a step given its possition.
	  *
	  * Preconditions:\begin{itemize}
	  * \item #1 <= step_poss && step_poss <= num_steps()# (throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual const std::string& get_step_name(poss_type step_poss) const;

	///
	/** Return the ref_count_ptr<...> object for the step object at step_poss.
	  *
	  * Preconditions:\begin{itemize}
	  * \item #1 <= step_poss && step_poss <= num_steps()# (throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual step_ptr_t& get_step(poss_type step_poss);

	///
	virtual const step_ptr_t& get_step(poss_type step_poss) const;

	//@}

	/** @name pre/post step information / access */
	//@{

	///
	/** Return the number of pre or post steps for the main step step_poss.
	  *
	  * Preconditions:\begin{itemize}
	  * \item #1 <= step_poss && step_poss <= num_steps()# (throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual int num_assoc_steps(poss_type step_poss, EAssocStepType type) const;

	///
	/** Return the possition of the pre or post step for the main step_poss.
	  *
	  * If a pre or post step does not exist with the name #assoc_step_name#
	  * then a value of DOES_NOT_EXIST will be retruned.
	  *
	  * Preconditions:\begin{itemize}
	  * \item #1 <= step_poss && step_poss <= num_steps()# (throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual poss_type get_assoc_step_poss(poss_type step_poss, EAssocStepType type
		,const std::string& assoc_step_name) const;

	///
	/** Return the name of the pre or post step at step_poss and at assoc_step_poss.
	  *
	  * Preconditions:\begin{itemize}
	  * \item #1 <= step_poss && step_poss <= num_steps()# (throw #DoesNotExist#)
	  * \item #1 <= assoc_step_poss && assoc_step_poss <= num_assoc_steps(step_poss,type)#
	  *		(throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual const std::string& get_assoc_step_name(poss_type step_poss, EAssocStepType type
		, poss_type assoc_step_poss) const;

	///
	/** Return the ref_count_ptr<...> object for the associated step object at step_poss
	  * and assoc_step_poss.
	  *
	  * Preconditions:\begin{itemize}
	  * \item #1 <= step_poss && step_poss <= num_steps()# (throw #DoesNotExist#)
	  * \item #1 <= assoc_step_poss && assoc_step_poss <= num_assoc_steps(step_poss,type)#
	  *		(throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual step_ptr_t& get_assoc_step(poss_type step_poss, EAssocStepType type
		, poss_type assoc_step_poss);

	///
	virtual const step_ptr_t& get_assoc_step(poss_type step_poss, EAssocStepType type
		, poss_type assoc_step_poss) const;

	//@}

	/** @name step manipulation */
	//@{

	///
	/** Insert a step object with the name #step_name# into the possition #step_poss#.
	  *
	  * All the steps at and after #step_poss# are pushed back one possition unless
	  * #step_poss == num_steps() + 1# in which case the new step is appended to the end.
	  * Initiaily this step will have no pre or post steps associated with it.
	  * 
	  * Preconditions:\begin{itemize}
	  * \item #running_state() != RUNNING]# (throw #InvalidRunningState#)
	  * \item #1 <= ste_poss && step_poss <= num_steps() + 1# (throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual void insert_step(poss_type step_poss, const std::string& step_name, const step_ptr_t& step);

	///
	/** Change the name of an existing step.
	  *
	  * None of the pre or post steps for the existing step are changes.
	  *
	  * Preconditions:\begin{itemize}
	  * \item #running_state() != RUNNING]# (throw #InvalidRunningState#)
	  * \item #1 <= poss && poss <= num_steps()# (throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual void change_step_name(poss_type step_poss, const std::string& new_name);

	///
	/** Replace the step object of an existing step.
	  *
	  * None of the pre or post steps for the existing step are changes.
	  *
	  * Preconditions:\begin{itemize}
	  * \item #running_state() != RUNNING]# (throw #InvalidRunningState#)
	  * \item #1 <= poss && poss <= num_steps()# (throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual void replace_step(poss_type step_poss, const step_ptr_t& step);

	///
	/** Remove an existing step object and all of its pre and post steps.
	  *
	  * All of the steps after #step_poss# will have thier possitions
	  * decreased by one.
	  *
	  * Preconditions:\begin{itemize}
	  * \item #running_state() != RUNNING]# (throw #InvalidRunningState#)
	  * \item #1 <= poss && poss <= num_steps()# (throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual void remove_step(poss_type step_poss);

	//@}

	/** @name pre/post step manipulation */
	//@{

	///
	/** Insert an pre or post step into for the main step step_poss into the possition
	  * assoc_step_poss.
	  *
	  * All of the pre or post steps at and after #assoc_step_poss# will be pushed back
	  * one possition.
	  *
	  * Preconditions:\begin{itemize}
	  * \item #running_state() != RUNNING]# (throw #InvalidRunningState#)
	  * \item #1 <= step_poss && step_poss <= num_steps() + 1# (throw #DoesNotExist#)
	  * \item #1 <= assoc_step_poss && assoc_step_poss <= num_assoc_steps(step_poss,type) + 1#
	  *		(throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual void insert_assoc_step(poss_type step_poss, EAssocStepType type, poss_type assoc_step_poss
		, const std::string& assoc_step_name, const step_ptr_t& assoc_step);

	///
	/** Remove an pre or post step for the main step step_poss in the possition
	  * assoc_step_poss.
	  *
	  * All of the pre or post steps after #assoc_step_poss# will be pushed forward
	  * one possition to fill in the hole.
	  *
	  * Preconditions:\begin{itemize}
	  * \item #running_state() != RUNNING]# (throw #InvalidRunningState#)
	  * \item #1 <= step_poss && step_poss <= num_steps()# (throw #DoesNotExist#)
	  * \item #1 <= assoc_step_poss && assoc_step_poss <= num_assoc_steps(step_poss,type)#
	  *		(throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual void remove_assoc_step(poss_type step_poss, EAssocStepType type, poss_type assoc_step_poss);

	//@}

	/** @name runtime configuration updating control */
	//@{

	///
	/** Changes from running_state() == RUNNING to running_state() == RUNNING_BEING_CONFIGURED.
	  *
	  * Must be called before the algorithm's configuration can be changed while it is running.
	  *
	  * Preconditions:\begin{itemize}
	  * \item #running_state() == RUNNING]# (throw #InvalidRunningState#)
	  * \end{itemize}
	  *
	  * Postconditions:\begin{itemize}
	  * \item #running_state() == RUNNING_BEING_CONFIGURED]#
	  * \end{itemize}
	  */
	virtual void begin_config_update();

	///
	/** Changes from running_state() == RUNNING_BEING_CONFIGURED to running_state() == RUNNING.
	  *
	  * Must be called after the algorithm's configuration can be changed while it is running.
	  *
	  * Preconditions:\begin{itemize}
	  * \item #running_state() == RUNNING_BEING_CONFIGURED]# (throw #InvalidRunningState#)
	  * \end{itemize}
	  *
	  * Postconditions:\begin{itemize}
	  * \item #running_state() == RUNNING]#
	  * \end{itemize}
	  */
	virtual void end_config_update();

	//@}

	/** @name algorithmic control */
	//@{

	///
	/** Called by step objects to set the step (given its name) that #this# will envoke the next time
	  * #this# calls a step.
	  *
	  * Preconditions:\begin{itemize}
	  * \item running_state() == RUNNING# (throw #InvalidRunningState#)
	  * \item #get_step_poss(step_name) != DOES_NOT_EXIST# (throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual void do_step_next(const std::string& step_name);

	///
	/** Called by step objects to set the step (given its possition) that #this# will envoke the next time
	  * #this# calls a step.
	  *
	  * Preconditions:\begin{itemize}
	  * \item running_state() == RUNNING# (throw #InvalidRunningState#)
	  * \item #1 <= step_poss && step_poss <= num_steps()# (throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual void do_step_next(poss_type step_poss);

	///
	/** Returns the name of the next step #this# will call the next time it calls a step.
	  *
	  * Preconditions:\begin{itemize}
	  * \item running_state() != NOT_RUNNING# (throw #InvalidRunningState#)
	  * \end{itemize}
	  */
	virtual const std::string& what_is_next_step_name() const;

	///
	/** Returns the possition of the next step #this# will call the next time it calls a step.
	  *
	  * Preconditions:\begin{itemize}
	  * \item running_state() != NOT_RUNNING# (throw #InvalidRunningState#)
	  * \end{itemize}
	  */
	virtual poss_type what_is_next_step_poss() const;

	///
	/** Calls #do_step(...)# on all of the pre step objects the step object and the post step objects
	  * in order for the step named #step_name#.
	  *
	  * This operation is called by step objects that need to take over control of the algorithm
	  * at some point.
	  *
	  * If any of the of the pre or post objects or the step object returns false, then this
	  * operation immediatly returns false.  It is assumed that if any step object returns
	  * false from do step that it has either also called #terminate(...)# or #do_step_next(...)#.
	  *
	  * Preconditions:\begin{itemize}
	  * \item running_state() == RUNNING# (throw #InvalidRunningState#)
	  * \item #get_step_poss(step_name) != DOES_NOT_EXIST# (throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual bool do_step(const std::string& step_name);

	///
	/** Call #do_step(...)# on all of the pre step objects the step object and the post step objects
	  * in order for the step in the possition #step_poss#.
	  *
	  * This operation is called by step objects that need to take over control of the algorithm
	  * at some point.
	  *
	  * If any of the of the pre or post objects or the step object returns false, then this
	  * operation immediatly returns false.  It is assumed that if any step object returns
	  * false from do step that it has either also called #terminate(...)# or #do_step_next(...)#. 
	  *
	  * Preconditions:\begin{itemize}
	  * \item running_state() == RUNNING# (throw #InvalidRunningState#)
	  * \item #1 <= step_poss && step_poss <= num_steps()# (throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual bool do_step(poss_type step_poss);

	///
	/** Called by step objects to terminate the algorithm.
	  *
	  * Calling with #success == true# cause #do_algorithm(...)# to completely return #TERMINATE_TRUE#
	  * and with #success == false# return  #TERMINATE_FALSE#.
	  *
	  * Preconditions:\begin{itemize}
	  * \item running_state() == RUNNING# (throw #InvalidRunningState#)
	  * \end{itemize}
	  */
	virtual void terminate(bool success);

	//@}

	/** @name start iterations */
	//@{

	///
	/** Called by clients to begin an algorithm.
	  *
	  * This operation acts as the central hub for the algorithm.  It calls the #do_step(i)#
	  * each #i# = 1,...,#num_steps()# and then loops around again for the major loop.  If
	  * #do_step(i)# returns false then it goes executes the step specified by the 
	  * #do_step_next(...)# operation which the step object supposivly called.  If a step
	  * object returns false but does not call #do_step_next(...)# to specify a step to
	  * jump to, then #this# will throw an #InvalidControlProtocal# exception.
	  *
	  * At the end of each iteration #this# calls #track().output_iteration(*this)# and
	  * #state().next_iteration()#.  It then checks if #state.k() - k_start# >= #max_iter()#.
	  * If it is then the #do_algorithm(...)# immediatly terminates with a value of
	  * #MAX_ITER_EXCEEDED#.
	  *
	  * Any step object can cause the algorithm to terminate by calling #terminate(success)#.
	  * This operation will then immediatly return #TERMINATE_TRUE# if #success == true#
	  * and #TERMINATE_FALSE#  if #success == false#.
	  *
	  * The algorithm starts on the step specified with #step_poss#.
	  *
	  * Preconditions:\begin{itemize}
	  * \item running_state() == NOT_RUNNING# (throw #InvalidRunningState#)
	  * \item #1 <= step_poss && step_poss <= num_steps()# (throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual EAlgoReturn do_algorithm(poss_type step_poss = 1);

	//@}

	/** @name algorithm information output */
	//@{

	///
	/**
	  *
	  */
	virtual void print_steps(std::ostream& out) const;

	///
	/**
	  *
	  */
	virtual void print_algorithm(std::ostream& out) const;

	//@}

	/** @name Algorithm Timing.
	  */
	//@{

	///
	/** Causes algorithm to be timed.
	  *
	  * Call with #algo_timing == true# before #do_algorithm(...)# to have the algorithm timed.
	  *
	  * Do not call when algorithm is running.
	  */
	virtual void set_algo_timing( bool algo_timing );

	///
	virtual bool algo_timing() const;

	///
	/** Outputs table of times for each step, cummulative times and other
	  * statistics.
	  *
	  * Call after #do_algorithm(...)# has executed to get a table
	  * of times.
	  *
	  * Do not call when algorithm is running.
	  */
	virtual void print_algorithm_times( std::ostream& out ) const;

	//@}

private:
	// /////////////////////////////////////////////////////
	// Private types

	///
	template<class T_Step_ptr>
	struct step_ptr_and_name {
		///
		step_ptr_and_name(const T_Step_ptr& _step_ptr
				, const std::string& _name )
			: step_ptr(_step_ptr), name(_name)
		{}
		///
		T_Step_ptr step_ptr;
		//
		std::string name;
	};	// end struct step_ptr_and_name

	///
	typedef step_ptr_and_name<step_ptr_t>			steps_ele_t;
	///
	typedef std::deque<steps_ele_t>					steps_t;

	///
	typedef step_ptr_and_name<step_ptr_t>			assoc_steps_ele_list_ele_t;
	///
	typedef std::list<assoc_steps_ele_list_ele_t>	assoc_steps_ele_list_t;
	///
	struct assoc_steps_ele_t {
		///
		assoc_steps_ele_list_t& operator[](int i)
		{ return assoc_steps_lists_[i]; }
		///
		const assoc_steps_ele_list_t& operator[](int i) const
		{ return assoc_steps_lists_[i]; }
	private:
		assoc_steps_ele_list_t assoc_steps_lists_[2];
	};

//	typedef assoc_steps_ele_list_t[2]				assoc_steps_ele_t; // PRE_STEP, POST_STEP
	///
	typedef std::deque<assoc_steps_ele_t>			assoc_steps_t;
	
	///
	enum ETerminateStatus { STATUS_KEEP_RUNNING, STATUS_TERMINATE_TRUE, STATUS_TERMINATE_FALSE };

	///
	template<class T_ele>
	class name_comp {
	public:
		///
		name_comp(const std::string& name) : name_(name) {}
		///
		bool operator()(const T_ele& ele) { return ele.name == name_; }
	private:
		const std::string& name_;
	};	// end class name_comp

	typedef std::vector<double> step_times_t;

	///
	enum { NUM_STEP_TIME_STATS = 5 };

	// /////////////////////////////////////////////////////
	// Private data members

	// aggregate members

	state_ptr_t				state_;
	// ref_count_ptr<...> object for the aggragate AlgorithmState object.

	track_ptr_t				track_;
	// ref_count_ptr<...> object for the aggragate AlgorithmTrack object.

	// algorithm control etc.
	
	ERunningState			running_state_;
	// The state this Algorithm object is in:
	//
	// NOT_RUNNING					do_algorithm() has not been called.
	// RUNNING						do_algorithm() has been called.
	// RUNNING_BEING_CONFIGURED		do_algorithm() is active and begin_config_update() has been called
	//								but end_config_update() has not.

	size_t					first_k_;
	// The first iteration from state().k().
	
	size_t					max_iter_;
	// The maximum number of iterations that #this# will execute.

	double					max_run_time_;
	// The maximum amount of time the algorithm is allowed to execute.

	ETerminateStatus		terminate_status_;
	// Flag for if it is time to terminate do_algorithm().

	poss_type				next_step_poss_;
	// The next step possition that #this# will execute when control is returned to do_algorithm().

	const std::string*		next_step_name_;
	// The name of the next step that #this will execute when control is returned to do_algorithm().
	
	bool					do_step_next_called_;
	// Flag for if do_step_next(...) was called so that #do_algorithm(...)# can validate
	// that if a step object returned #false# from its #do_step(...)# operation that it
	// must also call #do_step_next(...)# to specify a step to jump to.

	poss_type				curr_step_poss_;
	// The current step being executed in do_algorithm(...).
	// If the current step being executed is changed during the imp_do_step(...) operation, then
	// imp_do_step(...) must adjust to this step.

	std::string				saved_curr_step_name_;
	// The name of the current step that is saved when begin_config_update() is called
	// so that curr_step_poss_ can be reset when end_config_update() is called.

	std::string				saved_next_step_name_;
	// The name of the next step to call so that when begin_config_update() is called
	// so that next_step_poss_ and next_step_name_ can be reset when end_config_update()
	// is called.

	bool					reconfigured_;
	// A flag that is set to true when a runtime configuration has been preformed.  It
	// is used to flag this event for imp_do_assoc_steps(...).

	// step and associated step object data structures

	steps_t					steps_;
	// Array of std::pair<ref_count_ptr<step_ptr_t>,std::string> objects.
	//
	// *steps_[step_poss].first returns the step object for step_poss = 1...steps_.size().
	// steps_[step_poss].second returns the name of the step for step_poss = 1...steps_.size().

	assoc_steps_t			assoc_steps_;
	// Array of two lists of std::pair<step_ptr_t,std::string> objects
	//
	// *(assoc_steps_[step_poss][PRE_STEP].begin() + pre_step_poss).first gives the pre step object.
	// (assoc_steps_[step_poss][PRE_STEP].begin() + pre_step_poss).second gives the name of the pre step
	// *(assoc_steps_[step_poss][POST_STEP].begin() + post_step_poss).first gives the post step object.
	// (assoc_steps_[step_poss][POST_STEP].begin() + post_step_poss).second gives the name of the post step

	bool algo_timing_;
	// If true each step will be timed.

	mutable step_times_t step_times_;
	// Array of step times ( size (max_iter() + 1 + NUM_STEP_TIME_STATS) * (num_steps() + 1) ).
	//  The time in sec. for step step_i (one based)
	// for iteration iter_k (zero based) is:
	// 	step_times_[ iter_k + (step_i - 1) * (max_iter() + 1 + NUM_STEP_TIME_STATS) ].
	// Therefore the times for each step are stored by column (consecutive elements)
	// so that statistics will be easy to compute at the end.
	// The last five elements after max_iter() for each step are reserved for:
	// * total time for the step
	// * average time for the step
	// * min step time
	// * max step time
	// * percentage for each step to the total.
	// The last column is for the total times for each iteration with the last five
	// elements being for the statistics for each iteration.	 

	mutable bool time_stats_computed_;
	// A flag for if the timing statistics have already been computed or not.
	
	mutable double total_time_;
	// Records the total computed time for the algorithm.

	// /////////////////////////////////////////////////////
	// Private member functions

	/// Validate a step_poss and throw a DoesNotExist exception if it does not.
	poss_type validate(poss_type step_poss, int past_end = 0) const;

	/// Validate an assoc_step_poss and throw a DoesNotExist exception if it does not.
	poss_type validate(const assoc_steps_ele_list_t& assoc_list, poss_type assoc_step_poss, int past_end = 0) const;


	/// Validate that #this# is in a specific running state.
	void validate_in_state(ERunningState running_state) const;

	/// Validate that #this# is not in a specific running state.
	void validate_not_in_state(ERunningState running_state) const;

	/// Validate that the step_poss in not the current step.
	void validate_not_curr_step(poss_type step_poss) const;

	/// Validate that the step_name in not the next step.
	void validate_not_next_step(const std::string& step_name) const;

	///
	/** Find a step given its name and throw a DoesNotExist exception if not found.
	  */
	steps_t::iterator step_itr(const std::string& step_name);

	///
	steps_t::const_iterator step_itr(const std::string& step_name) const;

	///
	/** Find a step given its name and throw a DoesNotExist exception if not found.
	  */
	steps_t::iterator step_itr_and_assert(const std::string& step_name);

	///
	steps_t::const_iterator step_itr_and_assert(const std::string& step_name) const;

	///
	/** Find a an associated step given its name and throw a DoesNotExist exception if not found.
	  */
	static assoc_steps_ele_list_t::iterator assoc_step_itr(assoc_steps_ele_list_t& assoc_list
		, const std::string& assoc_step_name);

	///
	static assoc_steps_ele_list_t::const_iterator assoc_step_itr(const assoc_steps_ele_list_t& assoc_list
		, const std::string& assoc_step_name);

	///
	bool imp_do_step(poss_type step_poss);

	///
	bool imp_do_assoc_steps(EAssocStepType type);

	///
	void imp_print_algorithm(std::ostream& out, bool print_steps) const;

	/// EAssocStepType -> EDoStepType
	EDoStepType do_step_type(EAssocStepType assoc_step_type);

};	// end class Algorithm

// //////////////////////////////////////////////////////////////////////////////////////////////////
// Inline member function definitions for Algorithm

// «std comp» members for state 

inline
void Algorithm::set_state(const state_ptr_t& state)
{	state_ = state; }

inline
Algorithm::state_ptr_t& Algorithm::get_state()
{	return state_; }

inline
const Algorithm::state_ptr_t& Algorithm::get_state() const
{	return state_; }

inline
AlgorithmState& Algorithm::state()
{	assert(state_.get()); return *state_; }

inline
const AlgorithmState& Algorithm::state() const
{	assert(state_.get()); return *state_; }

// «std comp» members for track 

inline
void Algorithm::set_track(const track_ptr_t& track)
{	track_ = track; }

inline
Algorithm::track_ptr_t& Algorithm::get_track()
{	return track_; }

inline
const Algorithm::track_ptr_t& Algorithm::get_track() const
{	return track_; }

inline
AlgorithmTrack& Algorithm::track()
{	assert(track_.get()); return *track_; }

inline
const AlgorithmTrack& Algorithm::track() const
{	assert(track_.get()); return *track_; }

// running state

inline
Algorithm::ERunningState Algorithm::running_state() const
{	return running_state_; }

// validate poss

inline
Algorithm::poss_type Algorithm::validate(poss_type step_poss, int past_end) const
{
	if( step_poss < 1 || steps_.size() + past_end < step_poss ) {
		std::ostringstream omsg;
		omsg	<< "Algorithm::validate(step_poss) : The step_poss = " << step_poss
				<< " is not in range of 1 to " << steps_.size() + past_end;
		throw DoesNotExist(omsg.str());
	}
	return step_poss;
}	

inline
Algorithm::poss_type Algorithm::validate(const assoc_steps_ele_list_t& assoc_list
	, poss_type assoc_step_poss, int past_end) const
{
	if( assoc_step_poss < 1 || assoc_list.size() + past_end < assoc_step_poss ) {
		std::ostringstream omsg;
		omsg	<< "Algorithm::validate(assoc_list,assoc_step_poss) : The assoc_step_poss = "
				<< assoc_step_poss << " is not in range of 1 to " << assoc_list.size() + past_end;
		throw DoesNotExist(omsg.str());
	}
	return assoc_step_poss;
}

// lookup iterator given name

inline
Algorithm::steps_t::iterator Algorithm::step_itr(const std::string& step_name)
{
	return std::find_if( steps_.begin() , steps_.end()
		, name_comp<steps_ele_t>(step_name) );
}

inline
Algorithm::steps_t::const_iterator Algorithm::step_itr(const std::string& step_name) const
{
	return std::find_if( steps_.begin() , steps_.end()
		, name_comp<steps_ele_t>(step_name) );
}

inline
Algorithm::assoc_steps_ele_list_t::iterator Algorithm::assoc_step_itr(
	assoc_steps_ele_list_t& assoc_list, const std::string& assoc_step_name)
{
	return std::find_if( assoc_list.begin() , assoc_list.end()
		, name_comp<assoc_steps_ele_list_ele_t>(assoc_step_name) );
}

inline
Algorithm::assoc_steps_ele_list_t::const_iterator Algorithm::assoc_step_itr(
	const assoc_steps_ele_list_t& assoc_list, const std::string& assoc_step_name)
{
	return std::find_if( assoc_list.begin() , assoc_list.end()
		, name_comp<assoc_steps_ele_list_ele_t>(assoc_step_name) );
}

inline
EDoStepType Algorithm::do_step_type(EAssocStepType assoc_step_type) {
	switch(assoc_step_type) {
		case PRE_STEP	: return DO_PRE_STEP;
		case POST_STEP	: return DO_POST_STEP;
	}
	assert(true);
	return DO_PRE_STEP;	// will never execute.
}

}	// end namespace GeneralIterationPack 

#endif // ALGORITHM_H