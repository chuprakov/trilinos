// //////////////////////////////////////////////////////////
// rSQPppSolver.cpp
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

#include <fstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "ReducedSpaceSQPPack/Configurations/rSQPppSolver.h"
#include "ReducedSpaceSQPPack/Configurations/MamaJama/rSQPAlgo_ConfigMamaJama.h"
#include "ReducedSpaceSQPPack/include/rSQPSolverClientInterfaceSetOptions.h"
#include "ReducedSpaceSQPPack/include/rSQPAlgoClientInterface.h"
#include "ReducedSpaceSQPPack/include/rSQPAlgoContainer.h"
#include "ReducedSpaceSQPPack/include/rSQPState.h"
#include "ReducedSpaceSQPPack/include/std/rSQPTrackSummaryStd.h"
#include "ReducedSpaceSQPPack/include/std/rSQPTrackConsoleStd.h"
#include "ReducedSpaceSQPPack/include/std/rSQPTrackStatsStd.h"
#include "GeneralIterationPack/include/AlgorithmTrackComposite.h"
#include "NLPInterfacePack/include/NLPFirstOrderInfo.h"
#include "NLPInterfacePack/include/NLPFirstOrderDirect.h"
#include "NLPInterfacePack/test/test_nlp_first_order_info.h"
#include "NLPInterfacePack/test/test_nlp_first_order_direct.h"
#include "OptionsFromStream.h"
#include "stpwatch.h"
#include "StringToIntMap.h"
#include "StringToBool.h"
#include "WorkspacePack.h"
#include "ThrowException.h"

namespace ReducedSpaceSQPPack {

// Initialization and algorithm configuration

rSQPppSolver::rSQPppSolver()
	:reconfig_solver_(true)
	,workspace_MB_(-1.0)
	,obj_scale_(1.0)
	,test_nlp_(true)
	,print_algo_(true)
	,algo_timing_(true)
	,generate_stats_file_(false)
	,print_opt_grp_not_accessed_(true)
{}

void rSQPppSolver::set_nlp(const nlp_ptr_t& nlp)
{
	nlp_ = nlp;
	reconfig_solver_ = true;
}
	
const rSQPppSolver::nlp_ptr_t&
rSQPppSolver::get_nlp() const
{
	return nlp_;
}

void rSQPppSolver::set_track(const track_ptr_t& track)
{
	track_ = track;
	solver_.set_track(NULL); // Force the track objects to be rebuilt and added!
}
	
const rSQPppSolver::track_ptr_t&
rSQPppSolver::get_track() const
{
	return track_;
}
	
void rSQPppSolver::set_config( const config_ptr_t& config )
{
	config_ = config;
	solver_.set_config(NULL); // Must unset the config object.
	reconfig_solver_ = true;
}

const rSQPppSolver::config_ptr_t&
rSQPppSolver::get_config() const
{
	return config_;
}

void rSQPppSolver::set_options( const options_ptr_t& options )
{
	options_ = options;                  // Must totally free all of the references we
	const config_ptr_t                   // have to the current options.  That includes
		&config = solver_.get_config();  // removing the options object for the configuration
	if(config.get())                     // object.
		config->set_options(NULL);       // ...
	reconfig_solver_ = true;
}

const rSQPppSolver::options_ptr_t&
rSQPppSolver::get_options() const
{
	return options_;
}

void rSQPppSolver::set_console_out( const ostream_ptr_t& console_out )
{
	console_out_      = console_out;
	console_out_used_ = NULL;  // Remove every reference to this ostream object!
	solver_.set_track(NULL);
}

const rSQPppSolver::ostream_ptr_t&
rSQPppSolver::get_console_out() const
{
	return console_out_;
}

void rSQPppSolver::set_summary_out( const ostream_ptr_t& summary_out )
{
	summary_out_      = summary_out;
	summary_out_used_ = NULL;
	solver_.set_track(NULL);     // Remove every reference to this ostream object!
}
	
const rSQPppSolver::ostream_ptr_t&
rSQPppSolver::get_summary_out() const
{
	return summary_out_;
}

void rSQPppSolver::set_journal_out( const ostream_ptr_t& journal_out )
{
	journal_out_      = journal_out;
	journal_out_used_ = NULL;
	solver_.set_track(NULL);     // Remove every reference to this ostream object!
}
	
const rSQPppSolver::ostream_ptr_t&
rSQPppSolver::get_journal_out() const
{
	return journal_out_;
}

void rSQPppSolver::set_algo_out( const ostream_ptr_t& algo_out )
{
	algo_out_      = algo_out;
	algo_out_used_ = NULL;
}
	
const rSQPppSolver::ostream_ptr_t&
rSQPppSolver::get_algo_out() const
{
	return algo_out_;
}

// Solve the NLP

rSQPppSolver::ESolutionStatus rSQPppSolver::solve_nlp() const
{
	using std::endl;
	using std::setw;
	using StopWatchPack::stopwatch;
	namespace rcp = ReferenceCountingPack;
	using rcp::ref_count_ptr;

	update_solver();

	ESolutionStatus solve_return    = SOLVE_RETURN_EXCEPTION;
	bool            threw_exception = false;

	solver_.set_algo_timing(algo_timing_);

	//
	// Test the nlp if needed
	//
	
	if(test_nlp_) {
		
		// Set the scaling for the NLP
		nlp_->scale_f(obj_scale_);

		//
		// Test the nlp if needed
		//
		
		if( test_nlp_ ) {
			const char msg1[] = "\ntest_nlp = true: Testing the NLP!\n";
			*summary_out_used_ << msg1;
			*journal_out_used_ << msg1;
			if(NLPFirstOrderInfo* nlp_foi = dynamic_cast<NLPFirstOrderInfo*>(nlp_.get())) {
				const char msg[] = "\nTesting the supported NLPFirstOrderInfo interface ...\n";
				*summary_out_used_ << msg;
				*journal_out_used_ << msg;
				const bool
					result = NLPInterfacePack::test_nlp_first_order_info(nlp_foi,options_used_.get(),journal_out_used_.get());
				if(!result) {
					const char msg[] = "\nTest failed! exiting!\n";
					*summary_out_used_ << msg;
					*journal_out_used_ << msg;
					solve_return = SOLVE_RETURN_NLP_TEST_FAILED;
					return solve_return;
				}
			}
			if(NLPFirstOrderDirect* nlp_fod = dynamic_cast<NLPFirstOrderDirect*>(nlp_.get())) {
				const char msg[] = "\nTesting the supported NLPFirstOrderDirect interface ...\n";
				*summary_out_used_ << msg;
				*journal_out_used_ << msg;
				const bool
					result = NLPInterfacePack::test_nlp_first_order_direct(nlp_fod,options_used_.get(),journal_out_used_.get());
				if(!result) {
					const char msg[] = "\nTest failed! exiting!\n";
					*summary_out_used_ << msg;
					*journal_out_used_ << msg;
					solve_return = SOLVE_RETURN_NLP_TEST_FAILED;
					return solve_return;
				}
			}
			const char msg2[] = "\n... end testing of nlp\n";
			*summary_out_used_ << msg2;
			*journal_out_used_ << msg2;
		}

	}
	
	//
	// Solve the NLP
	//

	stopwatch timer;

	typedef ReducedSpaceSQPPack::rSQPSolverClientInterface    solver_interface_t;
	solver_interface_t::EFindMinReturn                        r_find_min;

	*journal_out_used_
		<< "\n*** Starting rSQP iterations ...\n\n";	
	
	try {
		timer.start();
		r_find_min = solver_.find_min();
	}
	catch(const std::exception& excpt) {
		*summary_out_used_ << "\nCaught a std::exception: " << excpt.what() << endl;
		*journal_out_used_ << "\nCaught a std::exception: " << excpt.what() << endl;
		std::cerr << "\nCaught a std::exception: " << excpt.what() << endl;
		threw_exception = true;
	}
	catch(...) {
		*summary_out_used_ << "\nCaught an unknown exception\n";
		*journal_out_used_ << "\nCaught an unknown exception\n";
		std::cerr << "\nCaught an unknown exception\n";
		threw_exception = true;
	}
	timer.stop();
	
	if(threw_exception) {
		*summary_out_used_	<< "\n\n****************************\n"
							<< "**** Threw an exception ****\n";
		std::cerr << "Threw an exception!\n";
		solve_return = SOLVE_RETURN_EXCEPTION;
	}
	else {
		switch( r_find_min ) {
		    case solver_interface_t::SOLUTION_FOUND: {
				*summary_out_used_	<< "\n\n************************\n"
									<< "**** Solution Found ****\n";
				std::cerr << "Solution Found!\n";
				solve_return = SOLVE_RETURN_SOLVED;
				break;
			}
		    case solver_interface_t::MAX_ITER_EXCEEDED: {
				*summary_out_used_	<< "\n\n**********************************************\n"
									<< "**** Maximun number of iteration exceeded ****\n";
				std::cerr << "Maximun number of iteration exceeded!\n";
				solve_return = SOLVE_RETURN_MAX_ITER;
				break;
			}
		    case solver_interface_t::MAX_RUN_TIME_EXCEEDED: {
				*summary_out_used_	<< "\n\n**********************************\n"
									<< "**** Maximun runtime exceeded ****\n";
				std::cerr << "Maximun runtime exceeded!\n";
				solve_return = SOLVE_RETURN_MAX_RUN_TIME;
				break;
			}
		    case solver_interface_t::ALGORITHMIC_ERROR: {
				*summary_out_used_	<< "\n\n*********************************************\n"
									<< "**** Some error occured in the algorithm ****\n";
				std::cerr << "Some algorithmic error occured!\n";
				solve_return = SOLVE_RETURN_EXCEPTION;
				break;
			}
		}
	}
	
	*summary_out_used_	<< "\n  total time = " << timer.read() << " sec.\n";

	if( solver_.algo_timing() )
		solver_.print_algorithm_times( *summary_out_used_ );
	
	// Print workspace usage statistics
	*summary_out_used_
		<< "\n*** Statistics for autmatic array workspace:"
		<< "\nNumber of megabytes of preallocated workspace                = "
		<< workspace_MB_
		<< "\nNumber of allocations using preallocated workspace           = "
		<< WorkspacePack::default_workspace_store->num_static_allocations()
		<< "\nNumber of dynamic allocations beyond preallocated workspace  = "
		<< WorkspacePack::default_workspace_store->num_dyn_allocations();
	
	// Print which options groups were not read
	if(print_opt_grp_not_accessed_) {
		*algo_out_used_
			<<	"\n***************************************************************\n"
			"Warning, the following options groups where not accessed.\n"
			"An options group may not be accessed if it is not looked for\n"
			"or if an \"optional\" options group was looked from and the user\n"
			"spelled it incorrectly:\n\n";
		options_used_->print_unaccessed_options_groups(*algo_out_used_);
	}

	algo_out_used_->flush();
	summary_out_used_->flush();
	journal_out_used_->flush();
	if(console_out_used_.get()) console_out_used_->flush();
	
	return solve_return;
	
}

// Get the underlying solver object

rSQPSolverClientInterface& rSQPppSolver::get_solver()
{
	update_solver();
	return solver_;
}

const rSQPSolverClientInterface& rSQPppSolver::get_solver() const
{
	update_solver();
	return solver_;
}

// private

void rSQPppSolver::update_solver() const
{
	using std::endl;
	using std::setw;
	using StopWatchPack::stopwatch;
	namespace rcp = ReferenceCountingPack;
	using rcp::ref_count_ptr;
	namespace ofsp = OptionsFromStreamPack;
	using ofsp::OptionsFromStream;
	using ofsp::StringToIntMap;
	using ofsp::StringToBool;

	//
	// Validate the input
	//

	THROW_EXCEPTION(
		nlp_.get() == NULL, std::logic_error
		,"rSQPSolver::update_solver() : Error, this->get_nlp().get() can not be NULL!" );

	//
	// Get the output streams if needed
	//
	
	if( summary_out_used_.get()==NULL ) {
		if( summary_out_.get() == NULL )
			summary_out_used_ = rcp::rcp(new std::ofstream("rSQPppSummary.out"));
		else
			summary_out_used_ = summary_out_;
	}
	if( journal_out_used_.get() == NULL ) {
		if( journal_out_.get() == NULL )
			journal_out_used_ = rcp::rcp(new std::ofstream("rSQPppJournal.out"));
		else
			journal_out_used_ = journal_out_;
	}
	if( algo_out_used_.get() == NULL ) {
		if( algo_out_.get() == NULL )
			algo_out_used_ = rcp::rcp(new std::ofstream("rSQPppAlgo.out"));
		else
			algo_out_used_ = algo_out_;
	}
	if( console_out_used_.get() == NULL ) {
		if( console_out_.get() != NULL )
			console_out_used_ = console_out_;
		else
			console_out_used_ = NULL;
	}

	//
	// Get the options (or lack of)
	//
	
	if( options_used_.get() == NULL ) {
		if( options_.get() == NULL ) {
			std::ifstream options_in("rSQPpp.opt");
			options_used_ = ( options_in
							  ? rcp::rcp<options_ptr_t::element_type>(new OptionsFromStream(options_in))
							  : rcp::rcp<options_ptr_t::element_type>(NULL)
				);
		}
		else
			options_used_ = options_;
	}

	//
	// Read in some options for "rSQPppSolver" if needed
	//
	
	if( reconfig_solver_ || solver_.get_track().get() == NULL ) {

		options_used_->reset_unaccessed_options_groups();

		const std::string optgrp_name = "rSQPppSolver";
		OptionsFromStream::options_group_t optgrp = options_used_->options_group( optgrp_name );
		if( OptionsFromStream::options_group_exists( optgrp ) ) {
			
			const int num_opt = 7;
			enum EOptions {
				WORKSPACE_MB
				,OBJ_SCALE
				,TEST_NLP
				,PRINT_ALGO
				,ALGO_TIMING
				,GENERATE_STATS_FILE
				,PRINT_OPT_GRP_NOT_ACCESSED
			};
			const char* SOptions[num_opt] = {
				"workspace_MB"
				,"obj_scale"
				,"test_nlp"
				,"print_algo"
				,"algo_timing"
				,"generate_stats_file"
				,"print_opt_grp_not_accessed"
			};
			StringToIntMap	opt_map( optgrp_name, num_opt, SOptions );
			
			OptionsFromStream::options_group_t::const_iterator itr = optgrp.begin();
			for( ; itr != optgrp.end(); ++itr ) {
				switch( (EOptions)opt_map( ofsp::option_name(itr) ) ) {
					case WORKSPACE_MB:
						workspace_MB_ = ::atof( ofsp::option_value(itr).c_str() );
						break;
					case OBJ_SCALE:
						obj_scale_ = ::atof( ofsp::option_value(itr).c_str() );
						break;
					case TEST_NLP:
						test_nlp_ = StringToBool( "test_nlp", ofsp::option_value(itr).c_str() );
						break;
					case PRINT_ALGO:
						print_algo_ = StringToBool( "print_algo", ofsp::option_value(itr).c_str() );
						break;
					case ALGO_TIMING:
						algo_timing_ = StringToBool( "algo_timing", ofsp::option_value(itr).c_str() );
						break;
					case GENERATE_STATS_FILE:
						generate_stats_file_ = StringToBool( "generate_stats_file", ofsp::option_value(itr).c_str() );
						break;
					case PRINT_OPT_GRP_NOT_ACCESSED:
						print_opt_grp_not_accessed_ = StringToBool( "algo_timing", ofsp::option_value(itr).c_str() );
						break;
					default:
						assert(0);	// this would be a local programming error only.
				}
			}
		}
		else {
			*algo_out_used_
				<< "\nWarning!  The options group \'rSQPppSolver\' was not found.\n"
				"Using a default set of options ...\n";
		}
	}

	bool did_reconfigured_solver = false;
	if(reconfig_solver_) {
		did_reconfigured_solver = true;

		//
		// Print the headers for the output files
		//
		
		int w = 15;
		int prec = 8;
		summary_out_used_->precision(prec);
		*journal_out_used_ << std::setprecision(prec) << std::scientific;

		*algo_out_used_
			<< "\n********************************************************************"
			<< "\n*** Algorithm information output                                 ***"
			<< "\n***                                                              ***"
			<< "\n*** Below, information about how the the rSQP++ algorithm is     ***"
			<< "\n*** setup is given and is followed by detailed printouts of the  ***"
			<< "\n*** contents of the algorithm state object (i.e. iteration       ***"
			<< "\n*** quantities) and the algorithm description printout           ***"
			<< "\n*** (if the option rSQPppSolver::print_algo = true is set).      ***"
			<< "\n********************************************************************\n";
		*summary_out_used_
			<< "\n********************************************************************"
			<< "\n*** Algorithm iteration summary output                           ***"
			<< "\n***                                                              ***"
			<< "\n*** Below, a summary table of the SQP iterations is given as     ***"
			<< "\n*** well as a table of the CPU times for each step (if the       ***"
			<< "\n*** option rSQPppSolver::algo_timing = true is set).             ***"
			<< "\n********************************************************************\n";
		*journal_out_used_
			<< "\n********************************************************************"
			<< "\n*** Algorithm iteration detailed journal output                  ***"
			<< "\n***                                                              ***"
			<< "\n*** Below, detailed information about the SQP algorithm is given ***"
			<< "\n*** while it is running.  The amount of information that is      ***"
			<< "\n*** produced can be specified using the option                   ***"
			<< "\n*** rSQPSolverClientInterface::journal_output_level (the default ***"
			<< "\n*** is PRINT_NOTHING and produces no output                      ***"
			<< "\n********************************************************************\n";

		// Echo options.
		*summary_out_used_ << "\n*** Echoing input options ...\n";
		options_used_->print_options( *summary_out_used_ );
		summary_out_used_->flush();
		*algo_out_used_ << "\n*** Echoing input options ...\n";
		options_used_->print_options( *algo_out_used_ );
		algo_out_used_->flush();
		*journal_out_used_ << "\n*** Echoing input options ...\n";
		options_used_->print_options( *journal_out_used_ );
		algo_out_used_->flush();

		//
		// Allocate the workspace
		//
		
		nlp_->initialize();
		const int default_ws_scale = 10;
		if( workspace_MB_ < 0.0 ) {
			workspace_MB_ = nlp_->n() * default_ws_scale * 1e-6 * sizeof(value_type);
			*algo_out_used_
				<< "\nworkspace_MB < 0.0:\n"
				<< "Setting workspace_MB = n * default_ws_scale * 1e-6 * sizeof(value_type) = "
				<< nlp_->n() << " * " << default_ws_scale << " * 1e-6 * " << sizeof(value_type)
				<< " = " << workspace_MB_ << " MB\n";
		}
		*summary_out_used_
			<< "\nAllocating workspace_MB = " << workspace_MB_ << " megabytes of temporary "
			"workspace for autmatic arrays only ...\n";
		WorkspacePack::default_workspace_store
			= rcp::rcp(new WorkspacePack::WorkspaceStoreInitializeable(1e+6*workspace_MB_));

		//
		// Reconfigure the algorithm
		//
		
		// Get and set up the configuration object
		config_ptr_t _config;
		if(config_.get() == NULL)
			_config = rcp::rcp(new rSQPAlgo_ConfigMamaJama());
		else
			_config = config_;
		_config->set_options( options_used_ );

		{
			std::ostringstream msg;
			msg << "\n*** Setting up to run rSQP++ on the NLP using a "
				<< "configuration object of type \'" << typeid(*_config).name() << "\' ...\n";
			*algo_out_used_
				<< msg.str();
			*summary_out_used_
				<< msg.str();
			*journal_out_used_
				<< msg.str();
		}
		
		// Set up the solver
	
		solver_.set_nlp(nlp_);           // Set the NLP object
		solver_.set_config(_config);     // Set the configuration object
		
		// Set the client interface options
		rSQPSolverClientInterfaceSetOptions
			solver_options_setter( &solver_ );
		solver_options_setter.set_options( *options_used_ );

		did_reconfigured_solver = true;
		reconfig_solver_    = false;
	}

	//
	// Set up the track objects if needed
	//

	if( solver_.get_track().get() == NULL ) {
		ref_count_ptr<AlgorithmTrackComposite>
			composite_track = rcp::rcp(new AlgorithmTrackComposite(journal_out_used_));
		composite_track->tracks().push_back(
			rcp::rcp(new rSQPTrackSummaryStd(summary_out_used_,journal_out_used_)) );
		if( console_out_used_.get() ) {
			composite_track->tracks().push_back(
				rcp::rcp(new rSQPTrackConsoleStd(console_out_used_,journal_out_used_)) );
		}
		if(generate_stats_file_) {
			ostream_ptr_t
				stats_out = rcp::rcp(new std::ofstream("rSQPppStats.out"));
			assert( !stats_out->eof() );
			composite_track->tracks().push_back(
				rcp::rcp(new rSQPTrackStatsStd(stats_out,stats_out)) );
		}
		if( track_.get() ) {
			track_->set_journal_out(journal_out_used_);
			composite_track->tracks().push_back( track_ );
		}
		solver_.set_track( composite_track );
	}

	//
	// Configure and print the algorithm if needed
	//

	if(did_reconfigured_solver) {
		solver_.configure_algorithm(algo_out_used_.get());
		if(print_algo_)
			solver_.print_algorithm(*algo_out_used_);
	}
	
}

} // end namespace ReducedSpaceSQPPack
