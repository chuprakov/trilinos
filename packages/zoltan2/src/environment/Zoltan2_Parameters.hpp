// @HEADER
// ***********************************************************************
//                Copyright message goes here.   TODO
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Parameters.hpp
    \brief Defines Parameter related enumerators, declares functions.
*/

#ifndef _ZOLTAN2_PARAMETERS_HPP_
#define _ZOLTAN2_PARAMETERS_HPP_

#include <Zoltan2_Standards.hpp>

namespace Zoltan2{

////////////////////////////////////////////////////////////////////
// Parameter-related namespace methods 

/*! \brief Create a ParameterList suitable for validating another list.
 *
 *   \param  plIn   a ParameterList that has not been validated yet.
 *   \param  plOut  on return, plOut has exactly the same parameters are
 *                      plIn, but it also has a validator for each parameter.
 *
 *   The validating ParameterList is required to verify that the
 *   input ParameterList is valid.
 */
void createValidatorList(
   const Teuchos::ParameterList &plIn, Teuchos::ParameterList &plOut);

/*! \brief Print out the ParameterList's built in documenation.
 *
 *   \param  pl      a ParameterList
 *   \param  os      the output stream to which the documentation should go
 *   \param listNames  leave this unset - it is set in recursive calls
 */
void printListDocumentation( const Teuchos::ParameterList &pl, ostream &os,
  string listNames=string(""));

////////////////////////////////////////////////////////////////////
// Parameter-related enumerated types.
//
//  If you change these enumerators, change their documentation
//  in Zoltan2_Parameters.cpp.
//

/*! \brief Level of error checking or assertions desired.
 *
 *  Each assertion in the code must have a level. Tests for
 *  logic errors should always be level DEBUG_MODE_ASSERTION.
 *  Quick tests are BASIC, longer tests for input errors are
 *  COMPLEX.
 *
 *  The user sets the assertion level with the parameter \c error_check_level.
 *
 *  If compiled with \b Z2_OMIT_ALL_ERROR_CHECKING, error checks don't happen.
 */

enum AssertionLevel {
  NO_ASSERTIONS,   /*!< \brief no assertion checks will be done */
  BASIC_ASSERTION, /*!< \brief fast typical checks for valid arguments */
  COMPLEX_ASSERTION,  /*!< \brief more involved, like validate a graph */
  DEBUG_MODE_ASSERTION, /*!< \brief checks for logic errors */
  NUM_ASSERTION_LEVELS
};

/*! \brief The amount of debugging or status output to print.
 *
 *  Each debug/status message must have an output level.  The
 *  user specfies the level desired with the \c debug_level parameter.
 *
 *  If Zoltan2 is compiled with \b Z2_OMIT_ALL_STATUS_MESSAGES, no
 *  messages will be processed. 
 */
 
enum MessageOutputLevel {
  NO_STATUS,         /*!< \brief don't display status/debug messages */
  BASIC_STATUS,      /*!< \brief the status at each high level step */
  DETAILED_STATUS,   /*!< \brief sub-steps, each method's entry and exit */
  VERBOSE_DETAILED_STATUS, /*!< \brief include more detail about sub-steps */
  NUM_STATUS_OUTPUT_LEVELS
};

/*! \brief The type of timers which should be active.
 *
 *  If timing is requested, \c timer_type can specify MACRO timing
 *  or MICRO timing.  For example, a MACRO timer would time an
 *  algorithm.  A MICRO timer would time steps in the algorithm.
 *  You can also ask for BOTH, but be aware that your MACRO timer
 *  is also timing the MICRO timers as well as the algorithm.
 *  
 *  If Zoltan2 is compiled with \b Z2_OMIT_ALL_PROFILING
 *  timing messages are ignored.
 */
 
enum TimerType {
  NO_TIMERS,    /*!< \brief No timing data will be collected (the default). */
  MACRO_TIMERS, /*!< \brief Time an algorithm (or other entity) as a whole. */
  MICRO_TIMERS, /*!< \brief Time the substeps of an entity. */
  BOTH_TIMERS,  /*!< \brief Run both MACRO and MICRO timers. */
  TEST_TIMERS,  /*!< \brief Timers added while testing, removed later. */
  NUM_TIMING_OPTIONS
};

/*! \brief Output stream types.
 */
 
enum OSType {
  COUT_STREAM,    /*!< \brief std::cout */
  CERR_STREAM,    /*!< \brief std::cerr */
  NULL_STREAM,    /*!< \brief /dev/null: do actions but don't output results */
  NUM_OUTPUT_STREAMS
};

/*!\brief  Enumerator used in code for multicriteria norm choice.
 */
enum multiCriteriaNorm{
  normMinimizeTotalWeight,   /*!< 1-norm = Manhattan norm */
  normBalanceTotalMaximum,   /*!< 2-norm = sqrt of sum of squares */
  normMinimizeMaximumWeight, /*!< inf-norm = maximum norm */
  normNumNorms
};

}  // end of namespace Zoltan2

#endif
