// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Exceptions.hpp

    \brief Exception handling macros
*/


#ifndef _ZOLTAN2_EXCEPTIONS_HPP_
#define _ZOLTAN2_EXCEPTIONS_HPP_

/*! \file Zoltan2_Exceptions.hpp

  Exception handling macros.  We throw 3 types of error:

   \li \c  std::runtime_error   for an error in input
   \li \c  std::bad_alloc       for failure to allocate memory
   \li \c  std::logic_error     for an apparent bug in the code

  The GLOBAL macros are for assertions that all processes in
  a communicator call.  They throw an error if any of the
  processes finds that the assertion fails.

  The LOCAL macros are local.  If the assertion fails, the
  process throws an error.
*/

#include <stdexcept>
#include <iostream>
#include <string>
#include <sstream>
#include <Zoltan2_Environment.hpp>
#include <Teuchos_CommHelpers.hpp>

#ifdef Z2_OMIT_ALL_ERROR_CHECKING

#define Z2_LOCAL_INPUT_ASSERTION(env, s, assertion, level) {}
#define Z2_GLOBAL_INPUT_ASSERTION(env, s, assertion, level) {}

#define Z2_LOCAL_BUG_ASSERTION(env, s, assertion, level) {}
#define Z2_GLOBAL_BUG_ASSERTION(env, s, assertion, level) {}

#define Z2_LOCAL_MEMORY_ASSERTION(env, nobj, assertion) {}
#define Z2_GLOBAL_MEMORY_ASSERTION(env, nobj , assertion) {}

#else

#define Z2_LOCAL_INPUT_ASSERTION(env, s, assertion, level) \
  if ((level <= (env).errorCheckLevel_) && !(assertion)){ \
    std::ostringstream assertion_msg; \
    assertion_msg<<(env).myRank_<<" "<<__FILE__<<","<<__LINE__; \
    assertion_msg<<", error: "<<s; \
    throw std::runtime_error(assertion_msg.str()); \
  }

#define Z2_LOCAL_BUG_ASSERTION(env, s, assertion, level) \
  if ((level <= (env).errorCheckLevel_) && !(assertion)){ \
    std::ostringstream assertion_msg; \
    assertion_msg<<(env).myRank_<<" "<<__FILE__<<","<<__LINE__; \
    assertion_msg<<", bug: "<<s; \
    throw std::logic_error(assertion_msg.str()); \
  }

/*! Memory assertions are always BASIC_ASSERTION level
 */
#define Z2_LOCAL_MEMORY_ASSERTION( env, nobj, assertion) { \
  if (!(assertion)){ \
    std::cerr<<(env).myRank_<<" "<<__FILE__<<","<<__LINE__<<","<<nobj<<" objects\n"; \
    throw std::bad_alloc(); \
  } \
}

#define Z2_GLOBAL_INPUT_ASSERTION( env, s, assertion, level) { \
  if (level <= (env).errorCheckLevel_) {  \
    int assertion_fail = 0, assertion_gfail=0;  \
    if (!(assertion)) assertion_fail = 1;  \
    Teuchos::reduceAll<int, int>(*(env).comm_, \
      Teuchos::REDUCE_MAX, 1, &assertion_fail, &assertion_gfail); \
    if (assertion_gfail > 0){  \
      std::ostringstream assertion_msg; \
      assertion_msg<<(env).myRank_<<" "<<__FILE__<<","<<__LINE__; \
      if (assertion_fail > 0) assertion_msg<<", error: "<<s; \
      else                    assertion_msg<<" exiting"; \
      throw std::runtime_error(assertion_msg.str()); \
    } \
  } \
}

#define Z2_GLOBAL_BUG_ASSERTION( env, s, assertion, level) { \
  if (level <= (env).errorCheckLevel_) {  \
    int assertion_fail = 0, assertion_gfail=0;  \
    if (!(assertion)) assertion_fail = 1;  \
    Teuchos::reduceAll<int, int>(*(env).comm_, \
      Teuchos::REDUCE_MAX, 1, &assertion_fail, &assertion_gfail); \
    if (assertion_gfail > 0){  \
      std::ostringstream assertion_msg; \
      assertion_msg<<(env).myRank_<<" "<<__FILE__<<","<<__LINE__; \
      if (assertion_fail > 0) assertion_msg<<", bug: "<<s; \
      else                    assertion_msg<<" exiting"; \
      throw std::logic_error(assertion_msg.str()); \
    } \
  } \
}

/*! Memory assertions are always BASIC_ASSERTION level
 */
#define Z2_GLOBAL_MEMORY_ASSERTION( env, nobj, assertion) {\
  int assertion_fail = 0, assertion_gfail=0;  \
  if (!(assertion)) assertion_fail = 1;  \
  Teuchos::reduceAll<int, int>(*(env).comm_, \
      Teuchos::REDUCE_MAX, 1, &assertion_fail, &assertion_gfail); \
  if (assertion_gfail > 0){  \
    int me = (env).myRank_; \
    if (assertion_fail > 0) \
      std::cerr<<me<<" "<<__FILE__<<","<<__LINE__<<","<<nobj<<" objects\n"; \
    else \
      std::cerr<<me<<" "<<__FILE__<<","<<__LINE__<<", exiting"; \
    throw std::bad_alloc(); \
  } \
}

#endif

/*! Throw an error returned from outside the Zoltan2 library.
 */
#define Z2_THROW_OUTSIDE_ERROR(env, e) { \
  std::cerr<<(env).myRank_<<" "<<__FILE__<<","<<__LINE__<<","<<e.what()<<std::endl; \
  throw e; \
}

//! Forward an exception back through call stack.
#define Z2_FORWARD_EXCEPTIONS \
  catch (std::runtime_error &e) { throw e; } \
  catch (std::logic_error   &e) { throw e; } \
  catch (std::bad_alloc     &e) { throw e; } 
   
#endif
