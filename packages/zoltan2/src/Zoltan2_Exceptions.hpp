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

  The bad_alloc exceptions are thrown in Zoltan2_Memory.hpp.
*/

#include <stdexcept>
#include <iostream>
#include <string>
#include <Zoltan2_Environment.hpp>
#include <Teuchos_CommHelpers.hpp>

#ifdef Z2_OMIT_ALL_ERROR_CHECKING

#define Z2_LOCAL_INPUT_ASSERTION(comm, env, s, assertion, level) {}
#define Z2_LOCAL_BUG_ASSERTION(comm,  env, s, assertion, level) {}
#define Z2_LOCAL_MEMORY_ASSERTION(comm,  env, requestSize, assertion) {}
#define Z2_GLOBAL_INPUT_ASSERTION( comm, env, s, assertion, level) {}
#define Z2_GLOBAL_BUG_ASSERTION( comm, env, s, assertion, level) {}
#define Z2_GLOBAL_MEMORY_ASSERTION( comm, env, requestSize, assertion) {}

#else

#define Z2_LOCAL_INPUT_ASSERTION(comm, env, s, assertion, level) { \
  if (level <= (env).errorCheckLevel_) { \
    if (!(assertion)){ \
      std::ostringstream oss; \
      oss << __FILE__ << ", " << __LINE__ << ", " << s << std::endl; \
      (env).dbg_->error(oss.str()); \
      throw std::runtime_error(oss.str()); \
    } \
  } \
}

#define Z2_LOCAL_BUG_ASSERTION( comm, env, s, assertion, level) { \
  if (level <= (env).errorCheckLevel_) { \
    if (!(assertion)){ \
      std::ostringstream oss; \
      oss << __FILE__ << ", " << __LINE__ << ", " << s << std::endl; \
      (env).dbg_->error(oss.str()); \
      throw std::logic_error(oss.str()); \
    } \
  } \
}

/*! We always check for success of memory allocation, regardless of ERROR_CHECK_LEVEL.
 */

#define Z2_LOCAL_MEMORY_ASSERTION( comm, env, requestSize, assertion) { \
  if (!(assertion)){ \
    std::ostringstream _msg; \
    _msg << __FILE__ << ", " << __LINE__ << ", size " << requestSize << std::endl; \
    (env).dbg_->error(_msg.str()); \
    throw std::bad_alloc(); \
  } \
}

#define Z2_GLOBAL_INPUT_ASSERTION( comm, env, s, assertion, level) { \
  if (level <= (env).errorCheckLevel_) {  \
    int fail = 0, gfail=0;  \
    if (!(assertion)) fail = 1;  \
    Teuchos::reduceAll<int, int>(comm, Teuchos::REDUCE_MAX, 1, &fail, &gfail); \
    if (gfail > 0){  \
      std::ostringstream _msg; \
      if (fail > 0){ \
        _msg << __FILE__ << ", " << __LINE__ << ", " << s << std::endl; \
        (env).dbg_->error(_msg.str()); \
      } \
      throw std::runtime_error(_msg.str()); \
    } \
  } \
}

#define Z2_GLOBAL_BUG_ASSERTION( comm, env, s, assertion, level) { \
  if (level <= (env).errorCheckLevel_) {  \
    int fail = 0, gfail=0;  \
    if (!(assertion)) fail = 1;  \
    Teuchos::reduceAll<int, int>(comm, Teuchos::REDUCE_MAX, 1, &fail, &gfail); \
    if (gfail > 0){  \
      std::ostringstream _msg; \
      if (fail > 0){ \
        _msg <<  __FILE__ << ", " << __LINE__ << ", " << s << std::endl; \
        (env).dbg_->error(_msg.str()); \
      } \
      throw std::logic_error(_msg.str()); \
    } \
  } \
}

/*! We always check for success of memory allocation, regardless of ERROR_CHECK_LEVEL.
 */

#define Z2_GLOBAL_MEMORY_ASSERTION( comm, env, requestSize, assertion) {\
  int fail = 0, gfail=0;  \
  if (!(assertion)) fail = 1;  \
  Teuchos::reduceAll<int, int>(comm, Teuchos::REDUCE_MAX, 1, &fail, &gfail); \
  if (gfail > 0){  \
    if (fail > 0){ \
      std::ostringstream _msg; \
      _msg <<  __FILE__ << ", " << __LINE__ << ", size " << requestSize << std::endl; \
      (env).dbg_->error(_msg.str()); \
    } \
    throw std::bad_alloc(); \
  } \
}

#endif

/*! Throw an error returned from outside the Zoltan2 library.
 */
#define Z2_THROW_OUTSIDE_ERROR(env, e) { \
   std::ostringstream oss; \
   oss << __FILE__ << ":" << __LINE__ << " " << e.what() << std::endl; \
   (env).dbg_->error(oss.str()); \
  throw e; \
}

//! Forward an exception back through call stack.
#define Z2_FORWARD_EXCEPTIONS \
  catch (std::runtime_error &e) { throw e; } \
  catch (std::logic_error   &e) { throw e; } \
  catch (std::bad_alloc     &e) { throw e; } 
   
#endif
