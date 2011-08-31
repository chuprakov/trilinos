// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Environment.hpp

    \brief The declarations for the Environment object.
*/


#ifndef _ZOLTAN2_ENVIRONMENT_DECL_HPP_
#define _ZOLTAN2_ENVIRONMENT_DECL_HPP_

/*! \file Zoltan2_Environment_decl.hpp
  
  \brief Declares the Zoltan2::Environment class.

*/

#include <Teuchos_ParameterList.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_DebugManager.hpp>

namespace Zoltan2 {

/*! Zoltan2::Environment
    \brief The problem parameters, library configuration, and other information.

  This is object is passed to almost every method in the library.  We may want
  to have a memory manager here as well.

  TODO: add teuchos validators, and all the other parameters
  TODO: Do we need to template on AppLID and AppGID for fixed vertices?
*/

// Definitions for ERROR_CHECK_LEVEL parameter.

/*!  We should always check basic assertions.
*/
#define Z2_BASIC_ASSERTION      0

/*!  Extra, more expensive level of checking.
 *
 * A parameter will state whether "extra" checking should be
 * done.  An example of extra checking is checking that an
 * input graph is valid.
 */
#define Z2_COMPLEX_ASSERTION    1

/*!  Even more extensive checking.
 *
 * This is extra checking we would do when debugging
 * a problem.
 *
 */
#define Z2_DEBUG_MODE_ASSERTION  2

#define Z2_MAX_CHECK_LEVEL Z2_DEBUG_MODE_ASSERTION

class Environment{

private:
  /*! The problem parameters
   */
  Teuchos::ParameterList _params;
  
  /*! The library configuration 
   */
  Teuchos::ParameterList _config;
  
#if 0
  /*! The node description is not yet implemented.
   */
  Kokkos::CUDANodeMemoryModel     _gpuNode;
  Kokkos::StandardNodeMemoryModel _standardNode;
  
  /*! The machine model is not yet implemented.  It will
      not necessarily be implemented as a ParameterList.
   */
  Teuchos::ParameterList _machine;
#endif

public:
  /*! Constructor 
      Because parameter lists are small, we save a
      copy of them instead of requiring an RCP.
   */
  Environment(Teuchos::ParameterList &prob, Teuchos::ParameterList &config);

  /*! Constructor
   */
  Environment();

  /*! Destructor */
  virtual ~Environment();

  /*! Copy Constructor */
  Environment(const Environment &env);

  /*! Assignment operator */
  Environment &operator=(const Environment &env);

  /*! Set or reset the problem parameters*/
  void setProblemParameters(Teuchos::ParameterList &problemParams);

  /*! Set or reset the library configuration */
  void setLibraryConfiguration(Teuchos::ParameterList &libraryConfig);

  /*! Get a reference to a read-only copy of the problem parameters. */
  const Teuchos::ParameterList &getProblemParameters() const;

  /*! Get a reference to a read-only copy of the library configuration. */
  const Teuchos::ParameterList &getLibraryConfiguration() const;

  /*! The debug manager, used by debug statements */
  Zoltan2::DebugManager _dbg;

  /*! The level of checking to do at runtime.  See Zoltan2_Exceptions.hpp. */
  int _errorCheckLevel;

  /*! The output stream for error messages. */
  std::ostream *_errorOStream;
};
  
  
}  //namespace Zoltan2
  
#endif
