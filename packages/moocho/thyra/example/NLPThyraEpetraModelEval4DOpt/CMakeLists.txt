INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)


SET(EPETRAEXT_4DOPT_DIR ${PACKAGE_SOURCE_DIR}/../epetraext/example/model_evaluator/4dopt)


INCLUDE_DIRECTORIES(${EPETRAEXT_4DOPT_DIR})


#
# NLPThyraEpetraModelEval4DOpt
#

MULTILINE_SET( XML_PARAMS_FWD_SOLVE
  "<ParameterList><Parameter name='Solve Mode' type='string' value='Forward Solve'/></ParameterList>"
  )

MULTILINE_SET( XML_PARAMS_DIRECT_LUMPED_PARAMS
  "<ParameterList>"
  "  <Parameter name='NLP Type' type='std::string' value='Direct'/>"
  "  <Parameter name='Use Parameter Lumping' type='bool' value='1'/>"
  "  <Parameter name='Show Model Evaluator Trace' type='bool' value='0'/>"
  "  <ParameterList name='Lumped Parameters Settings'>"
  "    <Parameter name='Number of Basis Columns' type='int' value='1'/>"
  "    <Parameter name='Ignore Parameter Bounds' type='bool' value='1'/>"
  "    <Parameter name='Dump Basis Matrix' type='bool' value='1'/>"
  "    <ParameterList name='VerboseObject'>"
  "      <Parameter name='Verbosity Level' type='std::string' value='high'/>"
  "    </ParameterList>"
  "  </ParameterList>"
  "  <Parameter name='Parameters Solution File Base Name' type='std::string' value='p.out'/>"
  "  <Parameter name='State Solution File Base Name' type='std::string' value='x.out'/>"
  "</ParameterList>"
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  NLPThyraEpetraModelEval4DOpt
  SOURCES
    NLPThyraEpetraModelEval4DOptMain.cpp
    ${EPETRAEXT_4DOPT_DIR}/EpetraModelEval4DOpt.cpp
  ARGS
    " "
    "--extra-moocho-thyra-params=\"${XML_PARAMS_FWD_SOLVE}\""
    "--x00=0.9"
    "--x00=0.9 --extra-moocho-thyra-params=\"${XML_PARAMS_FWD_SOLVE}\""
    "--xt0=0.9"
    "--xt0=0.9 --extra-moocho-thyra-params=\"${XML_PARAMS_FWD_SOLVE}\""
    "--x00=0.5 --xU0=0.6"
    "--x01=0.5 --xU1=0.6"
    "--x00=2.0 --xL0=1.1"
    "--x01=2.0 --xL1=1.1"
    "--x00=2.0 --xL0=1.1 --x01=2.0 --xL1=1.1"
    "--x00=0.5 --xU0=0.6 --x01=0.5 --xU1=0.6"
    "--moocho-thyra-params-file=moochoThyraFdNand.xml --no-support-derivs --p00=1.0"
    "--p00=4.0 --p01=-5.0 --extra-moocho-thyra-params=\"${XML_PARAMS_DIRECT_LUMPED_PARAMS}\""
  COMM serial mpi
  NUM_MPI_PROCS 1
  )


PACKAGE_COPY_FILES_TO_BINARY_DIR(ThyraNLPThyraEpetraModelEval4DOptCopyFiles
  DEST_FILES moochoThyraFdNand.xml
  SOURCE_PREFIX "_"
  EXEDEPS NLPThyraEpetraModelEval4DOpt
  )


ASSERT_DEFINED(PYTHONINTERP_FOUND)
IF (PYTHONINTERP_FOUND)

  ASSERT_DEFINED(PACKAGE_SOURCE_DIR)
  SET(abs_top_srcdir ${PACKAGE_SOURCE_DIR})

  ASSERT_DEFINED(PACKAGE_BINARY_DIR)
  SET(abs_top_builddir ${PACKAGE_BINARY_DIR})
  
  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/_test_lumped_param_model.py.stub.in
    ${CMAKE_CURRENT_BINARY_DIR}/test_lumped_param_model.py
    )

  PACKAGE_ADD_TEST(
    ${PYTHON_EXECUTABLE}
    NOEXEPREFIX
    NOEXESUFFIX
    NAME NLPThyraEpetraModelEval4DOpt_lumped_param_python
    ARGS test_lumped_param_model.py
    COMM serial
    STANDARD_PASS_OUTPUT
    )

ENDIF()


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_MPI)
IF (${PACKAGE_NAME}_ENABLE_MPI)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    NLPThyraEpetraMultiPointModelEval4DOptMain
    SOURCES
      NLPThyraEpetraMultiPointModelEval4DOptMain.cpp
      ${EPETRAEXT_4DOPT_DIR}/EpetraMultiPointModelEval4DOpt.cpp
    ARGS " "
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
  
ENDIF()
