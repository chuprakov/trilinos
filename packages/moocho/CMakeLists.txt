
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(MOOCHO)

#
# B) Set up package-specific options
#


#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Thyra)
IF (${PACKAGE_NAME}_ENABLE_Thyra)
  ADD_SUBDIRECTORY(thyra)
ENDIF()

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

PACKAGE_ADD_TEST_DIRECTORIES(test)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
