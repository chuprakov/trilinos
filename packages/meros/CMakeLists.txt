
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE(Meros)

#
# B) Set up package-specific options
#

# Nothing to do here

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(thyra)

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_EXCLUDE_FILES(
  README-MEROS
  doc/UsersGuide
  example/data/mac
  example/data/mac-vbr
  example/data/q1p0
  example/data/salsa
  example/data/tmac-vbr
 )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
