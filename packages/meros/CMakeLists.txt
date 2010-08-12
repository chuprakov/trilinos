
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Meros)

#
# B) Set up package-specific options
#

# Nothing to do here

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(thyra)

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

PACKAGE_ARCH_EXCLUDE_FILES(
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

PACKAGE_POSTPROCESS()
