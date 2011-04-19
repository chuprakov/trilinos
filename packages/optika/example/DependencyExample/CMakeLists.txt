INCLUDE(PackageAddExecutable)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE(
  DependencyExample
  SOURCES main.cpp
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(
  DependencyExampleXMLCopy
  SOURCE_FILES inputs.xml
  EXEDEPS DependencyExample
)
