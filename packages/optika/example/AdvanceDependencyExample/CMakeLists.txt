INCLUDE(PackageAddExecutable)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE(
  AdvanceDependencyExample
  SOURCES main.cpp
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(
  AdvDepFileCopy
  SOURCE_FILES inputs.xml
  EXEDEPS AdvanceDependencyExample
)
