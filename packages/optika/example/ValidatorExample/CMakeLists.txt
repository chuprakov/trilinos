INCLUDE(PackageAddExecutable)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE(
  ValidatorExample
  SOURCES main.cpp
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(
  ValidatorExampleFileCopy
  SOURCE_FILES inputs.xml
  EXEDEPS ValidatorExample
)

