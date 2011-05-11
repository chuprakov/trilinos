INCLUDE(PackageAddExecutable)
INCLUDE(PackageCopyFilesToBinaryDir)
PACKAGE_ADD_EXECUTABLE(
  PLConfigurator 
  SOURCES main.cpp
  NOEXEPREFIX
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(
  PL_ICON_COPY
  SOURCE_FILES "plicon.svg"
  EXEDEPS PLConfigurator
  NOEXEPREFIX
)
  

