
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  verySimple
  crsgraph
  )

IF (${PACKAGE_NAME}_ENABLE_Amesos)

  ADD_SUBDIRECTORIES(
    amesos
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_AztecOO)

  ADD_SUBDIRECTORIES(
    aztecoo
    )

ENDIF()
