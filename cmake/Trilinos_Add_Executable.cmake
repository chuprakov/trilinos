#2008/07/09: rabartl: ToDo::
#
# (*) Change the name of the current TRILINOS_ADD_TARGET(...) to 
#     TRILINOS_ADD_EXECUTABLE(...)
#
# (*) Add an optional DIRECTORY argument and put the executable in that directory
#
# (*) Add same COMM logic for deciding whether to add an executable depending
#     on serial or COMM
#
# (*) ??? 
#

INCLUDE(Parse_Variable_Arguments)

SET(CMAKE_EXECUTABLE_SUFFIX ".exe")

FUNCTION (TRILINOS_ADD_EXECUTABLE EXE_NAME)
   PARSE_ARGUMENTS(
     PARSE  #prefix
     "SOURCES;DIRECTORY;COMM" #lists
     "" #options
     ${ARGN} )
  

  SET(ADD_THE_EXE OFF)
  IF(TRILINOS_ENABLE_MPI)
   
    IF (NOT PARSE_COMM)
      # If no COMM is given assume we will add the exe
     
      SET(ADD_THE_EXE ON)
    ELSE()
      # Else, if COMM is defined we have to find 'mpi'
      LIST (FIND PARSE_COMM "mpi" DO_MPI_INDEX)
    ENDIF()
    IF(NOT ${DO_MPI_INDEX} EQUAL -1)
      
      SET(ADD_THE_EXE ON)    
    ENDIF()
  ELSE()
    
    IF (NOT PARSE_COMM)
      # If no COMM is given assume we will add the exe
      
      SET(ADD_THE_EXE ON)
    ELSE()
      # Else, if COMM is defined we have to find 'serial'
      LIST (FIND PARSE_COMM "serial" DO_SERIAL_INDEX)
    ENDIF()
    IF(NOT ${DO_SERIAL_INDEX} EQUAL -1)
      
      SET(ADD_THE_EXE ON)    
    ENDIF()
  ENDIF()
  
  IF(ADD_THE_EXE)  
    SET (EXE_SOURCES)
    SET(EXE_BINARY_NAME ${EXE_NAME})
    IF(PARSE_DIRECTORY ) #If exe is in subdirectory prepend that dir name to the source files
      FOREACH( SOURCE_FILE ${PARSE_SOURCES} )
        SET (EXE_SOURCES ${EXE_SOURCES} ${PARSE_DIRECTORY}/${SOURCE_FILE})
      ENDFOREACH( )
    ELSE()
      FOREACH( SOURCE_FILE ${PARSE_SOURCES} )
        SET (EXE_SOURCES ${EXE_SOURCES} ${SOURCE_FILE})
      ENDFOREACH( )
    ENDIF()
    IF("${VERBOSE_CONFIGURE}" STREQUAL "ON")
      MESSAGE("TRILINOS_ADD_EXECUTABLE: ${EXE_BINARY_NAME}")
    ENDIF()
    ADD_EXECUTABLE(${EXE_BINARY_NAME} ${EXE_SOURCES})
    IF(PARSE_DIRECTORY)
      SET_TARGET_PROPERTIES( ${EXE_BINARY_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PARSE_DIRECTORY} )
    ENDIF()
  ENDIF()
  
  
ENDFUNCTION()

# Setup include directories and library dependencies
INCLUDE_DIRECTORIES(${${PROJECT_NAME}_INCLUDE_DIRS})
LINK_LIBRARIES(${${PROJECT_NAME}_LIBRARIES})
