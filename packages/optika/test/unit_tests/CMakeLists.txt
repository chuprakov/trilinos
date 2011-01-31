if(UNIX AND (NOT APPLE) )
  message(STATUS "Target system is a non-apple, unix system. Testing ability to connect to X server")
  execute_process(
    COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/xtest.sh" 
    RESULT_VARIABLE CAN_CONNECT_TO_X
    OUTPUT_QUIET
    ERROR_QUIET)
  message(STATUS "Test Result: " ${CAN_CONNECT_TO_X})
  if(NOT ${CAN_CONNECT_TO_X} EQUAL 0)
    message(STATUS "Could not connect to an XServer. Disabling GUI unit tests.")
    SET(DO_GUI_UNIT_TESTS "FALSE")
  else(NOT ${CAN_CONNECT_TO_X} EQUAL 0)
    message(STATUS "Found X server. Enabling GUI unit tests")
    SET(DO_GUI_UNIT_TESTS "TRUE")
  endif(NOT ${CAN_CONNECT_TO_X} EQUAL 0)
else(UNIX AND (NOT APPLE) )
    message(STATUS "Not compling for a non-apple, unix system. Enableing GUI unit tests.")
    SET(DO_GUI_UNIT_TESTS "TRUE")
endif(UNIX AND (NOT APPLE) )

if(DO_GUI_UNIT_TESTS)

  INCLUDE(PackageAddExecutableAndTest)

  include(${QT_USE_FILE})
  
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
  QT4_GENERATE_MOC(GUI_UnitTests.cpp GUI_UnitTests.moc)
  
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
	  GUI_UnitTests
    SOURCES GUI_UnitTests.cpp GUI_UnitTests.moc
    FAIL_REGULAR_EXPRESSION "FAIL!"
    DEPLIBS teuchos ${TPL_QT_LIBRARIES}
    )
  
  
endif(DO_GUI_UNIT_TESTS)
