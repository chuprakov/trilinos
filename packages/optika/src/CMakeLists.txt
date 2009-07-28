INCLUDE(PackageLibraryMacros)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

FIND_PACKAGE(Qt4)

IF( NOT QT4_FOUND)
	MESSAGE(SEND_ERROR "Could not find Qt 4")
ENDIF( NOT QT4_FOUND)
IF(QT_VERSION_MINOR LESS 5)
	MESSAGE(SEND_ERROR "You must have Qt 4.5 or higher installed")
ENDIF(QT_VERSION_MINOR LESS 5)


INCLUDE( ${QT_USE_FILE} )

SET(HDRS Optika_Dependency.hpp 
Optika_DependencySheet.hpp Optika_ArrayHelperFunctions.hpp 
Optika_GUI.hpp Optika_InvalidDependencyException.hpp
Optika_SpecificParameterEntryValidators.hpp Optika_StandardDependencies.hpp
Optika_treeitem.hpp Optika_Types.hpp)


SET(MOC_HDRS Optika_metawindow.hpp Optika_treemodel.hpp
Optika_delegate.hpp Optika_ArrayWidget.hpp Optika_FilenameWidget.hpp
Optika_treeview.hpp)

SET(SRCS Optika_delegate.cpp Optika_Dependency.cpp 
Optika_DependencySheet.cpp Optika_ArrayHelperFunctions.cpp 
Optika_FilenameWidget.cpp Optika_GUI.cpp Optika_metawindow.cpp 
Optika_StandardDependencies.cpp Optika_treeitem.cpp 
Optika_treemodel.cpp Optika_treeview.cpp Optika_Types.cpp)
 
QT4_WRAP_CPP(MOC_SRCS ${MOC_HDRS} )

PACKAGE_ADD_LIBRARY(
	optika 
	HEADERS ${MOC_HDRS} ${HDRS} 
	SOURCES ${SRCS} ${MOC_SRCS} 
	DEPLIBS
	teuchos ${QT_LIBRARIES}
	)

