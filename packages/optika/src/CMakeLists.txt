INCLUDE(PackageLibraryMacros)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

include(${QT_USE_FILE})

SET(HDRS Optika_Dependency.hpp 
Optika_ArrayHelperFunctions.hpp 
Optika_GUI.hpp 
Optika_SpinBoxApplier.hpp
Optika_treeitem.hpp Optika_Types.hpp )


SET(MOC_HDRS Optika_metawindow.hpp Optika_treemodel.hpp
Optika_delegate.hpp Optika_ArrayWidget.hpp Optika_FileNameWidget.hpp
Optika_treeview.hpp)

SET(SRCS Optika_delegate.cpp 
Optika_ArrayHelperFunctions.cpp 
Optika_FileNameWidget.cpp Optika_GUI.cpp Optika_metawindow.cpp 
Optika_treeitem.cpp 
Optika_treemodel.cpp Optika_treeview.cpp Optika_Types.cpp
)
 
QT4_WRAP_CPP(MOC_SRCS ${MOC_HDRS} )


PACKAGE_ADD_LIBRARY(
	optika 
	HEADERS ${MOC_HDRS} ${HDRS} 
	SOURCES ${SRCS} ${MOC_SRCS} 
	DEPLIBS
	teuchos ${TPL_Qt_LIBRARIES}
	)

