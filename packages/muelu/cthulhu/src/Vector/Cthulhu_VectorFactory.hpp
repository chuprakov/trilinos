#ifndef CTHULHU_VECTOR_FACTORY_DECL_HPP
#define CTHULHU_VECTOR_FACTORY_DECL_HPP

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_Classes.hpp"

#include "Cthulhu_Map.hpp"
#include "Cthulhu_Vector.hpp"

#ifdef HAVE_CTHULHU_TPETRA
#include "Cthulhu_TpetraMap.hpp"
#include "Cthulhu_TpetraVector.hpp"
#endif
#ifdef HAVE_CTHULHU_EPETRA
#include "Cthulhu_EpetraMap.hpp"
#include "Cthulhu_EpetraVector.hpp"
#endif

#include "Cthulhu_Debug.hpp"

// This factory creates Cthulhu::Vector. User don't have to specify the exact class of object that he want to create (ie: a Cthulhu::TpetraVector or a Cthulhu::EpetraVector).
// Each Build() method takes at least one Cthulhu object in argument (a Map or another Vector) and Build() methods return a Vector created by using the same underlying library (Epetra or Tpetra).

namespace Cthulhu {
  
  template <class ScalarType, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType>

  class VectorFactory {
    
    typedef Map<LocalOrdinal, GlobalOrdinal, Node> Map;
    typedef Vector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> Vector;
#ifdef HAVE_CTHULHU_TPETRA
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMap;
    typedef TpetraVector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> TpetraVector;
#endif

  private:
    //! Private constructor. This is a static class. 
    VectorFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<Vector> Build(const Teuchos::RCP<const Map> &map, bool zeroOut=true) {
#ifdef HAVE_CTHULHU_TPETRA
      const RCP<const TpetraMap> &tMap = Teuchos::rcp_dynamic_cast<const TpetraMap>(map);
      if (tMap != null)
        return rcp( new TpetraVector(map, zeroOut) );
#endif
#ifdef HAVE_CTHULHU_EPETRA
      const RCP<const EpetraMap> &eMap = Teuchos::rcp_dynamic_cast<const EpetraMap>(map);
      if (eMap != null)
        return rcp( new EpetraVector(map, zeroOut) );
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"Cannot dynamically cast Cthulhu::Map to an EpetraMap or a TpetraMap. The exact type of the Map 'map' is unknown");
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    // Other constructors here
#endif
    
  };

}

// TODO: Only one factory for Vector and MultiVector ?? -> Yes
#endif
