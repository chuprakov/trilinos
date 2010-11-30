//TODO: rename or rm this file (Cthulhu_Example)
//      + get ride of Cthulhu_Classes.hpp

#ifndef CTHULHU_MYTYPES_HPP
#define CTHULHU_MYTYPES_HPP

#include "Cthulhu_ConfigDefs.hpp"

// The purpose of this file is to allow a user to write program without worrying about templates parameters.
// For example, you can replace:
//
// typedef int LocalOrdinal;
// typedef int GlobalOrdinal;
//
// GlobalOrdinal nx = 4;
// GlobalOrdinal ny = 4;
// RCP<const Map<LocalOrdinal,GlobalOrdinal> > map = rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal>(nx*ny, 0, comm) );
//
// by:
//
// #include <Cthulhu.hpp>
//
// GO nx = 4;
// GO ny = 4;
// RCP<const Map> map = rcp( new MyMap(nx*ny, 0, comm) ); // MyMap = TpetraMap or EpetraMap according to the macro CTHULHU_USE_EPETRA and CTHULHU_USE_TPETRA.
//
// The definition of template types can still be modified if needed.

// Declaration of Cthulhu classes.
#include "Cthulhu_Classes.hpp"

// Cthulhu_DefaultTypes.hpp defines the types: ScalarType, LocalOrdinal, GlobalOrdinal and Node.
// You can use your own definition of these types by defining them before including Cthulhu.hpp. You just have to define also the macro CTHULHU_DEFAULTTYPES_HPP.
#include "Cthulhu_UseDefaultTypes.hpp"

// Cthulhu_UseShortNames.hpp get ride of template type.
// For example, instead of using the type 'Map<LocalOrdinal,GlobalOrdinal,Node>', you can simply use the short type name 'Map'.
#include "Cthulhu_UseShortNames.hpp"

// You can optionally define which implementation of the linar algebra layer you want to use (ie: Epetra or Tpetra) by defining the macro CTHULHU_USE_EPETRA or CTHULHU_USE_TPETRA.
// Instead of using directly TpetraMap or EpetraMap when allocating a Cthulhu::Map, you can then use the type 'MyMap'. This type is define according to CTHULHU_USE_EPETRA and CTHULHU_USE_TPETRA.
// This allows you to switch from one implementation to another at compile time. Another way is to use a Cthulhu::MapFactory (selection at runtime).

#if defined (CTHULHU_USE_TPETRA) && !defined(CTHULHU_USE_EPETRA)

#ifndef HAVE_CTHULHU_TPETRA
#error You must enable Tpetra
#endif

#include <Cthulhu_TpetraMap.hpp>
#include <Cthulhu_TpetraCrsMatrix.hpp>
#include <Cthulhu_TpetraVector.hpp>
#include <Cthulhu_TpetraMultiVector.hpp>

//TODO: remove
typedef Cthulhu::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> MyMap;
typedef Cthulhu::TpetraCrsMatrix<ScalarType, LocalOrdinal, GlobalOrdinal, Node> MyCrsMatrix;
//typedef Cthulhu::TpetraVector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> MyVector; //TODO!
typedef Cthulhu::TpetraMultiVector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> MyMultiVector;


typedef Cthulhu::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> DefaultMap;
typedef Cthulhu::TpetraCrsMatrix<ScalarType, LocalOrdinal, GlobalOrdinal, Node> DefaultCrsMatrix;
//typedef Cthulhu::TpetraVector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> DefaultVector; //TODO!
typedef Cthulhu::TpetraMultiVector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> DefaultMultiVector;

#elif defined(CTHULHU_USE_EPETRA) && !defined(CTHULHU_USE_TPETRA)

#ifndef HAVE_CTHULHU_EPETRA
#error You must enable Epetra
#endif

#include <Cthulhu_EpetraMap.hpp>
#include <Cthulhu_EpetraCrsMatrix.hpp>
//#include <Cthulhu_EpetraVector.hpp>
#include <Cthulhu_EpetraMultiVector.hpp>

typedef Cthulhu::EpetraMap MyMap;
typedef Cthulhu::EpetraCrsMatrix MyCrsMatrix;
//typedef Cthulhu::EpetraVector MyVector;
typedef Cthulhu::EpetraMultiVector MyMultiVector;

typedef Cthulhu::EpetraMap DefaultMap;
typedef Cthulhu::EpetraCrsMatrix DefaultCrsMatrix;
//typedef Cthulhu::EpetraVector DefaultVector;
typedef Cthulhu::EpetraMultiVector DefaultMultiVector;

#endif // CTHULHU_USE_E/TPETRA

#endif // CTHULHU_MYTYPES_HPP
