// //////////////////////////////////////////
// TSFCoreEpetraTypes.hpp

#ifndef TSFCORE_EPETRA_TYPES_HPP
#define TSFCORE_EPETRA_TYPES_HPP

#include "TSFCoreTypes.hpp"

class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Operator;

namespace TSFCore {

class EpetraVectorSpace;
class EpetraVectorSpaceFactory;
class EpetraVector;
class EpetraVectorMultiVector;
class EpetraLinearOp;

} // namespace TSFCore

#endif // TSFCORE_EPETRA_TYPES_HPP
