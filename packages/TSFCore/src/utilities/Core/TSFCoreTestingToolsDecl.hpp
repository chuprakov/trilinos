// /////////////////////////////////////////////////////////////////////////
// TSFCoreTestingToolsDecl.hpp

#ifndef TSFCORE_TESTING_TOOLS_DECL_HPP
#define TSFCORE_TESTING_TOOLS_DECL_HPP

#include <iosfwd>

#include "TSFCoreTypes.hpp"

namespace TSFCore {

///
inline
const char* toString(ETransp transp) { return transp == NOTRANS ? "NOTRANS" : "TRANS"; }

///
inline
const char* passfail(bool pass) { return pass ? "passed" : "failed"; }

///
template<class Scalar>
std::ostream& operator<<( std::ostream& o, const Vector<Scalar>& v );

///
template<class Scalar>
std::ostream& operator<<( std::ostream& o, const LinearOp<Scalar>& M );

} // namespace TSFCore

#endif // TSFCORE_TESTING_TOOLS_DECL_HPP
