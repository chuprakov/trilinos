// /////////////////////////////////////////////
// RTOpCppToC.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#ifndef RTOP_CPP_TO_C_H
#define RTOP_CPP_TO_C_H

#include "RTOpCpp.hpp"

namespace RTOpPack {

///
/** Small utility class that will get (or create) a <tt>RTOp_RTOp</tt> interface
 * from an <tt>RTOpPack::RTOp</tt> operator object.
 *
 * ToDo: Finish documentation!
 */
class RTOpCppToC {
public:
	///
	RTOpCppToC( const RTOp &op_cpp );
	///
	~RTOpCppToC();
	///
	const RTOpPack::RTOp& op_cpp() const;
	///
	const RTOp_RTOp& op_c() const;
private:
	const RTOp        &op_cpp_;
	RTOp_RTOp         *op_c_;
	RTOp_RTOp_vtbl_t  *op_vtbl_;
	// Not defined and not to be called!
	RTOpCppToC();
	RTOpCppToC(const RTOpCppToC&);
	RTOpCppToC& operator=(const RTOpCppToC&);
}; // end class RTOpCppToC

// //////////////////////////
// Inline members

inline
const RTOpPack::RTOp& RTOpCppToC::op_cpp() const
{
	return op_cpp_;
}

inline
const RTOp_RTOp& RTOpCppToC::op_c() const
{
	return *op_c_;
}

} // end namspace RTOpPack

#endif // RTOP_CPP_TO_C_H
