// @HEADER
// ***********************************************************************
// 
//      TSFCoreUtils: Trilinos Solver Framework Utilities Package 
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

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
