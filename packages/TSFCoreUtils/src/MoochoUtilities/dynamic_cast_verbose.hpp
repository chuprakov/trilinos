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

#ifndef DYNAMIC_CAST_VERBOSE_H
#define DYNAMIC_CAST_VERBOSE_H

#include "TSFCoreUtils_ConfigDefs.hpp"

namespace DynamicCastHelperPack {

/** We create this class so that we may throw a bad_cast when appropriate
	and still use the TEST_FOR_EXCEPTION macro.  
	We recommend users try to catch a bad_cast.
*/
class m_bad_cast : public std::bad_cast {
	std::string msg;
public:
	explicit m_bad_cast(const std::string&  what_arg ) : msg(what_arg) {}
	virtual ~m_bad_cast() throw() {}
	virtual const char* what() const throw() { return msg.data(); }
};

// Throw the exception <tt>std::invalid_argument</tt> for below functions
void dyn_cast_throw_exception( const char type_from_name[], const char type_to_name[] );

///
/** Dynamic cast from an object of type <tt>T_From</tt> to type
 * <tt>T_To</tt> and if the cast fails at runtime throw a
 * <tt>std::invalid_argument</tt> exception with a good
 * error message.
 */
template <class T_To, class T_From>
inline
T_To& dyn_cast(T_From &from)
{
	T_To *to_ = dynamic_cast<T_To*>(&from);
	if(!to_)
		dyn_cast_throw_exception( typeid(from).name(), typeid(T_To).name() );
	return *to_;
}

}	// end namespace DynamicCastHelperPack

#endif // DYNAMIC_CAST_VERBOSE_H
