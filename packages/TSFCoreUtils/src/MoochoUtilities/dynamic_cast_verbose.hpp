// //////////////////////////////////////////////////////////////////////////
// dynamic_cast_verbose.hpp
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
