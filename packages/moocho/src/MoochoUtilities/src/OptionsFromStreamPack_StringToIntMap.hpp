// /////////////////////////////////////////////////////
// OptionsFromStreamPack_StringToIntMap.hpp
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

#ifndef STRING_TO_INT_MAP_H
#define STRING_TO_INT_MAP_H

#include <map>
#include <string>
#include <stdexcept>
#include <sstream>

#include "OptionsFromStreamPack_OptionsFromStreamExceptions.hpp"
#include "Teuchos_TestForException.hpp"

namespace OptionsFromStreamPack {

///
/** Map a string to an enumeration.
  *
  * The purpose of this class is to simplify mapping a standard string
  * to an integer which can be interpreted as an enumeration.
  *
  * Here is an example of its use.
  \verbatim

	const int n_opt = 3;
	enum MyOptEnum {
		OPT_ONE
		,OPT_TWO
		,OPT_THREE
	};	// must be 0, 1,..., n_opt - 1
	const char* MyOptStrings[n_opt] = {
		"OPT_ONE
		,"OPT_TWO"
		,"OPT_THREE"
	}; // parallels MyOptEnum
	StringToIntMap my_enum_map( "opt_map", n_opt, NyOptStrings );
	...
	switch( (MyEnum)my_enum_map( "OPT_ONE" ) ) {
		case OPT_ONE:
			// do stuff
		case OPT_TWO:
			// do stuff
		case OPT_THREE:
			// do stuff
		default:
			// ???
	}

  \endverbatim
  * 
  * The number of strings passed to the constructor must equal the number
  * of options in the enumeration.  If there are duplicate strings
  * (capitalization concidered) then the exception #AlreadyExists# is
  * throw.  If a string that was not passed in the
  * constructor if given to #operator()( const std::string& str )# then
  * the exception #DoesNotExist# is thrown.
  *
  * In the constructor, #name# is used in error messages in the exceptions
  * thrown to help make since out of the message.
  *
  * The default constructor is not defined and not to be called.
  */
class StringToIntMap {
public:

	///
	class AlreadyExists : public std::logic_error
	{public: AlreadyExists(const std::string& what_arg) : std::logic_error(what_arg) {}};

	///
	class DoesNotExist : public AccessException
	{public: DoesNotExist(const std::string& what_arg) : AccessException(what_arg) {}};

	///
	StringToIntMap( const std::string& name, int n, const char* strings[] );

	///
	int operator()( const std::string& str ) const;

	///
	const std::string& name() const;

private:
	typedef std::map< std::string, int > map_t;	// all share implementation.
	std::string name_;
	map_t map_;

	// not defined and not to be called.
	StringToIntMap();

};	// end class StringToIntMap

// ////////////////////////////////////////////
// Inline declarations

inline
const std::string& StringToIntMap::name() const
{
	return name_;
}

}	// end namespace OptionsFromStreamPack 

#endif	// STRING_TO_INT_MAP_H
