// /////////////////////////////////////////////////////////////
// OptionsFromStreamPack_StringToIntMap.cpp
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

#include "OptionsFromStreamPack_StringToIntMap.hpp"

namespace OptionsFromStreamPack {

StringToIntMap::StringToIntMap(  const std::string& name, int n, const char* strings[] )
	: name_(name)
{
	typedef map_t::value_type val_t;
	for( int i = 0; i < n; ++i ) {
		const bool unique = map_.insert( val_t( strings[i], i ) ).second;
		TEST_FOR_EXCEPTION(
			!unique,AlreadyExists
			,"StringToIntMap::StringToIntMap(...): "
			<< "Error, the option \"" << strings[i] << "\" is a duplicate for options_group \""
			<< name_ << "\"" );
	}
}

int StringToIntMap::operator()( const std::string& str ) const
{
	map_t::const_iterator itr = map_.find( str );
	TEST_FOR_EXCEPTION(
		itr == map_.end(), DoesNotExist
		,"StringToIntMap::operator(...): "
		<< "Error, the option \"" << str << "\" is not recongnised for options_group \""
		<< name_ << "\"" );
	return (*itr).second;	
}

} // end namespace OptionsFromStreamPack
