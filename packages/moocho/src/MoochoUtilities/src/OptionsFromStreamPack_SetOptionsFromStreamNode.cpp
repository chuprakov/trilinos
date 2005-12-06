// /////////////////////////////////////////////////////////////////////////
// OptionsFromStreamPack_SetOptionsFromStreamNode.cpp
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

#include <sstream>

#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_OptionsFromStream.hpp"
#include "Teuchos_TestForException.hpp"

namespace OptionsFromStreamPack {

SetOptionsFromStreamNode::SetOptionsFromStreamNode( const std::string& options_group
		, int num_options, const char* option_names[], bool exists_optional )
	: name_map_( options_group, num_options, option_names )
		, exists_optional_(exists_optional)
{}

void SetOptionsFromStreamNode::set_options( const OptionsFromStream& options )
{
	OptionsFromStream::options_group_t optgrp = options.options_group( name_map_.name() );
	if( !OptionsFromStream::options_group_exists( optgrp ) ) {
		TEST_FOR_EXCEPTION(
			!exists_optional_, std::invalid_argument
			,"SetOptionsFromStreamNode::set_options(...) : "
			<< "Error, The options group " << name_map_.name()
			<< " does not exist" );
		if(exists_optional_)
			return;
	}
		
	OptionsFromStream::options_group_t::const_iterator itr = optgrp.begin();
	for( ; itr != optgrp.end(); ++itr ) {
		setOption( name_map_( option_name(itr) ), option_value(itr).c_str() );
	}
}

}	// end namespace OptionsFromStreamPack
