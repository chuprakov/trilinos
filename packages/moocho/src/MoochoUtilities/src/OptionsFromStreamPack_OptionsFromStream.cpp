// /////////////////////////////////////////////////////////////////////////
// OptionsFromStreamPack_OptionsFromStream.cpp
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

#include <string>
#include <ostream>
#include <istream>

#include "OptionsFromStreamPack_OptionsFromStream.hpp"
#include "InputStreamHelperPack_EatInputComment.hpp"
#include "Teuchos_TestForException.hpp"

namespace {

// Remove white space from beginning and end.
inline
void clip_ws( std::string* str ) {
	// Remove the first non ' ' characters
	{
		std::string::iterator itr;
		for( itr = str->begin(); itr != str->end(); ++itr )
			if( *itr != ' ' ) break;
		str->erase( str->begin(), itr ); 
	}
	// Remove the last non ' ' characters
	{
		std::string::iterator itr;
		for( itr = str->end() - 1; itr > str->begin() - 1; --itr )
			if( *itr != ' ' ) break;
		str->erase( itr + 1, str->end() ); 
	}
}

}	// end namespace

namespace OptionsFromStreamPack {

namespace OptionsFromStreamUtilityPack {

// class OptionsGroup

std::string OptionsGroup::option_does_not_exist_;

}	// end namespace OptionsFromStreamUtilityPack

// class OptionsFromStream

OptionsFromStream::options_group_t
OptionsFromStream::options_group( const std::string& options_group_name )
{
	options_group_map_t::iterator itr = options_group_map_.find( options_group_name );
	if( itr != options_group_map_.end() ) {
		(*itr).second.second.set(true); // flag that we have accessed this options group
		return options_group_t(&(*itr).second.first);
	}
	return options_group_t(NULL);
}

const OptionsFromStream::options_group_t
OptionsFromStream::options_group( const std::string& options_group_name ) const
{
	options_group_map_t::const_iterator itr = options_group_map_.find( options_group_name );
	if( itr != options_group_map_.end() ) {
		const_cast<false_bool_t&>((*itr).second.second).set(true); // flag that we have accessed this options group
		return options_group_t(const_cast<option_to_value_map_t*>(&(*itr).second.first));
	}
	return options_group_t(NULL);
}

void OptionsFromStream::reset_unaccessed_options_groups()
{
	for( iterator og_itr = begin() ; og_itr != end(); ++og_itr ) {
		(*og_itr).second.second.set(false);	// rest to not accessed yet.
	}
}

void OptionsFromStream::print_unaccessed_options_groups( std::ostream& out ) const
{
	const_iterator og_itr = begin();
	for( ; og_itr != end(); ++og_itr ) {
		if( (*og_itr).second.second == false ) // if options group was not accessed
			out << "options_group " << options_group_name( og_itr ) << " {}\n";
	}
}

void OptionsFromStream::read_options( std::istream& in )
{
	using std::getline;
	
	using InputStreamHelperPack::eat_comment_lines;

	std::string curr_word;
	int ignore_len = 100;

	if(!in)
		return; // No options to read!

	// Eat words until you get to begin_options.
	while( curr_word != "begin_options" && !in.eof() )
		in >> curr_word;
	// Loop through each options group.
	while( !in.eof() ) {
	 	eat_comment_lines(in,'*');
		// read options_group
		in >> curr_word;
		if( curr_word == "end_options" )
			return;
		TEST_FOR_EXCEPTION(
			curr_word != "options_group", InputStreamError
			,"OptionsFromStream::read_options(...) : "
			"Error, curr_word = \'" << curr_word << " != \'options_group\'" );
		// read the name of the options group up to {
		in >> curr_word;
		if( curr_word[curr_word.length()-1] == '{' )
			curr_word.erase( curr_word.length()-1, 1 );	// remove '{' from the end.
		else
			in.ignore(ignore_len,'{');	// remove (hopefully whitespace) up to '{'
		clip_ws( &curr_word );
		// Access the options and values map for this options group.
		option_to_value_map_t& optval = options_group_map_[curr_word].first;
		// Loop through and add the options.
	 	eat_comment_lines(in,'*');
		in >> curr_word; // should be an option name.
		while( curr_word != "}" ) {
			std::string &value = optval[curr_word];
			in.ignore( ignore_len, '=' );
			// Read in the value.
			getline( in, value, ';' );
			clip_ws( &value );
			// Read the next option name.
		 	eat_comment_lines(in,'*');
			in >> curr_word;
		}
	}
}

void OptionsFromStream::print_options( std::ostream& out ) const {
	out << "\nbegin_options\n";
	const_iterator og_itr = begin();
	for( ; og_itr != end(); ++og_itr ) {
		out << "\noptions_group " << options_group_name( og_itr ) << " {\n";
		const options_group_t optgrp = OptionsFromStreamPack::options_group( og_itr );
		options_group_t::const_iterator itr = optgrp.begin();
		for( ; itr != optgrp.end(); ++itr ) {
			const std::string
				&name  = option_name(itr),
				&value = option_value(itr);
			out	<< "    " << name << " = " << value << ";\n";
		}
		out << "}\n";
	}
	out << "\nend_options\n\n";
}

}	// end namespace OptionsFromStreamPack
