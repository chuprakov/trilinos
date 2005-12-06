// /////////////////////////////////////////////////////////////////////////
// OptionsFromStreamPack_SetOptionsFromStreamNode.hpp
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

#ifndef SET_OPTIONS_FROM_STREAM_NODE_H
#define SET_OPTIONS_FROM_STREAM_NODE_H

#include "OptionsFromStreamPack_SetOptionsFromStream.hpp"
#include "OptionsFromStreamPack_StringToIntMap.hpp"

namespace OptionsFromStreamPack {

///
/** Node class for setting options from a stream.
  *
  * This class uses the template method pattern to
  * delegate the setting of options.
  */
class SetOptionsFromStreamNode: public SetOptionsFromStream {
public:

	///
	/** Constructs with the name of the options group and the names
	  * of the options.
	  *
	  *	@param	options_group	The name of the options group to access
	  *	@param	num_options		The number of options in the opitons
	  *							group.
	  *	@param	option_name		An array (length num_options) containing
	  *							the names of the options.
	  *	@param	exists_optional	Specifies if the options group must exist.
	  */
	SetOptionsFromStreamNode( const std::string& options_group
		, int num_options, const char* option_names[]
		, bool exists_optional = true );

	///
	/** Overridden from SetOptionsFromStream and calls setOption(...).
	  *
	  * The options group #options_group# is used.  If this options
	  * group does not exist and #exists_optional# == false then
	  * an #std::invalid_argument# exception will be thrown.
	  */
	void set_options( const OptionsFromStream& options );

protected:

	///
	/** To be overridden by the subclass to set an option given
	  * its integer position and the option value.
	  *
	  * The integer possition returned is the possition of the option
	  * in option_names[option_num] that was passed to the constructor.
	  */
	virtual void setOption( int option_num, const std::string& option_value ) = 0;

private:
	StringToIntMap	name_map_;
	bool			exists_optional_;

};	// end class SetOptionsFromStreamNode

}	// end namespace OptionsFromStreamPack

#endif	// SET_OPTIONS_FROM_STREAM_NODE_H
