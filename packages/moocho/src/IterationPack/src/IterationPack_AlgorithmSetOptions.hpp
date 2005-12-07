// ////////////////////////////////////////////////////////////////
// IterationPack_AlgorithmSetOptions.hpp
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

#ifndef ITERATION_PACK_ALGORITHM_SET_OPTIONS_H
#define ITERATION_PACK_ALGORITHM_SET_OPTIONS_H

#include "IterationPack_Algorithm.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace IterationPack {

///
/** Set options for <tt>Algorithm</tt> from an <tt>OptionsFromStream</tt> object.
 */
class AlgorithmSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			Algorithm >
{
public:

	///
	AlgorithmSetOptions(
		  Algorithm* target = 0
		, const char opt_grp_name[] = "IterationPack::Algorithm" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void setOption( int option_num, const std::string& option_value );

};	// end class AlgorithmSetOptions

}	// end namespace IterationPack

#endif	// ITERATION_PACK_ALGORITHM_SET_OPTIONS_H
