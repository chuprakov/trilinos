// //////////////////////////////////////////////
// dynamic_cast_verbose.cpp
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

#include <stdexcept>
#include <sstream>

#include "dynamic_cast_verbose.hpp"
#include "ThrowException.hpp"

void DynamicCastHelperPack::dyn_cast_throw_exception(
	const char T_from[], const char T_to[] )
{
	THROW_EXCEPTION(
		true, std::invalid_argument
		,"dyn_cast<" << T_to << ">(" << T_from
		<< ") : Error, the object with the concrete type \'"
		<< T_from << "\' does not support the interface \'"
		<< T_to << "\' and the dynamic cast failed!" );
}
