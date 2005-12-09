// ///////////////////////////////////////////////////////////////////////
// MeritFuncNLP.cpp
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

#include <typeinfo>

#include "ConstrainedOptPack_MeritFuncNLP.hpp"
#include "Teuchos_TestForException.hpp"

namespace ConstrainedOptPack {

MeritFuncNLP& MeritFuncNLP::operator=(const MeritFuncNLP& merit_func)
{
	TEST_FOR_EXCEPTION(
		this != &merit_func, std::logic_error
		,"MeritFuncNLP::operator=(merit_func) : Error, this is not assignment to self "
		"and the concreate subclass \'" << typeid(*this).name() << "\' has not overridden "
		"this method!" );
	return *this;
}

}	// end namespace ConstrainedOptPack
