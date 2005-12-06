// //////////////////////////////////////////////////////////////////////
// TestingHelperPack_update_success.cpp
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
#include <string>

#include "TestingHelperPack_update_success.hpp"
#include "Teuchos_TestForException.hpp"

bool TestingHelperPack::update_success(bool result_check, bool* success) {
	if(result_check == false) {
		TEST_FOR_EXCEPTION(
			throw_except_on_fail, std::runtime_error
			,"update_success(...) : Runtime check "
			"failed and throw_except_on_fail == false."	);	
		*success = false;
	}
	return result_check;
}

bool TestingHelperPack::throw_except_on_fail = false;
