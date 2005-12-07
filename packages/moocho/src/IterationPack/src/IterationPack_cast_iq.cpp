// ///////////////////////////////////////////////////////////
// cast_iq.cpp
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

#include "IterationPack_cast_iq.hpp"
#include "Teuchos_TestForException.hpp"

void IterationPack::imp_cast_iq_throw_error(
	const std::string&                 iq_name
	,const std::string&                iq_is_type_name
	,const std::string&                iq_want_type_name
	)
{
	TEST_FOR_EXCEPTION(
		true, IterationPack::InvalidTypeCastException
		,"cast_id<T>(state,iq_name) : Error, the iteration quantity \'"
		<< iq_name << "\' exists with type \'" << iq_is_type_name << "\' but does not "
		<< "support the \'IterQuantityAccess<" << iq_want_type_name << ">\' interface" );
}

void IterationPack::imp_cast_iq_throw_error(
	const AlgorithmState::iq_id_type   iq_id
	,const std::string&                iq_name
	,const std::string&                iq_is_type_name
	,const std::string&                iq_want_type_name
	)
{
	TEST_FOR_EXCEPTION(
		true, IterationPack::InvalidTypeCastException
		,"cast_id<T>(state,iq_id,iq_name) : Error, the iteration quantity \'"
		<< iq_name << "\' with iq_id = \'" << iq_id
		<< "\' exists with type \'" << iq_is_type_name << "\' but does not "
		<< "support the \'IterQuantityAccess<" << iq_want_type_name << ">\' interface" );
}
