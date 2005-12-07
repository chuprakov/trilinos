// ///////////////////////////////////////////////////////////
// CastIQMember.cpp
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
#include <sstream>
 
#include "IterationPack_CastIQMember.hpp"
#include "Teuchos_TestForException.hpp"

namespace IterationPack {

CastIQMemberBase::CastIQMemberBase( const std::string iq_name )
	:  iq_name_(iq_name), iq_id_(NOT_SET_YET)
{}

const std::string&
CastIQMemberBase::iq_name() const
{
	return iq_name_;
}

bool
CastIQMemberBase::exists_in( const AlgorithmState& s ) const
{
	cache_iq_id(s);
	return iq_id_ != AlgorithmState::DOES_NOT_EXIST;
}

void CastIQMemberBase::cache_iq_id( const AlgorithmState& s ) const
{
	if( iq_id_ == NOT_SET_YET ) {
		iq_id_ = s.get_iter_quant_id( iq_name_ );
	}
}

void CastIQMemberBase::throw_cast_error( const AlgorithmState::iq_id_type iq_id, const std::string& iqa_name ) const
{
	TEST_FOR_EXCEPTION(
		iq_id == AlgorithmState::DOES_NOT_EXIST, AlgorithmState::DoesNotExist
		,"CastIQMember<T>::operator()(...) : Error, the iteration quantity \""
		<< iq_name_ << "\" does not exist in this state object." );
	TEST_FOR_EXCEPTION(
		true, IterationPack::InvalidTypeCastException
		,"CastIQMember<T>::operator()(state) : Error, the iteration quantity \""
		<< iq_name_ << "\" exists but it is not of the type IterQuantityAccess<"
		<< iqa_name << ">" );
}

}	// namespace IterationPack
