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

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	
// disable VC 5.0 warnings about truncated identifier names (templates).
#pragma warning(disable : 4503)	

#include <string>
#include <sstream>
 
#include "GeneralIterationPack/include/CastIQMember.h"
#include "ThrowException.h"

namespace GeneralIterationPack {

CastIQMemberBase::CastIQMemberBase( const std::string iq_name )
	:  iq_name_(iq_name), iq_id_(NOT_SET_YET)
{}

const std::string&
CastIQMemberBase::iq_name() const
{
	return iq_name_;
}

void CastIQMemberBase::cache_iq_id( const AlgorithmState& s ) const
{
	if( iq_id_ == NOT_SET_YET ) {
		AlgorithmState::iq_id_type
			tmp_iq_id = s.get_iter_quant_id( iq_name_ );
		THROW_EXCEPTION(
			tmp_iq_id == AlgorithmState::DOES_NOT_EXIST, AlgorithmState::DoesNotExist
			,"CastIQMember<T>::operator()(...) : Error, the iteration quantity \""
			<< iq_name_ << "\" does not exist in this state object." );
		iq_id_ = tmp_iq_id;
	}
}

void CastIQMemberBase::throw_cast_error( const std::string& iqa_name ) const
{
	THROW_EXCEPTION(
		true, GeneralIterationPack::InvalidTypeCastException
		,"CastIQMember<T>::operator()(state) : Error, the iteration quantity \""
		<< iq_name_ << "\" exists but it is not of the type IterQuantityAccess<"
		<< iqa_name << ">" );
}

}	// namespace GeneralIterationPack
