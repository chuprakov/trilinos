// ///////////////////////////////////////////////////////////
// CastIQMember.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	
// disable VC 5.0 warnings about truncated identifier names (templates).
#pragma warning(disable : 4503)	

#include <string>
#include <sstream>
 
#include "GeneralIterationPack/include/CastIQMember.h"

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
		if( tmp_iq_id == AlgorithmState::DOES_NOT_EXIST ) {
			std::ostringstream omsg;
			omsg
				<< "CastIQMember<T>::operator()(...) : Error, the iteration quantity \""
				<< iq_name_ << "\" does not exist in this state object.";
			throw AlgorithmState::DoesNotExist( omsg.str() );
		}
		iq_id_ = tmp_iq_id;
	}
}

void CastIQMemberBase::throw_cast_error( const std::string& iqa_name ) const
{
	std::ostringstream omsg;
	omsg
		<< "CastIQMember<T>::operator()(state) : Error, the iteration quantity \""
		<< iq_name_ << "\" exists but it is not of the type IterQuantityAccess<"
		<< iqa_name << ">";
	throw GeneralIterationPack::InvalidTypeCastException( omsg.str() );
}

}	// namespace GeneralIterationPack
