// ///////////////////////////////////////////////////////////
// cast_iq.h
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

#ifndef CAST_IQ_H
#define CAST_IQ_H

#include <sstream>
#include <stdexcept>
#include <typeinfo>

#include "AlgorithmState.h"
#include "IterQuantityAccess.h"

namespace GeneralIterationPack {

///
/** Lookup an iteration quantity by name and cast it
  * to an IterQuantityAccess<T> of the given type T.
  * If the iteration quantity of that name does not
  * exist then a AlgorithmState::DoesNotExist exception
  * will be thrown.  If the type of the iteration quantity
  * is not of the type IterQuantityAcess<T> (as determined
  * by dynamic_cast<T>) then the exception InvalidTypeCastException:
  * will be thrown with a helpful error message.
  */
template<class T>
IterQuantityAccess<T>& cast_iq( AlgorithmState& state
	, const std::string& iq_name );

template<class T>
const IterQuantityAccess<T>& cast_iq( const AlgorithmState& state
	, const std::string& iq_name );

}	// namespace GeneralIterationPack

// ///////////////////
// Inline definitions

namespace {

template<class T>
inline
void imp_cast_iq_throw_error(  const std::string& iq_name )
{
	std::ostringstream omsg;
	omsg
		<< "cast_id<T>(state,iq_name) : Error, the iteration quantity \""
		<< iq_name << "\" exists but it is not of the type IterQuantityAccess<"
		<< typeid(T).name() << ">";
	throw GeneralIterationPack::InvalidTypeCastException( omsg.str() );
}

}	// end namespace

namespace GeneralIterationPack {

template<class T>
inline
IterQuantityAccess<T>& cast_iq( AlgorithmState& state
	, const std::string& iq_name )
{
	IterQuantityAccess<T>
		*p = dynamic_cast<IterQuantityAccess<T>*>( &state.iter_quant( iq_name ) );
			// will throw exception if iq_name does not exist
	if( !p )
		imp_cast_iq_throw_error<T>( iq_name );
	return *p;	
}

template<class T>
inline
const IterQuantityAccess<T>& cast_iq( const AlgorithmState& state
	, const std::string& iq_name )
{
	const IterQuantityAccess<T>
		*p = dynamic_cast<const IterQuantityAccess<T>*>( &state.iter_quant( iq_name ) );
			// will throw exception if iq_name does not exist
	if( !p )
		imp_cast_iq_throw_error<T>( iq_name );
	return *p;	
}

}	// namespace GeneralIterationPack

#endif	// CAST_IQ_H
