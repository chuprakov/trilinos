// //////////////////////////////////////////////////////////////////////////////
// LinAlgPackIO_NameLookups.hpp
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
//

#ifndef LINALGPACK_IO_NAME_LOOKUPS_H
#define LINALGPACK_IO_NAME_LOOKUPS_H

#ifdef _WINDOWS

//  MS VC++ 5.0 is not performing function look for these templated operator
//  functions as it should so I will do it for the compiler.  These
//	inline functions are injected into the local namespace.  These should
//  be removed once a standards conforming compiler is available.

namespace {

// Define some boiler plate macros

#define OPEATOR_FUNCTION(OPERATOR,STREAM_TYPE,FORMAT_TYPE,OBJECT_TYPE)								\
	inline STREAM_TYPE & OPERATOR ( STREAM_TYPE & s													\
		, LinAlgPack::LinAlgPackIO:: ## FORMAT_TYPE ## <LinAlgPack:: ## OBJECT_TYPE ## >& bf)	\
	{																								\
		return LinAlgPack:: ## OPERATOR ## (s,bf);													\
	}

#define INPUT_OPEATOR_FUNCTION(FORMAT_TYPE,OBJECT_TYPE)												\
	OPEATOR_FUNCTION( operator>> , std::istream , FORMAT_TYPE , OBJECT_TYPE )						\

#define OUTPUT_OPEATOR_FUNCTION(FORMAT_TYPE,OBJECT_TYPE)											\
	OPEATOR_FUNCTION( operator<< , std::ostream , FORMAT_TYPE , OBJECT_TYPE )						\


INPUT_OPEATOR_FUNCTION(		bound_format		,	Vector			)
OUTPUT_OPEATOR_FUNCTION(	bound_format		,	Vector			)
OUTPUT_OPEATOR_FUNCTION(	const_bound_format	,	Vector			)

INPUT_OPEATOR_FUNCTION(		bound_format		,	VectorSlice		)
OUTPUT_OPEATOR_FUNCTION(	bound_format		,	VectorSlice		)
OUTPUT_OPEATOR_FUNCTION(	const_bound_format	,	VectorSlice		)

INPUT_OPEATOR_FUNCTION(		bound_format		,	GenMatrix		)
OUTPUT_OPEATOR_FUNCTION(	bound_format		,	GenMatrix		)
OUTPUT_OPEATOR_FUNCTION(	const_bound_format	,	GenMatrix		)

INPUT_OPEATOR_FUNCTION(		bound_format		,	GenMatrixSlice	)
OUTPUT_OPEATOR_FUNCTION(	bound_format		,	GenMatrixSlice	)
OUTPUT_OPEATOR_FUNCTION(	const_bound_format	,	GenMatrixSlice	)

#undef OPEATOR_FUNCTION
#undef INPUT_OPEATOR_FUNCTION
#undef OUTPUT_OPEATOR_FUNCTION

}	// end namespace

#endif

#endif // LINALGPACK_IO_NAME_LOOKUPS_H
