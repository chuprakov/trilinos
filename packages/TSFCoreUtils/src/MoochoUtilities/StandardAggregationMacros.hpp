// ////////////////////////////////////////////////////////////
// StandardAggregationMacros.hpp
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

#ifndef STANDARD_AGGRAGATION_MACROS_H
#define STANDARD_AGGRAGATION_MACROS_H

#include "StandardCompositionRelationshipsPack.hpp"

///
/** \defgroup StandardAggregationMacros_grp Macros that add <<std aggr>> members for an association.
 * \ingroup Misc_grp
 *
 * For example, if you want to include a <<std aggr>> association
 * with an object of type MyClass of the name my_object you
 * would include the macro in the public section of YourClass
 * declaration as follows:
 *
 \verbatim
	class YourClass {
	public:
		STANDARD_AGGREGATION_MEMBERS( MyClass, my_object )
	};
 \endverbatim
 *
 * Note that the macro addes the private member #TYPE* NAME_#
 * to the class declaration and therefore the member NAME_ is
 * available for direct access (in a constructor for example).
 *
 * In order to have a const only association use:
 \verbatim
	class YourClass {
	public:
		STANDARD_CONST_AGGREGATION_MEMBERS( MyClass, my_object )
	};
 \endverbatim
*/
//@{

/// Insert class members for a non-const association
#define STANDARD_AGGREGATION_MEMBERS( TYPE, NAME )						\
public:																	\
	void set_ ## NAME ( TYPE* NAME )									\
	{	NAME ## _ = NAME; }												\
	TYPE* get_ ## NAME()												\
	{	return NAME ## _; }												\
	const TYPE* get_ ## NAME() const									\
	{	return NAME ## _; }												\
	TYPE& NAME()														\
	{																	\
		return StandardCompositionRelationshipsPack::role_name(			\
			NAME ## _, false, " ## NAME ## " );							\
	}																	\
	const TYPE& NAME() const											\
	{																	\
		return StandardCompositionRelationshipsPack::role_name(			\
			NAME ## _, false, " ## NAME ## " );							\
	}																	\
private:																\
	TYPE* NAME ## _;													\
public:

/// Insert class members for a constant association.
#define STANDARD_CONST_AGGREGATION_MEMBERS( TYPE, NAME )				\
public:																	\
	void set_ ## NAME ( const TYPE* NAME )								\
	{	NAME ## _ = NAME; }												\
	const TYPE* get_ ## NAME() const									\
	{	return NAME ## _; }												\
	const TYPE& NAME() const											\
	{																	\
		return StandardCompositionRelationshipsPack::const_role_name(	\
			NAME ## _, false, " ## NAME ## " );							\
	}																	\
private:																\
	const TYPE* NAME ## _;												\
public:
	
#endif	// STANDARD_AGGRAGATION_MACROS_H
