// /////////////////////////////////////////////////////////////////////////
// OptionsFromStreamPack_SetOptionsToTargetBase.hpp
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

#ifndef SET_OPTIONS_TO_TARGET_BASE_H
#define SET_OPTIONS_TO_TARGET_BASE_H

#include "StandardCompositionRelationshipsPack.hpp"

namespace OptionsFromStreamPack {

///
/** Templated node class manipulating a reference to a target object
  * who will have its options set..
  */
template< class T >
class SetOptionsToTargetBase {
public:

	/// 
	SetOptionsToTargetBase( T* target = 0 )
		: target_(target)
	{}

	/** @name <<std aggr>> stereotype members for target.
	  */
	//@{

	///
	void set_target(T* target)
	{	target_ = target; }
	///
	T* get_target()
	{	return target_; }
	///
	const T* get_target() const
	{	return target_; }
	///
	T& target()
	{	return StandardCompositionRelationshipsPack::role_name(target_, false, "target"); }
	///
	const T& target() const
	{	return StandardCompositionRelationshipsPack::role_name(target_, false, "target"); }

	//@}

private:
	T* target_;

};	// end class SetOptionsToTargetBase

}	// end namespace OptionsFromStreamPack

#endif	// SET_OPTIONS_TO_TARGET_BASE_H
