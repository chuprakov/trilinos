// ////////////////////////////////////////////////////////////
// StandardCompositionMacros.hpp
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

#ifndef STANDARD_COMPOSITION_MACROS_H
#define STANDARD_COMPOSITION_MACROS_H

#include "ref_count_ptr.hpp"

///
/** Macro that adds <<std comp>> members for a composition association.
  * \ingroup StandardContainmentMacros_grp
  *
  * This form is for when the object being held will have const attributes
  * the same as the <tt>this</tt> object.
  *
  * For example, if you want to include a <<std comp>> association
  * with an non-const object of type MyClass of the name my_object you
  * would include the macro in the public section of YourClass
  * declaration as follows:
  *
  \verbatim
	class YourClass {
	public:
		STANDARD_COMPOSITION_MEMBERS( MyClass, my_object )
	};
  \endverbatim
  *
  * Note that the macro adds the following type and data members
  * to the class declaration:<br>
  \verbatim
	public:
		MemMngPack::ref_count_ptr<TYPE> NAME_ptr_t;
	private:
		NAME_ptr_t NAME_;		
  \endverbatim
  *
  * Note that since this macro defines a public type it must be
  * included before any other declarations that may need to refere to
  * the type NAME_ptr_t (such as constuctors)
  */
#define STANDARD_COMPOSITION_MEMBERS( TYPE, NAME )					\
	typedef MemMngPack::ref_count_ptr<TYPE>				\
		NAME ## _ptr_t;												\
	void set_ ## NAME (const NAME ## _ptr_t& NAME )					\
	{	NAME ## _ = NAME ; }										\
	const NAME ## _ptr_t& get_ ## NAME() const						\
	{	return NAME ## _; }											\
	TYPE& NAME()													\
	{	return *NAME ## _; }										\
	const TYPE& NAME() const										\
	{	return *NAME ## _; }										\
private:															\
	NAME ## _ptr_t	NAME ## _;										\
public:


///
/** Macro that adds <<std comp>> members for a composition association.
  * \ingroup StandardContainmentMacros_grp
  *
  * This form is for when the object being held will have non-const attributes
  * irrespective of the const of <tt>this</tt>.
  *
  * For example, if you want to include a <<std comp>> association
  * with an non-const object of type MyClass of the name my_object you
  * would include the macro in the public section of YourClass
  * declaration as follows:
  *
  \verbatim
	class YourClass {
	public:
		STANDARD_NONCONST_COMPOSITION_MEMBERS( MyClass, my_object )
	};
  \endverbatim
  *
  * Note that the macro adds the following type and data members
  * to the class declaration:<br>
  \verbatim
	public:
		MemMngPack::ref_count_ptr<TYPE> NAME_ptr_t;
	private:
		NAME_ptr_t NAME_;		
  \endverbatim
  *
  * Note that since this macro defines a public type it must be
  * included before any other declarations that may need to refere to
  * the type NAME_ptr_t (such as constuctors)
  */
#define STANDARD_NONCONST_COMPOSITION_MEMBERS( TYPE, NAME ) \
	typedef MemMngPack::ref_count_ptr<TYPE> \
		NAME ## _ptr_t; \
	void set_ ## NAME (const NAME ## _ptr_t& NAME ) \
	{	NAME ## _ = NAME ; } \
	const NAME ## _ptr_t& get_ ## NAME() const \
	{	return NAME ## _; } \
	TYPE& NAME() const \
	{	return *NAME ## _; } \
private: \
	NAME ## _ptr_t	NAME ## _; \
public:

///
/** Macro that adds <<std comp>> members for a composition association.
  * \ingroup StandardContainmentMacros_grp
  *
  * This form is for when the object being held will have const attributes
  * irrespective of the const of <tt>this</tt>.
  *
  * For example, if you want to include a <<std comp>> association
  * with a const object of type MyClass of the name my_object you
  * would include the macro in the public section of YourClass
  * declaration as follows:
  *
  \verbatim
	class YourClass {
	public:
		STANDARD_CONST_COMPOSITION_MEMBERS( MyClass, my_object )
	};
  \endverbatim
  *
  * Note that the macro adds the following type and data members
  * to the class declaration:<br>
  \verbatim
	public:
		MemMngPack::ref_count_ptr<const TYPE> NAME_ptr_t;
	private:
		NAME_ptr_t NAME_;		
  \endverbatim
  *
  * Note that since this macro defines a public type it must be
  * included before any other declarations that may need to refere to
  * the type NAME_ptr_t (such as constuctors)
  */
#define STANDARD_CONST_COMPOSITION_MEMBERS( TYPE, NAME )			\
public:																\
	typedef MemMngPack::ref_count_ptr<const TYPE>		\
		NAME ## _ptr_t;												\
	void set_ ## NAME (const NAME ## _ptr_t& NAME )					\
	{	NAME ## _ = NAME ; }										\
	NAME ## _ptr_t& get_ ## NAME()									\
	{	return NAME ## _ ; }										\
	const NAME ## _ptr_t& get_ ## NAME() const						\
	{	return NAME ## _; }											\
	const TYPE& NAME() const										\
	{	return *NAME ## _; }										\
private:															\
	NAME ## _ptr_t	NAME ## _;										\
public:

#endif	// STANDARD_COMPOSITION_MACROS_H
