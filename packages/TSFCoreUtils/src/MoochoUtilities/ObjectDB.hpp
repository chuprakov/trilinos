// @HEADER
// ***********************************************************************
// 
//      TSFCoreUtils: Trilinos Solver Framework Utilities Package 
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// ///////////////////////////////////////////////
// ObjectDB.hpp

#ifndef MEM_MNG_PACK_OBJECT_DB_HPP
#define MEM_MNG_PACK_OBJECT_DB_HPP

#include <stack>
#include <deque>
#include <list>
#include <typeinfo>

#include "Teuchos_RefCountPtr.hpp"

namespace MemMngPack {

using Teuchos::RefCountPtr;

///
/** An object database based on <tt>RefCountPtr</tt>.
 *
 * The advantage of basing this database directly on
 * <tt>RefCountPtr</tt> is that any objects that are not
 * specifically removed and deleted by the client will
 * be automatically deleted by the object database itself.
 */
class ObjectDB {
public:

	///
	ObjectDB();

	///
	~ObjectDB();

	///
	/** Add an object to the object database and return its index.
	 *
	 * @param  obj_ptr   [in] Smart pointer to object to be stored.
	 *
	 * @return Returns the <tt>index</tt> that the object is given
	 * in the table (don't loose this index).
	 *
	 * Preconditions:<ul>
	 * <li> <tt>obj_ptr.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->get(index,(T*)NULL).get() == obj_ptr.get()</tt>
	 * </ul>
	 */
	template<class T>
	size_t add( const RefCountPtr<T> &obj_ptr );

	///
	/** Get the smart pointer for an object that was previously stored.
	 *
	 * @param  index     [in] Index to an object that was returned from
	 *                   <tt>this->add()</tt> when the object was
	 *                   added.
	 * @param  [noname]  [in] This is just a type selector that is used
	 *                   instantiate the correct function.  The argument
	 *                   (T*)NULL should be used.
	 *
	 * Preconditions:<ul>
	 * <li> This index and the type <tt>T</tt> must be the same as for
	 *      the initial call to <tt>this->add()</tt> or an
	 *      <tt>std::invalid_argument</tt> exception will be thrown.
	 * </ul>
	 */
	template<class T>
	RefCountPtr<T> get( size_t index, const T* ) const;
	
	///
	/** Remove an object from the object database and return its smart pointer.
	 *
	 * @param  index  [in] Index to an object that was returned from
	 *                <tt>this->add()</tt> when the object was
	 *                added.
	 *
	 * Preconditions:<ul>
	 * <li> This index must be the same as for
	 *      the initial call to <tt>this->add()</tt> or an
	 *      <tt>std::invalid_argument</tt> exception will be thrown.
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> The object is removed from this database and
	 *      <tt>this->get(index,(T*)NULL)</tt> throws an exception.
	 * </ul>
	 */
	void remove( size_t index );

public:

	typedef Teuchos::PrivateUtilityPack::RefCountPtr_node  rcp_node_t; // Must be public for object_entry_t to access!

private:

	// ////////////////////////////////
	// Private types

	struct object_entry_t {
		object_entry_t() : obj(NULL), type(NULL), rcp_node(NULL) {}
		object_entry_t(const void *_obj, const std::type_info *_type, rcp_node_t *_rcp_node)
			: obj(const_cast<void*>(_obj)), type(_type), rcp_node(_rcp_node) {}
		void                  *obj;
		const std::type_info  *type;
		rcp_node_t            *rcp_node;
	};

	typedef std::deque<object_entry_t>               objects_t;

	typedef std::stack<size_t,std::list<size_t> >    free_indexes_t;

	// ////////////////////////////////
	// Private data members

	objects_t           objects_;
	free_indexes_t      free_indexes_;

	// ///////////////////////////////
	// Private member functions

	size_t add_entry( const object_entry_t& entry );
	const object_entry_t& get_entry(size_t index, const char * typeid_name) const;
	void assert_types( const std::type_info &type_stored, const std::type_info &type_requested, size_t index ) const;

}; // class ObjectDB

// //////////////////////////////////
// Template definitions

template<class T>
inline
size_t ObjectDB::add( const RefCountPtr<T> &obj_ptr )
{
	rcp_node_t *rcp_node = obj_ptr.access_node();
	if(rcp_node) rcp_node->incr_count();
	return add_entry( object_entry_t( obj_ptr.get(), &typeid(T), rcp_node ) );
}

template<class T>
inline
RefCountPtr<T> ObjectDB::get( size_t index, const T* ) const
{
	const object_entry_t& entry = get_entry(index, typeid(T).name());
	assert_types( *entry.type, typeid(T), index );
	return RefCountPtr<T>( static_cast<T*>(entry.obj), entry.rcp_node );
}

} // namespace MemMngPack

#endif // MEM_MNG_PACK_OBJECT_DB_HPP
