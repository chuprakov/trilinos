// ///////////////////////////////////////////////
// ObjectDB.cpp

#include "ObjectDB.hpp"
#include "Teuchos_TestForException.hpp"

namespace MemMngPack {

ObjectDB::ObjectDB()
{}

ObjectDB::~ObjectDB()
{
	for( objects_t::iterator itr = objects_.begin(); itr != objects_.end(); ++itr ) {
                if((*itr).rcp_node && (*itr).rcp_node->deincr_count() == 0 ) delete (*itr).rcp_node;
	}
}

size_t ObjectDB::add_entry( const object_entry_t& entry )
{
	if( !free_indexes_.empty() ) {
		const size_t index = free_indexes_.top();
		free_indexes_.pop();
		objects_[index] = entry;
		return index;
	}
	objects_.push_back(entry);
	return objects_.size()-1;
}

const ObjectDB::object_entry_t& ObjectDB::get_entry(size_t index, const char * typeid_name) const
{
	TEST_FOR_EXCEPTION(
		index >= objects_.size(), std::invalid_argument
		,"ObjectDB::get_entry(index): Error, could not get object of type \'"
		<< typeid_name << "\'.  The index " << index
		<< " must be invalid since it does not fall in the range [0," << ((int)objects_.size())-1 << "]!"
		);
	const object_entry_t& entry = objects_[index];
	TEST_FOR_EXCEPTION(
		entry.obj == NULL, std::invalid_argument
		,"ObjectDB::get_entry(index): Error, the object of type \'"
		<< typeid_name << "\' with index " << index
		<< " has already been removed from the database!"
		);
	return entry;
}

void ObjectDB::remove(size_t index)
{
	TEST_FOR_EXCEPTION(
		index >= objects_.size(), std::invalid_argument
		,"ObjectDB::remove(index): Error, the index " << index
		<< " must be invalid since it does not fall in the range [0," << ((int)objects_.size()-1) << "]!"
		);
	const object_entry_t& entry = objects_[index];
	TEST_FOR_EXCEPTION(
		entry.obj == NULL, std::invalid_argument
		,"ObjectDB::remove(index): Error, the index " << index
		<< " is for an object that has been already removed from the database!"
		);
	free_indexes_.push(index);
	objects_[index] = object_entry_t();
}

void ObjectDB::assert_types( const std::type_info &type_stored, const std::type_info &type_requested, size_t index ) const
{
	TEST_FOR_EXCEPTION(
		!(type_stored==type_requested), std::invalid_argument
		,"ObjectDB::assert_types(...): Error, for object with index " << index
		<< ", the stored type \'" << type_stored.name() << "\' is not the same as the requested type \'"
		<< type_requested.name() << "\'!"
		);
}

} // namespace MemMngPack
