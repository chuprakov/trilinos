// /////////////////////////////////////////////////////////
// WorkspacePack.cpp
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

#include "WorkspacePack.hpp"
#include "ThrowException.hpp"

MemMngPack::ref_count_ptr<WorkspacePack::WorkspaceStore>
WorkspacePack::default_workspace_store(NULL);

namespace WorkspacePack {

// WorkspaceStore

WorkspaceStore::WorkspaceStore(size_t num_bytes)
	: workspace_begin_(NULL)
	, workspace_end_(NULL)
	, curr_ws_ptr_(NULL)
	, num_static_allocations_(0)
	, num_dyn_allocations_(0)
{
	if(num_bytes)
		protected_initialize(num_bytes);
}

WorkspaceStore::~WorkspaceStore() {
	if(workspace_begin_) delete [] workspace_begin_;
}

void WorkspaceStore::protected_initialize(size_t num_bytes)
{
	THROW_EXCEPTION(
		curr_ws_ptr_ != workspace_begin_, std::logic_error
		,"WorkspaceStore::set_workspace_size(...) : Error, "
		"You can not reset the workspace size when any RawWorkspace objects "
		"are using workspace!" );
	if(workspace_begin_) delete [] workspace_begin_;
	workspace_begin_        = ::new char[num_bytes];
	workspace_end_          = workspace_begin_ + num_bytes;
	curr_ws_ptr_            = workspace_begin_;
	num_static_allocations_ = 0;
	num_dyn_allocations_    = 0;
} 

// RawWorkspace

RawWorkspace::RawWorkspace(WorkspaceStore* workspace_store, size_t num_bytes)
{
	workspace_store_ = workspace_store;
	if( !workspace_store_ || workspace_store_->num_bytes_remaining() < num_bytes ) {
		workspace_begin_ = ::new char[num_bytes];
		workspace_end_   = workspace_begin_ + num_bytes;
		owns_memory_     = true;
		if(workspace_store_)
			workspace_store_->num_dyn_allocations_++;
	}
	else {
		workspace_begin_ = workspace_store_->curr_ws_ptr_;
		workspace_end_   = workspace_begin_ + num_bytes;
		owns_memory_     = false;
		workspace_store_->curr_ws_ptr_ += num_bytes;
		workspace_store_->num_static_allocations_++;
	}
}

RawWorkspace::~RawWorkspace()
{
	if(owns_memory_) {
		if(workspace_begin_) delete [] workspace_begin_;
	}
	else {
		THROW_EXCEPTION(
			workspace_store_->curr_ws_ptr_ != workspace_end_, std::logic_error
			,"RawWorkspace::~RawWorkspace(...): Error, "
			"Invalid usage of RawWorkspace class, corrupted WorspaceStore object!" );
		workspace_store_->curr_ws_ptr_ = workspace_begin_;
	}
}

#ifdef _PG_CXX // Should not have to define this since it should not be called!
void* RawWorkspace::operator new(size_t)
{
	assert(0);
	return NULL;
}
#endif

} // end namespace WorkspacePack
