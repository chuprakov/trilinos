// //////////////////////////////////////////////////////////
// ReleaseResource.hpp
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

#ifndef RELEASE_RESOURCE_H
#define RELEASE_RESOURCE_H

namespace MemMngPack {

///
/** Abstract interface for releasing an object when it is
 * not needed anymore {abstract}.
 *
 * The purpose of this object is so that a client can give another
 * peer a means to release needed resources when that object is
 * not bound anymore.
 */
class ReleaseResource {
public:

	///
	/** When object is deleted so will the resource if it is not
	 * needed anymore.
	 */
	virtual ~ReleaseResource() {}

	///
	/** Returns true if a resource is bound to this object.
	 */
	virtual bool resource_is_bound() const = 0;

};

}  // end namespace MemMngPack

#endif // RELEASE_RESOURCE_H
