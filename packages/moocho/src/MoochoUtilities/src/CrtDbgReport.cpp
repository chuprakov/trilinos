// ///////////////////////////////////////////////////////////////////////////
// CrtDbgReport.cpp
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
// Just create a empty function for the linker to find.
// Runtime version only.  This is needed because of a
// difficulty with MS VC++ 6.0 when a full release
// version can not be built.

extern "C" {
	int _CrtDbgReport()
	{ return 0; }	// do nothing and return no error
}
