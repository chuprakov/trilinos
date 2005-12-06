// //////////////////////////////////////////////////
// FortranTypes_f_open_file.hpp
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

#ifndef F_OPEN_FILE_H
#define F_OPEN_FILE_H

#include <stdexcept>

#include "Teuchos_F77_wrappers.h"

namespace FortranTypes {

/** @name Open a Fortran file.
  */
//@{

///
enum EOpenStatus { OPEN_OLD = 0, OPEN_NEW = 1, OPEN_SCRATCH = 2
	, OPEN_UNKNOWN = 3 };
///
enum EOpenForm { OPEN_FORMATTED = 0, OPEN_UNFORMATTED = 1 };
///
enum EOpenBlank { OPEN_NULL = 0, OPEN_ZERO = 1 };
///
enum EOpenAccess { OPEN_SEQUENTIAL = 0, OPEN_DIRECT = 1 };

/** Open a Fortran file given its name and unit number.
  *
  * If successful #iunit# is returned for the opened file.
  *
  * The standard options to Fortran OPEN(...) are included.  The
  * only mandatory option to set is for the file name in FILE.
  *
  * The length of the file name must be <= 100 characters long.
  *
  * Note that all of the options are optional but if
  * you set access == OPEN_DIRECT you must set recl to some
  * value greater than zero.  See Metcalf, 1990, p. 127.
  *
  * If #file# is not a valid ASCII string then the exception
  * InvalidFileNameException will be thrown
  *
  * If the file could not be opened for some reason then
  * the exception OpenException will be thrown.
  */
void f_open_file( const f_int iunit, const char file[]
	, EOpenStatus status = OPEN_UNKNOWN, EOpenForm form = OPEN_FORMATTED
	, EOpenBlank blank = OPEN_NULL, EOpenAccess access = OPEN_SEQUENTIAL
	, f_int recl = -1 );

/** Close a Fortran file given its unit number.
  *
  * If successful #iunit# is returned for the opened file.
  *
  * The standard options to Fortran OPEN(...) are included.  The
  * only mandatory option to set is for the file name in FILE.
  *
  * The length of the file name must be <= 100 characters long.
  *
  * Note that all of the options are optional but if
  * you set access == OPEN_DIRECT you must set recl to some
  * value greater than zero.  See Metcalf, 1990, p. 127.
  *
  * If #file# is not a valid ASCII string then the exception
  * InvalidFileNameException will be thrown
  *
  * If the file could not be opened for some reason then
  * the exception OpenException will be thrown.
  */
void f_close_file( const f_int iunit, bool keep = true );

//@}

/// Thrown if the file name is not a valid ASCII string
class InvalidFileNameException : public std::logic_error
{public: InvalidFileNameException(const std::string& what_arg) : std::logic_error(what_arg) {}};

/// Thrown if the open operation fails
class OpenException : public std::logic_error
{public: OpenException(const std::string& what_arg) : std::logic_error(what_arg) {}};

}	// end namespace FortranTypes 

#endif // F_OPEN_FILE_H
