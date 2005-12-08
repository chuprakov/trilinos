// //////////////////////////////////////////////////////////////////////////////
// DenseLinAlgPackIOFormat.cpp
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

#include "DenseLinAlgPack_IOFormat.hpp"

// //////////////////
//  ios_format_memento

DenseLinAlgPack::LinAlgPackIO::ios_format_memento
DenseLinAlgPack::LinAlgPackIO::ios_format_memento::save_format(const std::ios& s) {

	ios_format_memento m;

	m.flags_	= s.flags();
	m.prec_		= s.precision();
	m.wdt_		= s.width();
	m.fill_		= s.fill();

	return m;
}

void DenseLinAlgPack::LinAlgPackIO::ios_format_memento::set_format(std::ios& s) const {

	s.flags(flags_);
	s.precision(prec_);
	s.width(wdt_);
	s.fill(fill_);

}

// ///////////////////
// format

void DenseLinAlgPack::LinAlgPackIO::format::copy_format(const std::ios& s) {
	ios_base_flags().flags(s.flags());
	precision(s.precision());
	width(s.width());
	fill(s.fill());
}

void DenseLinAlgPack::LinAlgPackIO::format::set_format(std::ios& s) const {
	s.flags(ios_base_flags().flags());
	s.precision(precision());
	s.width(width());
	s.fill(fill());
}
