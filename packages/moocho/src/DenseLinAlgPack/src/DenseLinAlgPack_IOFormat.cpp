// //////////////////////////////////////////////////////////////////////////////
// LinAlgPackIOFormat.cpp

#include "../include/LinAlgPackIOFormat.h"

// //////////////////
//  ios_format_memento

LinAlgPack::LinAlgPackIO::ios_format_memento
LinAlgPack::LinAlgPackIO::ios_format_memento::save_format(const std::ios& s) {

	ios_format_memento m;

	m.flags_	= s.flags();
	m.prec_		= s.precision();
	m.wdt_		= s.width();
	m.fill_		= s.fill();

	return m;
}

void LinAlgPack::LinAlgPackIO::ios_format_memento::set_format(std::ios& s) const {

	s.flags(flags_);
	s.precision(prec_);
	s.width(wdt_);
	s.fill(fill_);

}

// ///////////////////
// format

void LinAlgPack::LinAlgPackIO::format::copy_format(const std::ios& s) {
	ios_base_flags().flags(s.flags());
	precision(s.precision());
	width(s.width());
	fill(s.fill());
}

void LinAlgPack::LinAlgPackIO::format::set_format(std::ios& s) const {
	s.flags(ios_base_flags().flags());
	s.precision(precision());
	s.width(width());
	s.fill(fill());
}
