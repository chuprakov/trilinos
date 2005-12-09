// ///////////////////////////////////////////////////////
// print_vector_change_stats.cpp
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

#include "ConstrainedOptPack_print_vector_change_stats.hpp"
#include "ConstrainedOptPack_vector_change_stats.hpp"

void ConstrainedOptPack::print_vector_change_stats(
	  const DVectorSlice& x, const char x_name[]
	, const DVectorSlice& d, const char d_name[], std::ostream& out )
{
	value_type	max_term,	min_term,	av_term;
	size_type	max_k,		min_k;
	vector_change_stats(
		  x, d
		, &max_term, &max_k
		, &min_term, &min_k
		, &av_term	);
	out	<< "\nmax(|"<<d_name<<"(i)|/(1+|"<<x_name<<"(i)|)"
			<< " => |"<<d_name<<"("<<max_k<<")|/(1+|"<<x_name<<"("<<max_k<<")| = "<< max_term
		<< "\nmin(|"<<d_name<<"(i)|/(1+|"<<x_name<<"(i)|)"
			<< " => |"<<d_name<<"("<<min_k<<")|/(1+|"<<x_name<<"("<<min_k<<")| = "<< min_term
		<< "\naverage(|"<<d_name<<"(i)|/(1+|"<<x_name<<"(i)|) = " << av_term << std::endl;
}
