// ///////////////////////////////////////////////////////
// print_vector_change_stats.cpp

#include "../include/print_vector_change_stats.h"
#include "../include/vector_change_stats.h"

void ConstrainedOptimizationPack::print_vector_change_stats(
	  const VectorSlice& x, const char x_name[]
	, const VectorSlice& d, const char d_name[], std::ostream& out )
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
