// ///////////////////////////////////////////////////////////
// assert_print_nan_inf.cpp

#include <ostream>
#include <sstream>
#include <iomanip>

#include "LinAlgPack/include/assert_print_nan_inf.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "Misc/include/check_nan_inf.h"

bool LinAlgPack::assert_print_nan_inf( const value_type& val, char name[]
	, bool throw_excpt, std::ostream* out )
{
	using NumericHelperPack::is_nan;
	using NumericHelperPack::is_inf;
	
	if( is_nan(val) || is_inf(val) ) {
		std::ostringstream omsg;
		omsg
			<< "The scalar \"" << name
			<< "\" = " << val << " is not a valid bounded number";
		if(out)
			*out << omsg.str() << std::endl;
		if( throw_excpt ) {
			if(out)
				out->flush();	
			throw NaNInfException( "assert_print_nan_inf(...) : Error, "
				+ omsg.str() );
		}
		return false;
	}
	return true;
}

bool LinAlgPack::assert_print_nan_inf( const VectorSlice& v, char name[]
	, bool throw_excpt, std::ostream* out )
{
	using NumericHelperPack::is_nan;
	using NumericHelperPack::is_inf;
	
	bool has_nan_or_inf = false;
	bool printed_header = false;

	for( VectorSlice::const_iterator v_itr = v.begin(); v_itr != v.end(); ++v_itr ) {
		if( is_nan(*v_itr) || is_inf(*v_itr) ) {
			if(out) {
				if(!printed_header) {
					*out
						<< "The vector \"" << name
						<< "\" has the following NaN or Inf entries\n";
					printed_header = true;
				}
				*out
					<< name << "(" << v_itr - v.begin() + 1 << ") = "
					<< *v_itr << std::endl;
			}
			has_nan_or_inf = true;
		}
	}
	if( has_nan_or_inf && throw_excpt ) {
		if(out)
			out->flush();	
		std::ostringstream omsg;
		omsg
			<< "assert_print_nan_inf(...) : Error, the vector named "
			<< name << " has at least one element which is NaN or Inf";
		throw NaNInfException( omsg.str() );
	}

	return !has_nan_or_inf;
}

bool LinAlgPack::assert_print_nan_inf( const GenMatrixSlice& m, char name[]
	, bool throw_excpt, std::ostream* out )
{
	using NumericHelperPack::is_nan;
	using NumericHelperPack::is_inf;
	
	bool has_nan_or_inf = false;
	bool printed_header = false;

	for( size_type j = 1; j <= m.cols(); ++j ) {
		const VectorSlice& v = m.col(j);
		for( VectorSlice::const_iterator v_itr = v.begin(); v_itr != v.end(); ++v_itr ) {
			if( is_nan(*v_itr) || is_inf(*v_itr) ) {
				if(out) {
					if(!printed_header) {
						*out
							<< "The matrix \"" << name
							<< "\" has the following NaN or Inf entries\n";
						printed_header = true;
					}
					*out
						<< name << "(" << v_itr - v.begin() + 1 << "," << j << ") = "
						<< *v_itr << std::endl;
				}
				has_nan_or_inf = true;
			}
		}
	}
	
	if( has_nan_or_inf && throw_excpt ) {
		if(out)
			out->flush();	
		std::ostringstream omsg;
		omsg
			<< "assert_print_nan_inf(...) : Error, the matrix named "
			<< name << " has at least one element which is NaN or Inf";
		throw NaNInfException( omsg.str() );
	}

	return has_nan_or_inf;
}
