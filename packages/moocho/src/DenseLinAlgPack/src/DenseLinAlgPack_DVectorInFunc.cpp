// /////////////////////////////////////////////////////////////////////////////////
// VectorInFunc.cpp

#include <sstream>

#include "../include/VectorInFunc.h"
#include "../include/VectorClass.h"
#include "../include/LinAlgPackAssertOp.h"

namespace {	// Local implementation
std::istream& input_vs(std::istream& is, LinAlgPack::VectorSlice* vs, const char func[]);
}

std::istream& LinAlgPack::input(std::istream& is, Vector* v, LinAlgPackIO::fmtflags extra_flags) {
	if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
		size_type n;
		is >> n;
		if(is.fail())
			throw LinAlgPackIO::InputException("LinAlgPack::input() {Vector}:  Input operation of vector dimension failed.  Check that the constant n is a valid integer.");
		if(is.bad())
			throw std::ios_base::failure("LinAlgPack::input() {Vector}: Input operation failed because the stream became currupted.");
		v->resize(n);
	}
	return input_vs(is,&(*v)(),"LinAlgPack::input() {Vector}");
}

std::istream& LinAlgPack::input(std::istream& is, VectorSlice* vs, LinAlgPackIO::fmtflags extra_flags) {
	if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
		size_type n;
		is >> n;
		if(is.fail())
			throw LinAlgPackIO::InputException("LinAlgPack::input() {VectorSlice}:  Input operation of vector dimension failed.  Check that the constant n is a valid integer.");
		if(is.bad())
			throw std::ios_base::failure("LinAlgPack::input() {VectorSlice}: Input operation failed because the stream became currupted.");
		LinAlgPack::Vp_V_assert_sizes( vs->size(), n );
	}
	return input_vs(is,vs,"LinAlgPack::input() {VectorSlice}");
}


// ///////////////////////////////////
// Local implementation

namespace {

// Read in a specified number of elements into a VectorSlice object.
// The dim of vs is not checked.  If an element input operation fails or the end of the file
// is reached before all of the elements are read in then a LinAlgPackIO::InputException is thrown.
// If the stream becomes currupted durring the input then a std::ios_base::failure exception
// is thrown.  The state of the input steam remains the same on return accept for the char's
// that have been extracted.
std::istream& input_vs(std::istream& is, LinAlgPack::VectorSlice* vs, const char func[]) {
	using std::ios_base;
	using LinAlgPack::VectorSlice;
	if(!vs->size()) return is;	// If there are no elements to read in just return
	ios_base::iostate old_state = is.exceptions();		// save the old state
	is.exceptions(ios_base::badbit | ios_base::failbit);
	try {
		// Read in the elements
		for(VectorSlice::iterator itr = vs->begin(); itr != vs->end(); ++itr)
			is >> *itr;
	}
	catch(std::ios_base::failure& excpt) {
		is.exceptions(old_state);
		if(is.bad()) throw;	// The stream was bad so rethrow the exception
		if(is.fail()) {
			std::ostringstream os;
			os << func << ":  An vector element input failed.  Check that the vector element is a valid C number.  "
			   << excpt.what();
			throw LinAlgPack::LinAlgPackIO::InputException(os.str());			
		}
		if(is.eof()) {
			std::ostringstream os;
			os << func << ":  Vector input failed.  The end of the file was found before all of the elements where read in.  "
			   << excpt.what();;
			throw LinAlgPack::LinAlgPackIO::InputException(os.str());			
		}
	}
	catch(...) {
		is.exceptions(old_state);
		throw;
	}
	is.exceptions(old_state);
	return is;
}

}	// end namespace