// //////////////////////////////////////////////////////////////////////
// QPKWIK_Output.h
//
// These are subroutines that are called by QPKWIK to output data for
// debugging purposes.

#include <ostream>

namespace QPKWIK_Output {

/// 
/** Class for setting the output stream (destructor unsets it).
  *
  * Create an object of this type to set a stream for outputing
  * before you call QPKWIK.
  * This is not tread safe.
  */
class set_output {
public:
	///
	set_output(std::ostream* out);
	///
	~set_output();
private:
	// not defined and not to be called
	set_output();
	set_output(const set_output&);
	set_output& operator=(const set_output&);
};	// end class set_output

// Output stream to use (default == 0, no output).
extern std::ostream* out;

}	// end namespace QPKWIK_Output
