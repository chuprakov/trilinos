// ///////////////////////////////////////////////////////////////////////
// VectorSpaceTester.h

#ifndef VECTOR_SPACE_TESTER_H
#define VECTOR_SPACE_TESTER_H

#include <iosfwd>

#include "AbstractLinAlgPack/include/AbstractLinAlgPackTypes.h"
#include "StandardMemberCompositionMacros.h"

namespace AbstractLinAlgPack {

///
/** Testing class for VectorSpace, VectorWithOp and VectorWithOpMutable.
 *
 * ToDo: Finish documentation.
 */
class VectorSpaceTester {
public:

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, print_all_tests )
	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, print_vectors )
	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, throw_exception )
	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_random_tests )
	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, warning_tol )
	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, error_tol )

	///	
	VectorSpaceTester(
		bool         print_all_tests  = false
		,bool        print_vectors    = false
		,bool        throw_exception  = true
		,size_type   num_random_tests = 4
		,value_type  warning_tol      = 1e-14
		,value_type  error_tol        = 1e-10
		);

	///
	virtual ~VectorSpaceTester() {}

	///
	/** Run a vector space through a set of comprehensive tets.
	 *
	 * ToDo: Finish documentation.
	 */
	virtual bool check_vector_space(
		const VectorSpace &space
		,std::ostream     *out
		) const;

private:

	///
	void check_test(value_type err, std::ostream* out, bool* success) const;

}; // end class VectorSpaceTester

} // end namespace AbstractLinAlgPack

#endif // VECTOR_SPACE_TESTER_H
