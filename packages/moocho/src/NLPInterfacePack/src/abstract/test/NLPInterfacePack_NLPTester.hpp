// ///////////////////////////////////////////////////////////
// NLPTester.h

#ifndef NLP_INTERFACE_PACK_NLP_TESTER_H
#define NLP_INTERFACE_PACK_NLP_TESTER_H

#include <iosfwd>

#include "NLPInterfacePack/include/NLPInterfacePackTypes.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace NLPInterfacePack {

///
/** Testing class for base NLP interface.
 *
 * This class is little more than a unit tester for the
 * \Ref{NLP} base interface.  This class will call all
 * of the #NLP# methods and print out quanities if asked to.
 * This class simply validates the pre and post conditions
 * for all of the methods.  In that it is useful.
 */
class NLPTester {
public:

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, print_all )
	
	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, throw_exception )

	///	
	NLPTester(
		bool     print_all        = false
		,bool    throw_exception  = true
		);
 
	///
	/** Test the NLP interface as the given base point xo.
	 *
	 * ToDo: Finish Documentation!
	 */
	bool test_interface(
		NLP                           *nlp
		,const VectorWithOp           &xo
		,bool                         print_all_warnings
		,std::ostream                 *out
		);

}; // end class NLPTester

} // end namespace NLPInterfacePack

#endif // NLP_INTERFACE_PACK_NLP_TESTER_H
