// //////////////////////////////////////////////////////
// Ifpack_PrecGenerator.hpp

#ifndef IFPACK_PREC_GENERATOR_HPP
#define IFPACK_PREC_GENERATOR_HPP

#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_RefCountPtr.hpp"

class Epetra_Operator;

namespace Ifpack {

///
/** Utility class that automates the creation of preconditions from Epetra_Operator objects.
 *
 * ToDo: Finish documentation!
 *
 * Note, the default copy constructor and assignment operators
 * functions are allowed since they have the correct semantics.
 * Copyies a <tt>PrecGenerator</tt> effectively just copies the
 * options (see the constructor <tt>PrecGenerator()</tt>).
 *
 * Note, this class maintains no specific state data except
 * for the options (see the constructor <tt>PrecGenerator()</tt>).
 * Therefore, the same <tt>PrecGenerator</tt> object can be used
 * to generate preconditioner objects for as many different types
 * and kinds of operators and the client would like.
 *
 * Also, once a preconditioner is created, it is a self contained
 * entity that exists indepenent of the operator it was created form
 * or from <tt>this</tt> (ToDo: Varify that this is indeed true).
 * This is the best that you could ever hope.
 */
class PrecGenerator {
public:

  ///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( int, levelFill );
  ///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( int, levelOverlap );
  ///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( double, absThreshold );
  ///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( double, relThreshold );

  ///
  /** Sets all of the adjustable options to default values.
   *
   * Note, these where the default values used in NOX as of 2004/1/16
   *
   * ToDo: Finish documentation!
   */
  PrecGenerator(
    const int        levelFill     = 1
    ,const int       levelOverlap  = 0
    ,const double    absThreshold  = 0.0
    ,const double    relThreshold  = 1.0
    );

  ///
  /** Setup (and create if not done so yet) a preconditioner.
   *
   * ToDo: Finish documentation!
   */
  void setupPrec(
    const Teuchos::RefCountPtr<Epetra_Operator>   &Op
    ,Teuchos::RefCountPtr<Epetra_Operator>        *Prec
    ) const;

}; // class PrecGenerator

} // namespace Ifpack

#endif // IFPACK_PREC_GENERATOR_HPP
