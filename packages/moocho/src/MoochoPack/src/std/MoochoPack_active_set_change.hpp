// //////////////////////////////////////////////////////////////////////////
// active_set_change.h

#ifndef ACTIVE_SET_CHANGE_H
#define ACTIVE_SET_CHANGE_H

#include <iosfwd>

#include "../ReducedSpaceSQPPackTypes.h"

namespace ReducedSpaceSQPPack {

///
/** Calculate the change in the active set and output change
  * if asked to.
  *
  * ToDo: Add more description of the output you get.
  *
  * @param	nu_k	[I]	Multipliers for variable bounds for iteration k.
  * @param	num_km1	[I] Multipliers for variable bounds for iteration k-1
  * @param	olevel	[I] Specifies the output level.  We have:\\
  *				PRINT_NOTHING : No output is sent to out\\
  *				PRINT_ALGORITHM_STEPS :
  *					Just the number of additions
  *					and deletions to the active set and the total number
  *					of active constraints is output.\\ 
  *				PRINT_ACTIVE_SET : Enumerates
  *					which variable were added and dropped from the active set.\\
  * @param	num_adds	[O] Gives the number of variables fixed at a bound
  *					added to the active set.
  * @param	num_drops	[O] Gives the number of variables freed from a 
  *					bound and dropped from the active set.
  * @param	out	[O] Target for output.
  */
void active_set_change( const SpVectorSlice& nu_k, const SpVectorSlice& nu_km1
	, EJournalOutputLevel olevel, size_type* num_adds, size_type* num_drops
	, std::ostream* out ); 

}	// end namespace ReducedSpaceSQPPack

#endif	// ACTIVE_SET_CHANGE_H