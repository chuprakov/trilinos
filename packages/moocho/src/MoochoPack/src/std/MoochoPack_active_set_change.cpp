// //////////////////////////////////////////////////////////////////////////
// active_set_change.cpp

#include <assert.h>

#include <iomanip>
#include <ostream>

#include "../../include/std/active_set_change.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"

void ReducedSpaceSQPPack::active_set_change(
	  const SpVectorSlice& nu_k, const SpVectorSlice& nu_km1
	, EJournalOutputLevel olevel, size_type* num_adds, size_type* num_drops
	, std::ostream* out )
{
	using std::setw;
	using std::endl;
	
	int w = 12;

	*num_adds = *num_drops = 0;

	if( !nu_k.nz() && !nu_km1.nz() )
		return;

	assert( nu_k.is_sorted() && nu_km1.is_sorted() );

	bool dump_change = (int)olevel >= (int)PRINT_ACTIVE_SET;

	if( dump_change ) {
		*out
			<< "\n*** Changes in active set\n\n"
			<< setw(w) << "i"
			<< setw(w) << "uplo"
			<< setw(w) << "change\n"
			<< setw(w) << "----------"
			<< setw(w) << "----------"
			<< setw(w) << "----------\n";
	}

	SpVectorSlice::const_iterator
		nu_k_itr	= nu_k.begin(),
		nu_k_end	= nu_k.end(),
		nu_km1_itr	= nu_km1.begin(),
		nu_km1_end	= nu_km1.end();

	while( nu_k_itr != nu_k_end || nu_km1_itr != nu_km1_end ) {
		if( nu_k_itr != nu_k_end && ( nu_km1_itr == nu_km1_end || ( nu_k_itr != nu_k_end
				&& nu_k_itr->indice()+nu_k.offset() < nu_km1_itr->indice()+nu_km1.offset() ) ) )
		{
			// *nu_k_itr was added to active set.
			(*num_adds)++;
			if(dump_change)
				*out
					<< setw(w) << nu_k_itr->indice() + nu_k.offset()
					<< setw(w) << ( nu_k_itr->value() >= 0.0 ? "upper" : "lower" )
					<< setw(w) << "added" << endl;
			nu_k_itr++;
		}
		else if( nu_km1_itr != nu_km1_end && ( nu_k_itr == nu_k_end || ( nu_km1_itr != nu_km1_end
				&& nu_k_itr->indice()+nu_k.offset() > nu_km1_itr->indice()+nu_km1.offset() ) ) )
		{
			// *nu_km1_itr was removed from the active set.
			(*num_drops)++;
			if(dump_change)
				*out
					<< setw(w) << nu_km1_itr->indice()+nu_km1.offset()
					<< setw(w) << ( nu_km1_itr->value() >= 0.0 ? "upper" : "lower" )
					<< setw(w) << "dropped" << endl;
			nu_km1_itr++;
		}
		else {
			// same variable (but the bound may have changed)
			if( nu_k_itr->value() * nu_km1_itr->value() < 0.0 ) {
				// Switch bounds.
				(*num_adds)++;
				(*num_drops)++;
			if(dump_change)
				*out
					<< setw(w) << nu_k_itr->indice()+nu_k.offset()
					<< setw(w) << ( nu_k_itr->value() >= 0.0 ? "upper" : "lower" )
					<< setw(w) << "switch bnd" << endl;
			}
			nu_k_itr++;
			nu_km1_itr++;
		}
	}

	// Output summary
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		*out
			<< "\n*** Active set change summary\n"
			<< "nact_old  = "		<< nu_km1.nz()		<< endl
			<< "nact_new  = "		<< nu_k.nz()		<< endl
			<< "num_adds  = "		<< *num_adds		<< endl
			<< "num_drops = "		<< *num_drops		<< endl;
	}
}
