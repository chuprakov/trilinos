// //////////////////////////////////////////////////////
// Ifpack_PrecGenerator.cpp

#include "Ifpack_PrecGenerator.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Ifpack_CrsRiluk.h"
#include "dynamic_cast_verbose.hpp"

// RAB: 2004/01/19: It turns out that as of this writting the Ifpack
// preconditioner class Ifpack_CrsRiluk does not allow reuse with
// InitValues(...) and Factor(...)  which is unfortunate!
#define IFPACK_PRECGENERATOR_CREATE_PREC_EVERY_TIME  1

namespace Ifpack {

PrecGenerator::PrecGenerator(
  const int        levelFill
  ,const int       levelOverlap
  ,const double    absThreshold
  ,const double    relThreshold
  )
  :levelFill_(levelFill)
  ,levelOverlap_(levelOverlap)
  ,absThreshold_(absThreshold)
  ,relThreshold_(relThreshold)
{}

void PrecGenerator::setupPrec(
  const Teuchos::RefCountPtr<Epetra_Operator>   &Op
  ,Teuchos::RefCountPtr<Epetra_Operator>        *Prec_in
  ) const
{
  using DynamicCastHelperPack::dyn_cast;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::set_extra_data;
  using Teuchos::get_extra_data;

  RefCountPtr<Epetra_Operator> &Prec = *Prec_in;

  // Determine what type of Epetra_Operator we have for DcDy
  const Epetra_RowMatrix *Op_rm = dynamic_cast<const Epetra_RowMatrix*>(&*Op);
  if(Op_rm) {
    // Yea, We can use AztecOO or Ifpack preconditioners!  Now
    // determine if Op_rm is actually a a Crs or a Vbr matrix.
    const Epetra_CrsMatrix *Op_crs
      = dynamic_cast<const Epetra_CrsMatrix*>(Op_rm);
    const Epetra_VbrMatrix *Op_vbr
      = ( Op_crs ? (const Epetra_VbrMatrix*)NULL : dynamic_cast<const Epetra_VbrMatrix*>(Op_rm) );
    if( Op_crs || Op_vbr ) { 
      //
      // Op is a Crs or a Vbr matrix!
      //
      // Create the preconditioner if it has not been created already.
			const bool createPrec
#ifdef IFPACK_PRECGENERATOR_CREATE_PREC_EVERY_TIME
				= true
#else
				= !Prec.get()
#endif
				;
			if( createPrec ) {
        RefCountPtr<Ifpack_IlukGraph>  Op_iluk_graph;
#ifdef IFPACK_PRECGENERATOR_CREATE_PREC_EVERY_TIME
				if(Prec.get()) {
					// Get the precreated graph
					Op_iluk_graph = get_extra_data<RefCountPtr<Ifpack_IlukGraph> >(Prec,0);
				}
				else {
#endif
					// Create the graph
					Op_iluk_graph = rcp(
            new Ifpack_IlukGraph(
              ( Op_crs
                ? Op_crs->Graph()
                : Op_vbr->Graph() )
              ,levelFill()
              ,levelOverlap()
              )
            );
					Op_iluk_graph->ConstructFilledGraph();
#ifdef IFPACK_PRECGENERATOR_CREATE_PREC_EVERY_TIME
				}
#endif
				// Create the preconditioner
        Prec = rcp(new Ifpack_CrsRiluk(*Op_iluk_graph));
        // Give ownership of the graph to the preconditioner.  Note,
        // if multiple operator and multiple preconditioner objects
        // that have the same structure are to be used, this
        // implementation will create a new graph for each of these
        // preconditioner objects.  However, as long as these
        // preconditioner objects are reused over and over again then
        // this approach should not be a performance problem and this
        // greatly simplifies how this class is used.
        set_extra_data( Op_iluk_graph, &Prec );
			}
      // Get a Ifpack_CrsRiluk subclass pointer for prec
      Ifpack_CrsRiluk *Prec_crs_riluk = &dyn_cast<Ifpack_CrsRiluk>(*Prec);
      // Now initialize the values
      if(Op_crs)
        Prec_crs_riluk->InitValues(*Op_crs);
      else if(Op_vbr)
        Prec_crs_riluk->InitValues(*Op_vbr);
      else
        assert(0); // Should never get here!
      // Set diagonal perturbations
      Prec_crs_riluk->SetAbsoluteThreshold(absThreshold());
      Prec_crs_riluk->SetRelativeThreshold(relThreshold());
      // Finally, complete the factorization
      Prec_crs_riluk->Factor();
    }
    else {
      // It turns out that Op only supports the Epetra_RowMatrix interface!
      TEST_FOR_EXCEPTION(
        true,std::logic_error
        ,"PrecGenerator::setupPrec(...): Have not implemented support for only Epetra_RowMatrix interface yet!");
    }
  }
  else {
    TEST_FOR_EXCEPTION(
      true,std::logic_error
      ,"PrecGenerator::setupPrec(...):Error, can't handle operator that does not support Epetra_RowMatrix!");
  }
}

} // namespace Ifpack
