#ifndef PETRAVECTORTYPE_H
#define PETRAVECTORTYPE_H

#include "TSFDefs.h"
#include "TSFVectorTypeBase.h"


#define PETRA_BOOL_SUPPORTED
#include "TSFVectorSpaceBase.h"
#include "Epetra_Map.h"

namespace TSF
{
  class TSFMatrixOperator;
  /** \ingroup VectorSpace
   *
   */

  class PetraVectorType : public TSFVectorTypeBase
    {
    public:
      /** empty ctor */
      PetraVectorType(){;}
      /** TUVD */
      virtual ~PetraVectorType(){;}

      /** create a vector space */
      virtual TSFVectorSpace createSpace(int dimension) const ;

      /** create a vector space */
      virtual TSFVectorSpace createSpace(int dimension, int nLocal,
                                         int firstLocal) const ;

      /** create a vector space */
      virtual TSFVectorSpace createSpace(int dimension, int nLocal,
                                         const int* localIndices) const ;

      /** create a vector space */
      virtual TSFVectorSpace createSpace(int dimension, int nLocal,
                                         const int* localIndices,
                                         int nGhost,
                                         const int* ghostIndices) const ;

      /** */
      virtual TSFMatrixOperator* createMatrix(const TSFVectorSpace& domain,
                                              const TSFVectorSpace& range) const ;


      /** write to a stream */
      ostream& print(ostream& os) const {return os << "PetraVectorType";}

      /** */
      virtual TSFLinearSolver defaultSolver() const ;

    private:
      static Epetra_Comm& comm();
    };
}


#endif
