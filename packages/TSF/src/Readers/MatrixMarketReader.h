#ifndef MATRIXMARKETREADER_H
#define MATRIXMARKETREADER_H

#include "TSFDefs.h"
#include "TSFMatrixReaderBase.h"
#include "TSFMatrixView.h"
#include <string>
#include <iostream>

namespace TSF
{
  using std::string;
  using std::istream;

  enum MMStructure {MMGeneral, MMSymmetric, MMSkew};

  /** \ingroup MatrixReaders
   * Reads a matrix in
   * <A HREF="http://math.nist.gov/MatrixMarket"> Matrix Market </A> format.
   * The class is restricted to reading real-valued matrices only. Matrix Market sparse
   * and dense, general, symmetric, and skew formats are supported.
   *
   * */

  class MatrixMarketReader : public TSFMatrixReaderBase
    {
    public:
      /** construct with the name of the file that stores the matrix */
      MatrixMarketReader(const string& filename);

      /** virtual dtor */
      virtual ~MatrixMarketReader(){;}

      /** read a matrix */
      virtual TSFLinearOperator read(const TSFVectorType& vectorType) const ;

    private:

      void readSparse(FILE* fp, TSFMatrixView& matrix, int nRows, MMStructure structure,
                      int nnz) const ;

      void readDense(FILE* fp, TSFMatrixView& matrix, MMStructure structure) const ;

    };
}


#endif
