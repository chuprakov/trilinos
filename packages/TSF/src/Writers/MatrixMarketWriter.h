#ifndef MATRIXMARKETWRITER_H
#define MATRIXMARKETWRITER_H

#include "TSFDefs.h"
#include "TSFMatrixWriterBase.h"
#include "mmio.h"
#include <string>
#include <iostream>

namespace TSF
{
  using std::string;
  using std::ostream;

  /** \ingroup MatrixWriters
   * Writes a matrix in Matrix Market format.
   *
   */
  class MatrixMarketWriter : public TSFMatrixWriterBase
    {
    public:
      /** construct with name of file to write */
      MatrixMarketWriter(const string& filename);

      /** destructor */
      ~MatrixMarketWriter();

      /** write a matrix */
      void write(const TSFLinearOperator& A) const;

      /** write a matrix */
      void write(const string& name, const TSFLinearOperator& A) const;

    };
};

#endif
