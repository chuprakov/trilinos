#ifndef TSFMATRIXWRITERBASE_H
#define TSFMATRIXWRITERBASE_H

#include "TSFDefs.h"
#include "TSFLinearOperatorBase.h"
#include "TSFLinearOperator.h"
#include <string>

namespace TSF
{

  using std::string;

  /**
   * Base class for matrix writer implementations.
   *
   */
  class TSFMatrixWriterBase
    {
    public:

      /** \name Constructor and Destructor */
      //@{
      /** construct with a given filename */
      TSFMatrixWriterBase(const string& filename);

      /** the usual virtual dtor */
      virtual ~TSFMatrixWriterBase();
      //@}


      /** \name Write methods */
      //@{
      /** write and name a matrix (as in Matlab format) */
      virtual void write(const string& name, const TSFLinearOperator& A) const = 0;

      /** write a matrix */
      virtual void write(const TSFLinearOperator& A) const = 0;
      //@}

    protected:
      string filename_;

    };
};

#endif
