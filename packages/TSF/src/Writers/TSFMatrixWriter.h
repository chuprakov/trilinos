#ifndef TSFMATRIXWRITER_H
#define TSFMATRIXWRITER_H

#include "TSFConfig.h"
#include "TSFMatrixWriterBase.h"
#include <string>

namespace TSF
{

  using std::string;

  /**
   * TSFMatrixWriter: User-level handle class for writing matrices.
   *
   */
  class TSFMatrixWriter
    {
    public:
      /** \name Constructor and Destructor */
      //@{
      /** construct with a pointer to a derived writer type. */
      TSFMatrixWriter(TSFMatrixWriterBase* ptr)
	: ptr_(ptr) {;}

      //@}
      /** \name Write methods */
      //@{
      /** Write and name a matrix (as in Matlab format) */
      void write(const string& name, const TSFLinearOperator A) const
	{ptr_->write(name, A);}

      /** Write a matrix */
      void write(const TSFLinearOperator& A) const
	{ptr_->write(A);}
      //@}
      

    private:
      TSFSmartPtr<TSFMatrixWriterBase> ptr_;
    };
}
      

#endif
