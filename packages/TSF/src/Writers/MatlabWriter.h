#ifndef MATLABWRITER_H
#define MATLABWRITER_H

#include "TSFConfig.h"
#include "TSFMatrixWriterBase.h"
#include <string>
#include <iostream>
#include <fstream>

namespace TSF
{
  using std::string;
  using std::ostream;

  /** \ingroup MatrixWriters
   * Writes a matrix in Matlab format.
   *
   */
  class MatlabWriter : public TSFMatrixWriterBase
    {
    public:
      /** construct with name of file to write */
      MatlabWriter(const string& filename);

      /** destructor closes file stream */
      ~MatlabWriter();

      /** write a matrix to a given name */
      void write(const string& name, const TSFLinearOperator& A) const;

      /** write a matrix */
      void write(const TSFLinearOperator& A) const;

      /** write a vector to a given name */
      void write(const string& name, const TSFVector& v) const;

    private:

      void writeRow(const string& name, int i, const TSFArray<int>& indices, 
		    const TSFArray<TSFReal>& values) const;

      mutable ofstream fout_;

    };
};

#endif
