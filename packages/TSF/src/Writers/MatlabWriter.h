#ifndef MATLABWRITER_H
#define MATLABWRITER_H

#include "TSFDefs.h"
#include "TSFMatrixWriterBase.h"
#include <string>
#include <iostream>
#include <fstream>

namespace TSF
{
  using std::string;
  using std::ostream;

  /** \ingroup MatrixWriters
   * Writes a matrix in Matlab format.  The particular format of the
   * file is as follows: multiple lines of the form
   *
   *  i   j   value
   *
   * where i is the row index and j the column index of a nonzero value.
   *
   *  Notes:
   *
   *     1. The first line contains
   *             numRows  numCols  0.0
   *        so that the right amount of space will be allocated by
   *        Matlab when it loads the file.  This allows for a zero
   *        block at the end of the matrix.
   *
   *     2. You should only write one object (matrix or vector) to any
   *     given file.  If more is written Matlab will not be able to
   *     load it properly.
   *
   *     3. To load this file into Matlab, do the following: (Assume
   *     that MatlabWriter has been constructed with the string
   *     "myMat.dat".)
   *
   *           temp = load('myMat.dat');
   *           A = sparse(temp(:,1), temp(:,2), temp(:,3));
   *     This creates a sparse matrix A with the i, j, values data as above.
   *     For more information, RTFM for Matlab.
   *
   * */
  class MatlabWriter : public TSFMatrixWriterBase
    {
    public:
      /** construct with name of file to write */
      MatlabWriter(const string& filename);

      /** destructor closes file stream */
      ~MatlabWriter();

      /** write a matrix to the file */
      void write(const TSFLinearOperator& A) const;

      /** \name Write methods */
      //@{
      /** write and name a matrix (as in Matlab format) */
      virtual void write(const string& name, const TSFLinearOperator& A) const
        {
          write(A);
        }

      /** write a vector to the file */
      void write(const TSFVector& v) const;

    private:

      void writeRow(int i, const TSFArray<int>& indices,
                    const TSFArray<TSFReal>& values) const;

      mutable ofstream fout_;

    };
};

#endif
