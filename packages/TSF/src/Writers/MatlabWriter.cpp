#include "MatlabWriter.h"
#include "TSFError.h"
#include "TSFMatrixView.h"

using namespace TSF;

MatlabWriter::MatlabWriter(const string& filename)
  : TSFMatrixWriterBase(filename), fout_(filename.c_str())
{;}

MatlabWriter::~MatlabWriter()
{
  fout_.close();
}

void MatlabWriter::write(const TSFLinearOperator& A) const
{
  int nRows = A.range().dim();
  int nCols = A.domain().dim();
  //fout_ << name << " = sparse(" << nRows << " , " << nCols << ");" << endl;
  fout_ << nRows << "  " << nCols << "  0.0" << endl;

  //TSFMatrixView mv(A);

  for(int i=0; i<nRows; i++)
    {
      TSFArray<int> ind;
      TSFArray<TSFReal> aij;
      A.getRow(i, ind, aij);
      //cerr << "writing row " << i << endl;
      writeRow(i, ind, aij);
    }
}


void MatlabWriter::write(const TSFVector& v) const
{
  fout_ << v.space().dim() << "  1  0.0" << endl;
  for(int i=0; i<v.space().dim(); i++)
    {
      if(v[i] != 0.0)
        {
          fout_ << i+1 << "  1  " << v[i] << endl;
        }
    }
}

void MatlabWriter::writeRow(int i, const TSFArray<int>& indices,
			    const TSFArray<TSFReal>& values) const
{
  int length = indices.length();

  for(int j=0; j<length; j++)
    {
      //      fout_ << name << "(" << i+1 << " , " << indices[j]+1 << ") = " << values[j] <<  ";" << endl;
      fout_ << i+1 << "  " << indices[j]+1  << "  " << values[j] << endl;
    }
}
