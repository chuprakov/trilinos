#include "MatlabWriter.h"
#include "TSFError.h"
#include "TSFMatrixView.h"

using namespace TSF;

MatlabWriter::MatlabWriter(const string& filename)
  : TSFMatrixWriterBase(filename), fout_(filename.c_str(), ios::out|ios::app)
{;}

MatlabWriter::~MatlabWriter()
{
  fout_.close();
}

void MatlabWriter::write(const string& name, const TSFLinearOperator& A) const
{
  int nRows = A.domain().dim();
  TSFMatrixView mv(A);

  for(int i=0; i<nRows; i++)
    {
      TSFArray<int> ind;
      TSFArray<TSFReal> aij;
      mv.getRow(i, ind, aij);
      
      writeRow(name, i, ind, aij);
    }
}

void MatlabWriter::write(const TSFLinearOperator& A) const
{
  TSFError::raise("Must use MatlabWriter::write(string, TSFLinearOperator)");
}

void MatlabWriter::write(const string& name, const TSFVector& v) const
{
  for(int i=0; i<v.space().dim(); i++)
    {
      if(v[i] != 0.0)
	{
	  fout_ << name << "[" << i << "] = " << v[i] << endl;
	}
    }
}

void MatlabWriter::writeRow(const string& name, int i, const TSFArray<int>& indices,
			    const TSFArray<TSFReal>& values) const
{
  int length = indices.length();

  for(int j=0; j<length; j++)
    {
      fout_ << name << "[" << i << "][" << indices[j] << "] = " << values[j] <<  endl;
    }
}
