#include "MatrixMarketWriter.h"
#include "TSFError.h"
#include "TSFMatrixView.h"

using namespace TSF;

MatrixMarketWriter::MatrixMarketWriter(const string& filename)
  : TSFMatrixWriterBase(filename)
{;}

MatrixMarketWriter::~MatrixMarketWriter()
{
}

void MatrixMarketWriter::write(const TSFLinearOperator& A) const
{
  
  MM_typecode t;
  mm_initialize_typecode(&t);
  mm_set_matrix(&t);
  mm_set_coordinate(&t);
  mm_set_real(&t);
  mm_set_general(&t);

  int m = A.domain().dim();
  int n = A.range().dim();
  
  TSFMatrixView mv(A);
  int nnz = 0;

  FILE* fp = fopen(filename_.c_str(), "w");

  mm_write_banner(fp, t);

  for(int i=0; i<m; i++)
    {
      TSFArray<int> indices;
      TSFArray<TSFReal> ai;
      mv.getRow(i, indices, ai);
      nnz += indices.length();
    }

  mm_write_mtx_crd_size(fp, m, n, nnz);

  for(int i=0; i<m; i++)
    {
      TSFArray<int> indices;
      TSFArray<TSFReal> ai;
      mv.getRow(i, indices, ai);

      int length = indices.length();
      for(int j=0; j<length; j++)
	{
	  fprintf(fp, "%d %d %f \n", i, indices[j], ai[j]);
	}
    }

  fclose(fp);
}

void MatrixMarketWriter::write(const string& name, const TSFLinearOperator& A) const
{
  write(A);
}
