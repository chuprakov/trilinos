#include "MatrixMarketReader.h"
#include "TSFError.h"

#include "TSFMatrixOperator.h"
#include "TSFMatrixView.h"

using namespace TSF;

#include "mmio.h"

MatrixMarketReader::MatrixMarketReader(const string& filename)
	: TSFMatrixReaderBase(filename)
{}


TSFLinearOperator MatrixMarketReader::read(const TSFVectorType& vectorType) const
{
	FILE* fp = fopen(filename_.c_str(), "r");
	
	MM_typecode mmType;

	int ierr = mm_read_banner(fp, &mmType);

	switch(ierr)
		{
		case MM_PREMATURE_EOF:
			TSFError::raise("unexpected EOF in MatrixMarketReader::read");
			break;
		case MM_NO_HEADER:
			TSFError::raise("bad header in MatrixMarketReader::read");
			break;
		case MM_UNSUPPORTED_TYPE:
			TSFError::raise("unsupported type in MatrixMarketReader::read");
			break;
		default:
			break;
		}

	if (!mm_is_real(mmType))
		{
			TSFError::raise("attempted to read non-real matrix in MatrixMarketReader::read");
		}

	TSFLinearOperator rtn;

	MMStructure structure = MMGeneral;
	if (mm_is_general(mmType)) structure = MMGeneral;
	else if (mm_is_symmetric(mmType)) structure = MMSymmetric;
	else if (mm_is_skew(mmType)) structure = MMSkew;
	else TSFError::raise("unknown structure in MatrixMarketReader::read");

	if (mm_is_coordinate(mmType))
		{
			int nRows;
			int nCols;
			int nnz;
			int ierr = mm_read_mtx_crd_size(fp, &nRows, &nCols, &nnz);
			switch(ierr)
				{
				case MM_PREMATURE_EOF:
					TSFError::raise("unexpected EOF in MatrixMarketReader::read");
					break;
				default:
					break;
				}
			TSFVectorSpace rowSpace = vectorType.createSpace(nRows);
			TSFVectorSpace colSpace = vectorType.createSpace(nCols);

			rtn = vectorType.createMatrix(colSpace, rowSpace);

			TSFMatrixView A(rtn);
			readSparse(fp, A, nRows, structure, nnz);
		}
	else
		{
			int nRows;
			int nCols;
			int ierr = mm_read_mtx_array_size(fp, &nRows, &nCols);
			switch(ierr)
				{
				case MM_PREMATURE_EOF:
					TSFError::raise("unexpected EOF in MatrixMarketReader::read");
					break;
				default:
					break;
				}
			TSFVectorSpace rowSpace = vectorType.createSpace(nRows);
			TSFVectorSpace colSpace = vectorType.createSpace(nCols);

			rtn = vectorType.createMatrix(colSpace, rowSpace);

			TSFMatrixView A(rtn);
			readDense(fp, A, structure);
		}

	return rtn;
		
}

void MatrixMarketReader::readSparse(FILE* fp, TSFMatrixView& matrix, int nRows,
																			 MMStructure structure, int nnz) const
{
	TSFArray<int> bandwidth(nRows, 0);
	TSFArray<TSFArray<int> > colIndices(nRows);
	TSFArray<TSFArray<double> > aij(nRows);

	for (int k=0; k<nnz; k++)
		{
			int i;
			int j;
			double val;
			fscanf(fp, "%d %d %lg", &i, &j, &val);
			/* convert from fortran to C indexing */
			i--;
			j--;
			/* put the indices and values into a temporary sparse matrix */
			bandwidth[i]++;
			colIndices[i].append(j);
			aij[i].append(val);
			/* add elements implied by symmetry */
			switch(structure)
				{
				case MMGeneral:
					break;
				case MMSymmetric:
					if (i != j)
						{
						        ++k;
							bandwidth[j]++;
							colIndices[j].append(i);
							aij[j].append(val);
						}
					break;
				case MMSkew:
					bandwidth[j]++;
					colIndices[j].append(i);
					aij[j].append(-val);
				}
		}

	if (matrix.requiresBandwidth())
		{
			matrix.setBandwidth(bandwidth.size(), &(bandwidth[0]));
		}

	matrix.freezeStructure();

	for (int row=0; row<nRows; row++)
		{
			matrix.setRowStructure(row, bandwidth[row], &(colIndices[row][0]));
			matrix.addToRow(row, bandwidth[row], &(colIndices[row][0]), &(aij[row][0]));
		}

	matrix.freezeValues();
}

void MatrixMarketReader::readDense(FILE* /* fp */, TSFMatrixView& /* matrix */, 
																	 MMStructure /* structure */) const 
{
	TSFError::raise("MatrixMarketReader::readDense not yet implemented");
}

