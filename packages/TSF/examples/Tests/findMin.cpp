#include "TSF.h"
#include "MatrixMarketReader.h"
#include "MatlabReader.h"
#include "TSFBlockVector.h"
#include "TSFProductSpace.h"
#include "TSFBlockLinearOperator.h"
#include <iostream>
#include <fstream>

const double TOL12 = 1.0e-12;
const double TOL10 = 1.0e-10;

/**
 * Minimize c'x + 0.5x'Qx s.t. Ax = b
 */

TSFVector schurComplement(const TSFLinearOperator &blockA, const TSFVector &blockB);
TSFReal   getError(const TSFVector a, const TSFVector b);

int main(int argc, void** argv)
{
  TSFMPI::init(&argc, &argv);

  ofstream fout("findMin.csv");

  try
    {
      const string AFILE = "matrices/findMinA.mtx";
      const string QFILE = "matrices/tenbyten.m";


      TSFVectorType vt = new DenseSerialVectorType();

      // Read matrices from files
      MatlabReader readQ(QFILE);
      MatrixMarketReader readA(AFILE);
      TSFLinearOperator A = readA.read(vt);
      TSFLinearOperator Q = readQ.read(vt);

      // Create domain and range for blocks
      TSFVectorSpace domain = new TSFProductSpace(Q.domain(), A.adjoint().domain());
      TSFVectorSpace range = new TSFProductSpace(Q.range(), A.range());

      // Create block operator, leaving zero block empty.
      TSFLinearOperator blockA = new TSFBlockLinearOperator(domain, range);
      blockA.setBlock(0, 0, Q);
      blockA.setBlock(0, 1, A.adjoint());
      blockA.setBlock(1, 0, A);

      TSFVector x1 = Q.domain().createMember();
      TSFVector lambda1 = A.adjoint().domain().createMember();
      x1.randomize();
      lambda1.randomize();
      TSFVector blockSoln = domain.createMember();
      blockSoln.setBlock(0, x1);
      blockSoln.setBlock(1, lambda1);

      TSFVector blockB = range.createMember();

      // blockA.apply(blockSoln, blockB);
      TSFVector negC = Q*x1 + A.adjoint()*lambda1;
      TSFVector b = A * x1;
      blockB.setBlock(0, negC);
      blockB.setBlock(1, b);

      TSFVector blockX = schurComplement(blockA, blockB);

      cout << "Error: " << getError(blockX, blockSoln) << endl;

      TSFReal tol = 1.0e-10;
      if(getError(blockX, blockSoln) < tol)
	{
	  fout << "  PASSED" << endl;
	}
      else
	{
	  fout << " FAILED error>tolerance" << endl;
	}
    }
  catch(exception e)
    {
      cerr << e.what() << endl;
      fout << " EXCEPTION" << endl;
    }

  fout.close();
}

TSFVector schurComplement(const TSFLinearOperator &blockA, const TSFVector &blockB)
{
  TSFVector ret(blockA.domain().createMember());

  TSFLinearSolver innerCG = new CGSolver(TOL12);
  TSFLinearSolver outerCG = new CGSolver(TOL10);

  TSFLinearOperator Q = blockA.getBlock(0, 0);
  TSFLinearOperator A = blockA.getBlock(1, 0);
  TSFVector c = -1 * blockB.getBlock(0);
  TSFVector b = blockB.getBlock(1);

  TSFLinearOperator QInv = Q.inverse(innerCG);
  TSFLinearOperator S = -1 * A * QInv * A.adjoint();
  TSFVector QInvC = QInv * c;
  TSFVector lambda = S.inverse(outerCG) * (b + A*QInvC);
  TSFVector x = -1 * (QInvC + QInv*A.adjoint()*lambda);

  ret.setBlock(0, x);
  ret.setBlock(1, lambda);

  return ret;
}

TSFReal   getError(const TSFVector a, const TSFVector b)
{
  TSFVector z = a.getBlock(0);
  TSFVector y = b.getBlock(0);

  return (z-y).norm2();
}
