#include "TSF.h"
#include "MatrixMarketReader.h"
#include "BICGSTABSolver.h"
#include <iostream>
#include <fstream>

int main(int argc, void** argv)
{
  TSFMPI::init(&argc, &argv);

  ofstream fout("fitSine.csv");

  try
    {
      const string MFILE = "matrices/approxSinM7.mtx";
      // const double TOL   = 1.0e-12;

      TSFVectorType vt = new PetraVectorType();

      MatrixMarketReader reader(MFILE);
      TSFLinearOperator M = reader.read(vt);

      TSFVector b = M.range().createMember();
      b[0] = 0.63661977236758;
      b[1] = 0.24960135916919;
      b[2] = 0.13082166937241;
      b[3] = 0.079909504118639;

      // TSFLinearSolver solv = new CGSolver(TOL);
      TSFLinearSolver solv = new BICGSTABSolver();

      TSFLinearOperator MInv = M.inverse(solv);
      TSFVector f = M.inverse(solv) * b;

      cout << "f: " << endl << f << endl;

      fout << " PASSED" << endl;
    }
  catch(exception e)
    {
      cerr << e.what() << endl;
      fout << " EXCEPTION" << endl;
    }

  fout.close();
}
