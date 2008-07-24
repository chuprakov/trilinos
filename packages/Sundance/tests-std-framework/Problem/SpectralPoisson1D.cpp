/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "Sundance.hpp"


/** 
 * Solves the Poisson equation in 1D
 */


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})

  int main(int argc, char** argv)
{
  try
  {
    Sundance::init(&argc, &argv);
    int np = MPIComm::world().getNProc();

    /* We will do our linear algebra using Epetra */
    VectorType<double> vecType = new EpetraVectorType();

    /* Create a mesh. It will be of type BasisSimplicialMesh, and will
     * be built using a PartitionedLineMesher. */
    MeshType meshType = new BasicSimplicialMeshType();
    MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10, meshType);
    Mesh mesh = mesher.getMesh();

    /* Create a cell filter that will identify the maximal cells
     * in the interior of the domain */
    CellFilter interior = new MaximalCellFilter();
    CellFilter points = new DimensionalCellFilter(0);
    CellFilter leftPoint = points.subset(new LeftPointTest());

    /* Create the Spectral Basis */
    int ndim = 2;
    int order = 2;

    SpectralBasis sbasis = new HermiteSpectralBasis(ndim, order); 

    cout << "created the spectral basis" << endl;

    /* Create the Spectral Unknown and test functions */

    Expr u = new UnknownFunction(new Lagrange(2),sbasis, "u");
    Expr v = new TestFunction(new Lagrange(2), sbasis, "v");

    cout <<"Unknown and Test Functions " << endl; 

    /* Create differential operator and coordinate function */
    Expr dx = new Derivative(0);
    Expr x = new CoordExpr(0);

    /* We need a quadrature rule for doing the integrations */
    QuadratureFamily quad = new GaussianQuadrature(2);


    /* Define the Stochastic RHS */

    Array<Expr> Coeff(sbasis.nterms());  /* array of coefficients for the spectra expression */
    Coeff[0] = 2.0;
    for(int i=1; i<sbasis.nterms(); i++)
    {

      cout << i<< endl;
      Coeff[i] = 0.0;
    }

    Expr F = new SpectralExpr(sbasis, Coeff);

    cout << "created F" << endl ;

    /* Define the weak form */
    Expr eqn = Integral(interior, -(dx*v)*(dx*u) - F*v, quad);
    /* Define the Dirichlet BC */
    Expr bc = EssentialBC(leftPoint, v*u, quad);

      
    cout << "done eq and bc " << endl;


    /* We can now set up the linear problem! */

    LinearProblem prob(mesh, eqn, bc, v, u, vecType); 

    ParameterXMLFileReader reader(searchForFile("SolverParameters/bicgstab.xml"));
    ParameterList solverParams = reader.getParameters();
    cout << "params = " << solverParams << endl;


    LinearSolver<double> solver 
      = LinearSolverBuilder::createSolver(solverParams);

    Expr soln = prob.solve(solver);

    TEST_FOR_EXCEPTION(!(prob.solveStatus().finalState() == SolveConverged),
      RuntimeError,
      "solve failed");

    Expr zz = Integral(interior, soln[0], quad);
    double totalU = evaluateIntegral(mesh, zz);
    cerr << "integral(u[0]) = " << totalU << endl;

    Expr exactSoln = x*(x-2.0);
    Expr err = pow(soln[0] - exactSoln, 2.0);
    for (int i=1; i<sbasis.nterms(); i++)
    {
      err = err + pow(soln[i], 2.0);
    }
        

    Expr errExpr = Integral(interior, 
      err,
      new GaussianQuadrature(4));

    Expr derivErr = dx*(exactSoln-soln[0]);
    Expr derivErrExpr = Integral(interior, 
      pow(dx*(soln[0]-exactSoln), 2),
      new GaussianQuadrature(2));

    double errorSq = evaluateIntegral(mesh, errExpr);
    cout << "error norm = " << sqrt(errorSq) << endl << endl;



    double derivErrorSq = evaluateIntegral(mesh, derivErrExpr);
    cout << "deriv error norm = " << sqrt(derivErrorSq) << endl << endl;

    double tol = 1.0e-12;
    Sundance::passFailTest(sqrt(errorSq + derivErrorSq), tol);
  }
	catch(exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); return Sundance::testStatus(); 
}
