#include "Epetra_config.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_LAPACK.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RefCountPtr.hpp"
// #include "includes.h"
// #include "usr_par.h"

using Teuchos::RefCountPtr;
using Teuchos::rcp;

class Epetra_BLAS;
void gfctn(const Epetra_SerialDenseVector & , Epetra_SerialDenseVector & );
int compproduct(Epetra_SerialDenseVector &, double *, double *);
int compproduct(Epetra_SerialDenseVector &, double *, double *, double *);
double determinant(const Epetra_SerialDenseMatrix &);
int quadrature(const int, const int, Epetra_SerialDenseMatrix &,
               Epetra_SerialDenseVector &);


/**  \brief Performs finite-element assembly of the nonlinear reaction term.

  \param  Comm      [in]  - The Epetra (MPI) communicator.
  \param  ipindx    [in]  - Vector of NUMIP indices of nodes that are `unique' to a subdomain
                            (i.e. owned by the corresponding processor).
  \param  ipcoords  [in]  - Matrix (NUMIP x 2) of x-y-coordinates of the vertices cooresponding
                            to indices ipindx: \n
                            ipcoords(i,0) x-coordinate of node i, \n
                            ipcoords(i,1) y-coordinate of node i.
  \param  pindx     [in]  - Vector of NUMP indices of ALL nodes in a subdomain, including
                            the shared nodes.
  \param  pcoords   [in]  - Matrix (NUMP x 2) of x-y-coordinates of the vertices cooresponding
                            to indices pindx: \n
                            pcoords(i,0) x-coordinate of node i, \n
                            pcoords(i,1) y-coordinate of node i.
  \param  t         [in]  - Matrix (ELE x 3) of indices of the vertices in a triangle: \n
                            t(i,j) index of the j-th vertex in triangle i, where i = 1, ..., ELE
  \param  y         [out] - Reference-counting pointer to the Epetra_MultiVector at which the nonlinear
                            form is evaluated.
  \param  g         [out] - Reference-counting pointer to the Epetra_FEVector containing the value
                            of the nonlinear form.
  \return 0                 if successful.

  \par Detailed Description:

  Assembles the nonlinear term \e g, represented by
     \f[
       \{N(y)\}_{j} = \langle g(y_h),\phi_j \rangle =  \int_{\Omega} g(y_h(x)) \phi_j(x) dx,
     \f]
  where \f$ g(y_h) \f$ is given in the local function \e gfctn, and \f$\{ \phi_j \}_{j = 1}^{m}\f$ is the
  piecewise linear nodal basis for the state space.
*/
int nonlinvec(const Epetra_Comm & Comm,
              const Epetra_IntSerialDenseVector & ipindx,
              const Epetra_SerialDenseMatrix & ipcoords,
              const Epetra_IntSerialDenseVector & pindx,
              const Epetra_SerialDenseMatrix & pcoords,
              const Epetra_IntSerialDenseMatrix & t,
              const Teuchos::RefCountPtr<const Epetra_MultiVector> & y,
              Teuchos::RefCountPtr<Epetra_FEVector> & g)
{

  int myPID = Comm.MyPID();
  int numProcs = Comm.NumProc();

  int numLocNodes     = pindx.M();
  int numMyLocNodes   = ipindx.M();
  int numLocElems     = t.M();
  int numNodesPerElem = 3;

  int indexBase = 1;

  Epetra_Map standardmap(-1, numMyLocNodes, (int*)ipindx.A(), indexBase, Comm);
  Epetra_Map overlapmap(-1, numLocNodes, (int*)pindx.A(), indexBase, Comm);

  g = rcp(new Epetra_FEVector(standardmap));

  int* nodes = new int[numNodesPerElem];
  int i=0, j=0, err=0;
  
  // get quadrature nodes and weights
  Epetra_SerialDenseMatrix Nodes;
  Epetra_SerialDenseVector Weights;
  quadrature(2,3,Nodes,Weights);
  int numQuadPts = Nodes.M();

  // Evaluate nodal basis functions and their derivatives at quadrature points
  // N(i,j) = value of the j-th basis function at quadrature node i.
  Epetra_SerialDenseMatrix N;
  N.Shape(numQuadPts,3);
  for (int i=0; i<numQuadPts; i++) {
    N(i,0) = 1.0 - Nodes(i,0) - Nodes(i,1);
    N(i,1) = Nodes(i,0);
    N(i,2) = Nodes(i,1);
  }

  // Declare quantities needed for the call to the local assembly routine.
  Epetra_IntSerialDenseVector epetra_nodes(View, nodes, numNodesPerElem);
  Epetra_SerialDenseMatrix vertices(numNodesPerElem, pcoords.N());

  Epetra_SerialDenseVector ly;        // local entries of y
  Epetra_SerialDenseVector Nly;       // N*ly
  Epetra_SerialDenseVector lgfctn;    // gfctn(Nly)
  Epetra_SerialDenseVector lgfctnNi;  // lgfctn.*N(:,i)
  Epetra_SerialDenseVector lg;        // local contribution
  // Size and init to zero.
  ly.Size(numNodesPerElem);
  Nly.Size(numQuadPts);
  lgfctn.Size(numQuadPts);
  lgfctnNi.Size(numQuadPts);
  lg.Size(numNodesPerElem);
  
  Epetra_SerialDenseMatrix B(2,2);
  double adB;
  
  for(i=0; i<numLocElems; i++) {

    nodes[0] = t(i,0); nodes[1] = t(i,1); nodes[2] = t(i,2);
    for (j=0; j<numNodesPerElem; j++) {
      vertices(j,0) = pcoords(overlapmap.LID(nodes[j]), 0);
      vertices(j,1) = pcoords(overlapmap.LID(nodes[j]), 1);
    }

    // Construct affine transformation matrix.
    for(int i=0; i<2; i++) {
      B(i,0) = vertices(1,i)-vertices(0,i);
      B(i,1) = vertices(2,i)-vertices(0,i);
    }
    adB  = abs(determinant(B));

    // Construct local (to each processor) element view of y. 
    for (j=0; j<numNodesPerElem; j++) {
      ly(j) = (*((*y)(0)))[overlapmap.LID(nodes[j])];
    }

    Nly.Multiply('N', 'N', 1.0, N, ly, 0.0);
    gfctn(Nly, lgfctn);
    
    for (int i=0; i<numNodesPerElem; i++) {
      compproduct(lgfctnNi, lgfctn.A(), N[i]);
      lg(i) = adB*lgfctnNi.Dot(Weights);
    }
    
    err = g->SumIntoGlobalValues(epetra_nodes, lg);
    if (err<0) return(err);
  }

  // Call global assemble.

  err = g->GlobalAssemble();
  if (err<0) return(err);

  delete [] nodes;

  return(0);
}


/**  \brief Componentwise evaluation of the nonlinear reaction term.
  \param  v   [in]  - Vector at which the nonlinear function is evaluated.
  \param  gv  [out] - Vector value.
*/
void gfctn(const Epetra_SerialDenseVector & v, Epetra_SerialDenseVector & gv) {
  for (int i=0; i<v.M(); i++) {
    gv(i) = pow(v(i),3)-v(i);
  }  
}

