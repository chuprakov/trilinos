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
void g2pfctn(const Epetra_SerialDenseVector & , Epetra_SerialDenseVector & );
int compproduct(Epetra_SerialDenseVector &, double *, double *);
int compproduct(Epetra_SerialDenseVector &, double *, double *, double *);
double determinant(const Epetra_SerialDenseMatrix &);
int quadrature(const int, const int, Epetra_SerialDenseMatrix &,
               Epetra_SerialDenseVector &);


/**  \brief Performs finite-element assembly of the Hessian of the nonlinear form times an arbitrary vector.

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
  \param  y         [in]  - Reference-counting pointer to the Epetra_MultiVector at which the nonlinear
                            term is evaluated.
  \param  s         [in]  - Reference-counting pointer to the Epetra_MultiVector that is multiplied by the
                            Hessian of the nonlinear form.
  \param  lambda    [in]  - Reference-counting pointer to the Epetra_MultiVector of Lagrange Multipliers.
  \param  hessvec   [out] - Reference-counting pointer to the Epetra_FECrsMatrix containing Hessian of
                            the nonlinear form times vector product.
  \return 0                 if successful.

  \par Detailed Description:

  Assembles the nonlinear term \e hessvec, represented by
     \f[
       \{N''(y,\lambda)s\}_{j} = \langle g''(y_h) \lambda s,\phi_j \rangle
                               = \int_{\Omega} g''(y_h(x)) \lambda(x) s(x) \phi_j(x) dx,
     \f]
  where \f$ g''(y_h) \f$ is given in the local function \e g2pfctn, and \f$\{ \phi_j \}_{j = 1}^{m}\f$ is the
  piecewise linear nodal basis for the state space.
*/
int nonlinhessvec(const Epetra_Comm & Comm,
                  const Epetra_IntSerialDenseVector & ipindx,
                  const Epetra_SerialDenseMatrix & ipcoords,
                  const Epetra_IntSerialDenseVector & pindx,
                  const Epetra_SerialDenseMatrix & pcoords,
                  const Epetra_IntSerialDenseMatrix & t,
                  const Teuchos::RefCountPtr<const Epetra_MultiVector> & y,
                  const Teuchos::RefCountPtr<const Epetra_MultiVector> & s,
                  const Teuchos::RefCountPtr<const Epetra_MultiVector> & lambda,
                  Teuchos::RefCountPtr<Epetra_FEVector> & hessvec)
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

  hessvec = rcp(new Epetra_FEVector(standardmap));

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
  Epetra_SerialDenseVector ls;        // local entries of s
  Epetra_SerialDenseVector Nls;       // N*ls
  Epetra_SerialDenseVector llambda;   // local entries of lambda
  Epetra_SerialDenseVector Nllambda;  // N*llambda
  Epetra_SerialDenseVector lgfctn;    // gfctn(Nly)
  Epetra_SerialDenseVector lgfctnNi;  // lgfctn.*N(:,i)
  Epetra_SerialDenseVector lgfctnNls; // lgfctn.*N(:,i).*llambda.*ls
  Epetra_SerialDenseVector lhessvec;  // local contribution
  // Size and init to zero.
  ly.Size(numNodesPerElem);
  Nly.Size(numQuadPts);
  ls.Size(numNodesPerElem);
  Nls.Size(numQuadPts);
  llambda.Size(numNodesPerElem);
  Nllambda.Size(numQuadPts);
  lgfctn.Size(numQuadPts);
  lgfctnNi.Size(numQuadPts);
  lgfctnNls.Size(numQuadPts);
  lhessvec.Size(numNodesPerElem);
  
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

    // Construct local (to each processor) element view of y, s, l.
    for (j=0; j<numNodesPerElem; j++) {
      ly(j) = (*((*y)(0)))[overlapmap.LID(nodes[j])];
      ls(j) = (*((*s)(0)))[overlapmap.LID(nodes[j])];
      llambda(j) = (*((*lambda)(0)))[overlapmap.LID(nodes[j])];
    }

    Nly.Multiply('N', 'N', 1.0, N, ly, 0.0);            //N*ly
    Nls.Multiply('N', 'N', 1.0, N, ls, 0.0);            //N*ls
    Nllambda.Multiply('N', 'N', 1.0, N, llambda, 0.0);  //N*llambda
    g2pfctn(Nly, lgfctn);
    
    for (int i=0; i<numNodesPerElem; i++) {
      compproduct(lgfctnNi, lgfctn.A(), N[i]);
      compproduct(lgfctnNls, lgfctnNi.A(), Nls.A(), Nllambda.A());  // g2pfctn(Nly).*Nls.*Nllambda.*N(:,i)
      lhessvec(i) = adB*lgfctnNls.Dot(Weights);
    }

    err = hessvec->SumIntoGlobalValues(epetra_nodes, lhessvec);
    if (err<0) return(err);
  }

  // Call global assemble.

  err = hessvec->GlobalAssemble();
  if (err<0) return(err);

  delete [] nodes;

  return(0);
}



/**  \brief Componentwise evaluation of the second derivative of the nonlinear reaction term.
  \param  v   [in]  - Vector at which the second derivative is evaluated.
  \param  gv  [out] - Vector value.
*/
void g2pfctn(const Epetra_SerialDenseVector & v, Epetra_SerialDenseVector & gv) {
  for (int i=0; i<v.M(); i++) {
    gv(i) = 6.0*v(i);
  }  
}

