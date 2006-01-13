#include "Epetra_config.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RefCountPtr.hpp"
#include "includes.h"
#include "usr_par.h"

using Teuchos::RefCountPtr;
using Teuchos::rcp;
using Teuchos::set_extra_data;

class Epetra_BLAS;
int compproduct(Epetra_SerialDenseVector &, double *, double *);
int compproduct(Epetra_SerialDenseVector &, double *, double *, double *);

/**  \brief Performs finite-element assembly of volume stiffness and mass matrices,
            and the right-hand side (forcing term).

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
  \param  e         [in]  - Matrix (EDGE x 3) of edges. \n
                            e(i,1:2) contains the indices of the endpoints of edge i, where i = 1, ..., EDGE \n
                            e(i,3) contains the boundary marker \n
                            e(i,3) = 1  Dirichlet boundary conditions on edge i \n
                            e(i,3) = 2  Neumann boundary conditions on edge i \n
                            e(i,3) = 3  Robin boundary conditions on edge i
  \param  A         [out] - Reference-counting pointer to the Epetra_FECrsMatrix describing the FE volume
                            stiffness matrix for the advection-diffusion equation. Includes advection,
                            diffusion, and reaction terms, and modifications that account for the boundary
                            conditions.
  \param  M         [out] - Reference-counting pointer to the Epetra_FECrsMatrix describing the FE volume
                            mass matrix.
  \param  b         [out] - Reference-counting pointer to the Epetra_FEVector describing the discretized
                            forcing term in the advection-diffusion equation. Includes modifications that
                            account for the boundary conditions.
  \return 0                 if successful.

  \par Detailed Description:

  -# Assembles the FE volume stiffness matrix \e A and the right-hand side \e b for the
  solution of an advection-diffusion equation using piecewise linear finite elements.
  The advection-diffusion equation is 
  \f{align*}
       - \nabla \cdot \left( K(x) \nabla y(x) \right) + \mathbf{c}(x) \cdot \nabla y(x) + r(x) y(x)
       &= f(x), & x &\in \Omega,\;  \\
       y(x) &= d(x),    & x &\in {\partial \Omega}_d, \\
       K(x) \frac{\partial}{\partial \mathbf{n}}  y(x) &= g(x), & x &\in {\partial \Omega}_n, \\
       \sigma_0 K(x) \frac{\partial}{\partial \mathbf{n}} y(x)
       + \sigma_1 y(x) &= g(x), & x &\in {\partial \Omega}_r,
  \f}
  where \f$ K \f$ represents scalar diffusion, \f$ \mathbf{c} \f$ is the advection vector field,
  \f$ r \f$ is the reaction term, \f$ d \f$ and  \f$ g \f$ are given functions, \f$sigma_0\f$ and
  \f$ sigma_1 \f$ are given quantities, \f$ {\partial \Omega}_d \f$ is the Dirichlet boundary,
  \f$ {\partial \Omega}_n \f$ is the Neumann boundary, and \f$ {\partial \Omega}_r \f$ is the Robin
  boundary. The quantities \f$ K \f$, \f$ \mathbf{c} \f$, \f$ r \f$, \f$ d \f$, and \f$ g \f$ are
  assumed to be piecewise linear. Currently, they are to be hard-coded inside this function.
  -# Assembles the FE volume mass matrix \e M.
*/
int assemble(const Epetra_Comm & Comm,
             const Epetra_IntSerialDenseVector & ipindx,
             const Epetra_SerialDenseMatrix & ipcoords,
             const Epetra_IntSerialDenseVector & pindx,
             const Epetra_SerialDenseMatrix & pcoords,
             const Epetra_IntSerialDenseMatrix & t,
             const Epetra_IntSerialDenseMatrix & e,
             RefCountPtr<Epetra_FECrsMatrix> & A,
             RefCountPtr<Epetra_FECrsMatrix> & M,
             RefCountPtr<Epetra_FEVector> & b)
{

  int myPID = Comm.MyPID();
  int numProcs = Comm.NumProc();
  Usr_Par usr_par;

  int numLocElems = t.M();
  int numNodesPerElem = 3;

  int indexBase = 1;

  Epetra_Map standardmap(-1, ipindx.M(), (int*)ipindx.A(), indexBase, Comm);
  Epetra_Map overlapmap(-1, pindx.M(), (int*)pindx.A(), indexBase, Comm);

  int* nodes = new int[numNodesPerElem];
  int i=0, j=0, err=0;

  A = rcp(new Epetra_FECrsMatrix(Copy, standardmap, 0));
  M = rcp(new Epetra_FECrsMatrix(Copy, standardmap, 0));
  b = rcp(new Epetra_FEVector(standardmap));

  // Declare quantities needed for the call to the local assembly routine.
  int format = Epetra_FECrsMatrix::COLUMN_MAJOR;
  Epetra_IntSerialDenseVector epetra_nodes(View, nodes, numNodesPerElem);


  /* ************************  First, we build A and b.  ************************ */
  Epetra_SerialDenseMatrix At;
  Epetra_SerialDenseVector bt;
  Epetra_SerialDenseMatrix vertices(numNodesPerElem, pcoords.N());
  
  Epetra_SerialDenseVector k(numNodesPerElem);
  for (i=0; i< numNodesPerElem; i++) k(i)=1.0;
  Epetra_SerialDenseMatrix c(numNodesPerElem,2);
  for (i=0; i< numNodesPerElem; i++) { c(i,0)=0.0; c(i,1)=0.0; }
  Epetra_SerialDenseVector r(numNodesPerElem);
  for (i=0; i< numNodesPerElem; i++) r(i)=0.0;
  Epetra_SerialDenseVector f(numNodesPerElem);
  for (i=0; i< numNodesPerElem; i++) f(i)=0.0;
  Epetra_SerialDenseVector g(2);
  g(0)=0.0; g(1)=0.0;
  Epetra_SerialDenseVector sigma(2);
  sigma(0)=0.0; sigma(1)=0.0;
  for(i=0; i<numLocElems; i++) {
    nodes[0] = t(i,0); nodes[1] = t(i,1); nodes[2] = t(i,2);
    for (j=0; j<numNodesPerElem; j++) {
      // get vertex
      vertices(j,0) = pcoords(overlapmap.LID(nodes[j]), 0);
      vertices(j,1) = pcoords(overlapmap.LID(nodes[j]), 1);
      // set rhs function
      f(j) = cos(M_PI*vertices(j,0))*cos(M_PI*vertices(j,1)) * 
               (2*M_PI*M_PI + pow(cos(M_PI*vertices(j,0)),2)*pow(cos(M_PI*vertices(j,1)),2) - 1.0);
    }
    lassembly(vertices, k, c, r, f, usr_par, At, bt);
    err = A->InsertGlobalValues(epetra_nodes, At, format);
    if (err<0) return(err);
    err = b->SumIntoGlobalValues(epetra_nodes, bt);
    if (err<0) return(err);
  }

  // MAKE ADJUSTMENTS TO A and b TO ACCOUNT FOR BOUNDARY CONDITIONS.

  // Get Neumann boundary edges.
  Epetra_IntSerialDenseMatrix neumann(e.M(),2);
  j = 0;
  for (i=0; i<e.M(); i++) {
    if (e(i,2)==2) {
      neumann(j,0) = e(i,0);  neumann(j,1) = e(i,1);
      j++;
    }
  }
  neumann.Reshape(j,2);
  // Adjust for Neumann BC's.
  if (neumann.M() != 0) {
    Epetra_BLAS blas;
    Epetra_SerialDenseMatrix quadnodes;
    Epetra_SerialDenseVector quadweights;
    Epetra_SerialDenseMatrix N;
    Epetra_SerialDenseMatrix NN;
    Epetra_SerialDenseVector product(2);

    // Get quadrature weights and points.
    quadrature(1, 2, quadnodes, quadweights);
    
    // Evaluate nodal basis functions at quadrature points
    // N(i,j) value of basis function j at quadrature node i
    N.Shape(quadnodes.M(),2);
    for (i=0; i<quadnodes.M(); i++) {
      N(i,0) = 1.0 - quadnodes(i,0);
      N(i,1) = quadnodes(i,0);
    }

    // Evaluate integrals of 4 products N(i)*N(j).
    NN.Shape(2,2);
    for (i=0; i<2; i++) {
      for (j=0; j<2; j++) {
        compproduct(product, N[i], N[j]);
        NN(i,j) = blas.DOT(quadweights.M(), quadweights.A(), product.A());
      }
    }

    Epetra_IntSerialDenseVector neumann_nodes(2);
    Epetra_SerialDenseVector neumann_add(2);
    double l;
    for (i=0; i<neumann.M(); i++) {
      neumann_nodes(0) = neumann(i,0); neumann_nodes(1) = neumann(i,1);
      neumann_add(0) = pcoords(overlapmap.LID(neumann_nodes(0)),0)
                      -pcoords(overlapmap.LID(neumann_nodes(1)),0);
      neumann_add(1) = pcoords(overlapmap.LID(neumann_nodes(0)),1)
                      -pcoords(overlapmap.LID(neumann_nodes(1)),1);
      l = blas.NRM2(neumann_add.M(), neumann_add.A());
      neumann_add.Multiply('N', 'N', 1.0, NN, g, 0.0);
      neumann_add.Scale(l);
      err = b->SumIntoGlobalValues(neumann_nodes, neumann_add); 
      if (err<0) return(err);
    }
  }

  // Get Robin boundary edges.
  Epetra_IntSerialDenseMatrix robin(e.M(),2);
  j = 0;
  for (i=0; i<e.M(); i++) {
    if (e(i,2)==3) {
      robin(j,0) = e(i,0);  robin(j,1) = e(i,1);
      j++;
    }
  }
  robin.Reshape(j,2);
  // Adjust for Robin BC's.
  if (robin.M() != 0) {
    Epetra_BLAS blas;
    Epetra_SerialDenseMatrix quadnodes;
    Epetra_SerialDenseVector quadweights;
    Epetra_SerialDenseMatrix N;
    Epetra_SerialDenseMatrix NN;
    Epetra_SerialDenseMatrix * NNN;
    Epetra_SerialDenseVector product(2);

    // Get quadrature weights and points.
    quadrature(1, 2, quadnodes, quadweights);
    
    // Evaluate nodal basis functions at quadrature points
    // N(i,j) value of basis function j at quadrature node i
    N.Shape(quadnodes.M(),2);
    for (i=0; i<quadnodes.M(); i++) {
      N(i,0) = 1.0 - quadnodes(i,0);
      N(i,1) = quadnodes(i,0);
    }

    // Evaluate integrals of 4 products N(i)*N(j).
    NN.Shape(2,2);
    for (i=0; i<2; i++) {
      for (j=0; j<2; j++) {
        compproduct(product, N[i], N[j]);
        NN(i,j) = blas.DOT(quadweights.M(), quadweights.A(), product.A());
      }
    }

    // Evaluate integrals of 8 products N(i)*N(j)*N(k).
    NNN = new Epetra_SerialDenseMatrix[2];
    NNN[0].Shape(2,2); NNN[1].Shape(2,2);
    for (i=0; i<2; i++) {
      for (j=0; j<2; j++) {
        for (int k=0; k<2; k++) {
          compproduct(product, N[i], N[j], N[k]);
          NNN[k](i,j) = blas.DOT(quadweights.M(), quadweights.A(),
                                 product.A());
        }
      }
    }

    Epetra_IntSerialDenseVector robin_nodes(2);
    Epetra_SerialDenseVector robin_b_add(2);
    Epetra_SerialDenseMatrix robin_A_add(2,2);
    double l;
    for (i=0; i<robin.M(); i++) {
      robin_nodes(0) = robin(i,0); robin_nodes(1) = robin(i,1);
      
      robin_b_add(0) = pcoords(overlapmap.LID(robin_nodes(0)),0)
                      -pcoords(overlapmap.LID(robin_nodes(1)),0);
      robin_b_add(1) = pcoords(overlapmap.LID(robin_nodes(0)),1)
                      -pcoords(overlapmap.LID(robin_nodes(1)),1);
      l = blas.NRM2(robin_b_add.M(), robin_b_add.A());
      robin_b_add.Multiply('N', 'N', 1.0, NN, g, 0.0);
      robin_b_add.Scale(l);
      err = b->SumIntoGlobalValues(robin_nodes, robin_b_add); 
      if (err<0) return(err);

      NNN[0].Scale(sigma(0)); NNN[1].Scale(sigma(1));
      robin_A_add += NNN[0]; robin_A_add += NNN[1];
      robin_A_add.Scale(l);
      err = A->InsertGlobalValues(robin_nodes, robin_A_add, format);
      if (err<0) return(err);
    }

    delete [] NNN;
  }

  // Get Dirichlet boundary edges.
  Epetra_IntSerialDenseMatrix dirichlet(e.M(),2);
  j = 0;
  for (i=0; i<e.M(); i++) {
    if (e(i,2)==1) {
      dirichlet(j,0) = e(i,0);  dirichlet(j,1) = e(i,1);
      j++;
    }
  }
  dirichlet.Reshape(j,2);
  // DIRICHLET NOT DONE! DO THIS LATER!!!!

  /* ************************  Done building A and b.  ************************ */



  /* ************************  Second, we build M.  ************************ */

  Epetra_SerialDenseMatrix Mt;

  for (i=0; i< numNodesPerElem; i++) k(i)=0.0;
  for (i=0; i< numNodesPerElem; i++) { c(i,0)=0.0; c(i,1)=0.0; }
  for (i=0; i< numNodesPerElem; i++) r(i)=1.0;
  for (i=0; i< numNodesPerElem; i++) f(i)=0.0;
  g(0)=0.0; g(1)=0.0;
  sigma(0)=0.0; sigma(1)=0.0;
  for(i=0; i<numLocElems; i++) {
    nodes[0] = t(i,0); nodes[1] = t(i,1); nodes[2] = t(i,2);
    for (j=0; j<numNodesPerElem; j++) {
      vertices(j,0) = pcoords(overlapmap.LID(nodes[j]), 0);
      vertices(j,1) = pcoords(overlapmap.LID(nodes[j]), 1);
    }
    lassembly(vertices, k, c, r, f, usr_par, Mt, bt);
    err = M->InsertGlobalValues(epetra_nodes, Mt, format);
    if (err<0) return(err);
  }

  /* ************************  Done building M.  ************************ */



  // Call global assemble and FillComplete at the same time.

  err = A->GlobalAssemble();
  if (err<0) return(err);

  err = b->GlobalAssemble();
  if (err<0) return(err);

  err = M->GlobalAssemble();
  if (err<0) return(err);

  delete [] nodes;

  return(0);
}

