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
//#include "usr_par.h"
#include <stdlib.h>
#include <algorithm>

using Teuchos::RefCountPtr;
using Teuchos::rcp;
using Teuchos::set_extra_data;

class Epetra_BLAS;
int compproduct(Epetra_SerialDenseVector &, double *, double *);
int compproduct(Epetra_SerialDenseVector &, double *, double *, double *);


/**  \brief Performs finite-element assembly of face mass matrices.

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
                            e(i,3) contains the boundary marker
  \param  B         [out] - Reference-counting pointer to the Epetra_FECrsMatrix describing the FE
                            state/control face mass matrix.
  \param  R         [out] - Reference-counting pointer to the Epetra_FECrsMatrix describing the FE
                            control/control volume mass matrix.
  \return 0                 if successful.

  \par Detailed Description:

  -# Assembles the FE state/control mass matrix \e B, given by
     \f[
       \mathbf{B}_{jk} = b(\mu_k,\phi_j) = - \int_{\partial \Omega}  \mu_k(x) \phi_j(x) dx,
     \f]
     where \f$\{ \phi_j \}_{j = 1}^{m}\f$ is the piecewise linear nodal basis for the finite-dimensional
     state space, and \f$\{ \mu_j \}_{j = 1}^{n}\f$ is the piecewise linear nodal basis for the
     finite-dimensional control space.
  -# Assembles the FE control/control mass matrix \e R, given by
     \f[
       \mathbf{R}_{jk} = b(\mu_k,\mu_j) = - \int_{\partial \Omega}  \mu_k(x) \mu_j(x) dx,
     \f]
     where \f$\{ \mu_j \}_{j = 1}^{n}\f$ is the piecewise linear nodal basis for the finite-dimensional
     control space.
*/
int assemble_bdry(const Epetra_Comm & Comm,
                  const Epetra_IntSerialDenseVector & ipindx,
                  const Epetra_SerialDenseMatrix & ipcoords,
                  const Epetra_IntSerialDenseVector & pindx,
                  const Epetra_SerialDenseMatrix & pcoords,
                  const Epetra_IntSerialDenseMatrix & t,
                  const Epetra_IntSerialDenseMatrix & e,
                  RefCountPtr<Epetra_FECrsMatrix> & B,
                  RefCountPtr<Epetra_FECrsMatrix> & R)
{

  int myPID = Comm.MyPID();
  int numProcs = Comm.NumProc();
  //Usr_Par usr_par;

  int numLocElems = t.M();
  int numLocEdges = e.M();

  int indexBase = 1;

  /* Determine ALL boundary vertices in a subdomain. */
  int * lastindx = 0;
  Epetra_IntSerialDenseVector BdryNode(2*numLocEdges);
  for (int j=0; j<numLocEdges; j++) {
    BdryNode[j] = e(j,0);
    BdryNode[j+numLocEdges] = e(j,1);
  }
  sort(BdryNode.Values(), BdryNode.Values()+2*numLocEdges);
  lastindx  = unique(BdryNode.Values(), BdryNode.Values()+2*numLocEdges);
  const int numBdryNodes = lastindx - BdryNode.Values();
  BdryNode.Resize(numBdryNodes);  

  /* Determine the boundary vertices that belong to this processor. */
  Epetra_IntSerialDenseVector MyBdryNode(numBdryNodes);
  lastindx  = set_intersection(BdryNode.Values(), BdryNode.Values()+numBdryNodes,
                               ipindx.Values(), ipindx.Values()+ipindx.M(),
                               MyBdryNode.Values());
  const int numMyBdryNodes = lastindx - MyBdryNode.Values();
  MyBdryNode.Resize(numMyBdryNodes);

  /* Define data maps. */
  Epetra_Map standardmap(-1, ipindx.M(), const_cast<int*>(ipindx.A()), indexBase, Comm);
  Epetra_Map overlapmap(-1, pindx.M(), const_cast<int*>(pindx.A()), indexBase, Comm);
  Epetra_Map mybdryctrlmap(-1, MyBdryNode.M(), const_cast<int*>(MyBdryNode.A()), indexBase, Comm);

  const int numNodesPerEdge = 2;
  int nodes[numNodesPerEdge];
  int i=0, j=0, err=0;

  /* Define desired matrices. */
  B = rcp(new Epetra_FECrsMatrix(Copy, standardmap, 0));
  R = rcp(new Epetra_FECrsMatrix(Copy, standardmap, 0));
  // NOTE: The data map is the same as for the volume matrices. Later, when
  // FillComplete is called, we will fix the proper domain and range maps. 

  // Declare quantities needed for the call to the local assembly routine.
  int format = Epetra_FECrsMatrix::COLUMN_MAJOR;
  Epetra_IntSerialDenseVector epetra_nodes(View, nodes, numNodesPerEdge);
  Epetra_SerialDenseMatrix vertices(numNodesPerEdge, pcoords.N());

  // Local contribution matrix.
  Epetra_SerialDenseMatrix Bt;

  for(i=0; i<numLocEdges; i++) {
    nodes[0] = e(i,0); nodes[1] = e(i,1);
    for (j=0; j<numNodesPerEdge; j++) {
      vertices(j,0) = pcoords(overlapmap.LID(nodes[j]), 0);
      vertices(j,1) = pcoords(overlapmap.LID(nodes[j]), 1);
    }
    
    double l = sqrt(pow(vertices(0,0)-vertices(1,0),2) + pow(vertices(0,1)-vertices(1,1),2));

    // We have an explicit formula for Bt.
    Bt.Reshape(2,2);
    double sixth = 1.0/6.0;
    Bt(0,0) = l * sixth * 2.0;
    Bt(0,1) = l * sixth * 1.0;
    Bt(1,0) = l * sixth * 1.0;
    Bt(1,1) = l * sixth * 2.0;

    err = B->InsertGlobalValues(epetra_nodes, Bt, format);
    if (err<0) return(err);
    err = R->InsertGlobalValues(epetra_nodes, Bt, format);
    if (err<0) return(err);
  }

  err = B->GlobalAssemble(false);
  err = R->GlobalAssemble(false);

  err = B->FillComplete(mybdryctrlmap, standardmap);
  if (err<0) return(err);
  err = R->FillComplete(mybdryctrlmap, mybdryctrlmap);
  if (err<0) return(err);

  return(0);
}

