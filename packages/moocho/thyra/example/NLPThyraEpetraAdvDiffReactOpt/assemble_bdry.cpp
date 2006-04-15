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
#include "Teuchos_VerboseObject.hpp"
#include "includes.h"
//#include "usr_par.h"
#include <stdlib.h>
#include <algorithm>

class Epetra_BLAS;
int compproduct(Epetra_SerialDenseVector &, double *, double *);
int compproduct(Epetra_SerialDenseVector &, double *, double *, double *);

//#define GLPAPP_SHOW_BOUNDARY_ASSEMBLY

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
                            t(i,j) index of the 0-based j-th vertex in triangle i, where i = 0, ..., numElements-1
  \param  e         [in]  - Matrix (EDGE x 3) of edges. \n
                            e(i,1:2) contains the indices of the endpoints of edge i, where i = 0, ..., numEdges-1 \n
                            e(i,3) contains the boundary marker
  \param  B_out     [out] - Reference-counting pointer to the Epetra_FECrsMatrix describing the FE
                            state/control face mass matrix.
  \param  R_out     [out] - Reference-counting pointer to the Epetra_FECrsMatrix describing the FE
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
int assemble_bdry(
  const Epetra_Comm                                &Comm
  ,const Epetra_IntSerialDenseVector               &ipindx
  ,const Epetra_SerialDenseMatrix                  &ipcoords
  ,const Epetra_IntSerialDenseVector               &pindx
  ,const Epetra_SerialDenseMatrix                  &pcoords
  ,const Epetra_IntSerialDenseMatrix               &t
  ,const Epetra_IntSerialDenseMatrix               &e
  ,Teuchos::RefCountPtr<Epetra_FECrsMatrix>        *B_out
  ,Teuchos::RefCountPtr<Epetra_FECrsMatrix>        *R_out
  )
{

  using Teuchos::rcp;

#ifdef GLPAPP_SHOW_BOUNDARY_ASSEMBLY
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  *out << "\nEntering assemble_bdry(...) ...\n";
#endif

  int numLocElems = t.M();
  int numLocEdges = e.M();

  int indexBase = 1;

  //
  // Create a sorted (by global ID) list of boundry nodes
  //
  int * lastindx = 0;
  Epetra_IntSerialDenseVector BdryNode(2*numLocEdges);
  for (int j=0; j<numLocEdges; j++) {
    BdryNode[j] = e(j,0);
    BdryNode[j+numLocEdges] = e(j,1);
  }
  std::sort(BdryNode.Values(), BdryNode.Values()+2*numLocEdges);
  lastindx  = std::unique(BdryNode.Values(), BdryNode.Values()+2*numLocEdges);
  const int numBdryNodes = lastindx - BdryNode.Values();
  BdryNode.Resize(numBdryNodes); // RAB: This does not overwrite?

  //
  // Determine the boundary vertices that belong to this processor.
  //
  Epetra_IntSerialDenseVector MyBdryNode(numBdryNodes);
  lastindx = std::set_intersection(
    BdryNode.Values(), BdryNode.Values()+numBdryNodes,  // (Sorted) set A
    ipindx.Values(), ipindx.Values()+ipindx.M(),        // (Sorted) set B
    MyBdryNode.Values()                                 // Intersection
    );
  const int numMyBdryNodes = lastindx - MyBdryNode.Values();
  MyBdryNode.Resize(numMyBdryNodes); // RAB: This does not overwrite?
  
  //
  // Define the maps for the various lists
  //
  Epetra_Map standardmap(-1, ipindx.M(), const_cast<int*>(ipindx.A()), indexBase, Comm);
  Epetra_Map overlapmap(-1, pindx.M(), const_cast<int*>(pindx.A()), indexBase, Comm);
  Epetra_Map mybdryctrlmap(-1, numMyBdryNodes, const_cast<int*>(MyBdryNode.A()), indexBase, Comm);
  // Above, it is important to note what mybndyctrlmap represents.  It is the
  // a sorted list of global vertex node IDS for nodes on a boundary that are
  // uniquely owned by the local process.

#ifdef GLPAPP_SHOW_BOUNDARY_ASSEMBLY
  *out << "\nstandardmap:\n";
  standardmap.Print(*Teuchos::OSTab(out).getOStream());
  *out << "\nmybdryctrlmap:\n";
  mybdryctrlmap.Print(*Teuchos::OSTab(out).getOStream());
#endif

  //
  // Allocate matrices to fill
  //
  Teuchos::RefCountPtr<Epetra_FECrsMatrix>
    B = rcp(new Epetra_FECrsMatrix(Copy,standardmap,0)),
    R = rcp(new Epetra_FECrsMatrix(Copy,standardmap,0));
  // NOTE: The data map is the same as for the volume matrices. Later, when
  // FillComplete is called, we will fix the proper domain and range maps. 

  // Declare quantities needed for the call to the local assembly routine.
  const int numNodesPerEdge = 2;
  Epetra_IntSerialDenseVector epetra_nodes(numNodesPerEdge);

  //
  // Load B and R by looping through the edges
  //

  Epetra_SerialDenseMatrix Bt(2,2);
  int err=0;

  for( int i=0; i < numLocEdges; i++ ) {

    const int
      global_id_0 = e(i,0),
      global_id_1 = e(i,1),
      local_id_0  = overlapmap.LID(global_id_0), // O(log(numip)) bindary search
      local_id_1  = overlapmap.LID(global_id_1); // O(log(numip)) bindary search

    epetra_nodes(0) = global_id_0;
    epetra_nodes(1) = global_id_1;

    const double
      x0 = pcoords(local_id_0,0),
      y0 = pcoords(local_id_0,1),
      x1 = pcoords(local_id_1,0),
      y1 = pcoords(local_id_1,1);
    
    const double l = sqrt(pow(x0-x1,2) + pow(y0-y1,2));  // Length of this edge
    
    // We have an explicit formula for Bt.
    const double l_sixth = l * (1.0/6.0);
    Bt(0,0) = l_sixth * 2.0;
    Bt(0,1) = l_sixth * 1.0;
    Bt(1,0) = l_sixth * 1.0;
    Bt(1,1) = l_sixth * 2.0;

#ifdef GLPAPP_SHOW_BOUNDARY_ASSEMBLY
    *out
      << "\nEdge global nodes = ("<<global_id_0<<","<<global_id_1<<"),"
      << " local nodes = ("<<local_id_0<<","<<local_id_1<<"),"
      << " Bt = ["<<Bt(0,0)<<","<<Bt(0,1)<<";"<<Bt(1,0)<<","<<Bt(1,1)<<"]\n";
#endif


    const int format = Epetra_FECrsMatrix::COLUMN_MAJOR;
    err = B->InsertGlobalValues(epetra_nodes,Bt,format);
    if (err<0) return(err);
    err = R->InsertGlobalValues(epetra_nodes,Bt,format);
    if (err<0) return(err);
    
  }

/*

  err = B->GlobalAssemble(false);
  if (err<0) return(err);
  err = R->GlobalAssemble(false);
  if (err<0) return(err);

  err = B->FillComplete(mybdryctrlmap,standardmap);
  if (err<0) return(err);
  err = R->FillComplete(mybdryctrlmap,mybdryctrlmap);
  if (err<0) return(err);

*/

  err = B->GlobalAssemble(mybdryctrlmap,standardmap,true);
  if (err<0) return(err);
  err = R->GlobalAssemble(mybdryctrlmap,mybdryctrlmap,true);
  if (err<0) return(err);

  if(B_out) *B_out = B;
  if(R_out) *R_out = R;

#ifdef GLPAPP_SHOW_BOUNDARY_ASSEMBLY
  *out << "B =\n";
  B->Print(*Teuchos::OSTab(out).getOStream());
  *out << "\nLeaving assemble_bdry(...) ...\n";
#endif

  return(0);

}

