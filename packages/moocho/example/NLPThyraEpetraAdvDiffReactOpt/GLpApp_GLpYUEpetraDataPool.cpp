#include "GLpApp_GLpYUEpetraDataPool.hpp"
//#include "Amesos.h"
//#include "Ifpack.h"
#include "Amesos_Klu.h"
//#include "Amesos_Umfpack.h"
//#include "Amesos_Superludist.h"
#include "Teuchos_ParameterList.hpp"
#include "EpetraExt_Reindex_LinearProblem.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_Transpose_RowMatrix.h"
//#include "GLpApp_SchurOp.hpp"

namespace GLpApp {

GLpYUEpetraDataPool::GLpYUEpetraDataPool( Teuchos::RefCountPtr<const Epetra_Comm> const& commptr, double beta, const char myfile[] )
  : commptr_(commptr),
    beta_(beta)
{
  ipcoords_ = Teuchos::rcp( new Epetra_SerialDenseMatrix() );
  ipindx_ = Teuchos::rcp( new Epetra_IntSerialDenseVector() );
  pcoords_ = Teuchos::rcp( new Epetra_SerialDenseMatrix() );
  pindx_ = Teuchos::rcp( new Epetra_IntSerialDenseVector() );
  t_ = Teuchos::rcp( new Epetra_IntSerialDenseMatrix() );
  e_ = Teuchos::rcp( new Epetra_IntSerialDenseMatrix() );
  
  // Read subdomain info.
  meshreader(*commptr_, *ipindx_, *ipcoords_, *pindx_, *pcoords_, *t_, *e_, myfile);

  // Assemble volume and boundary mass and stiffness matrices, and the right-hand side of the PDE.
  assemble(*commptr, *ipindx_, *ipcoords_, *pindx_, *pcoords_, *t_, *e_, A_, H_, b_);
  assemble_bdry(*commptr, *ipindx_, *ipcoords_, *pindx_, *pcoords_, *t_, *e_, B_, R_);

  // Set desired state q.
  Epetra_Map standardmap(A_->DomainMap());
  q_ = Teuchos::rcp(new Epetra_FEVector(standardmap));
  int * qintvalues = new int[standardmap.NumMyElements()];
  double * qdvalues = new double[standardmap.NumMyElements()];
  standardmap.MyGlobalElements(qintvalues);
  for (int i = 0; i < standardmap.NumMyElements(); i++)
      qdvalues[i]=cos( M_PI* ((*ipcoords_)(i,0)) ) * cos( M_PI* ((*ipcoords_)(i,1)) );
  q_->ReplaceGlobalValues(standardmap.NumMyElements(), qintvalues, qdvalues);
  q_->GlobalAssemble();
}

void GLpYUEpetraDataPool::computeAll( const GenSQP::Vector &x )
{

  // Dynamic cast back to Epetra vectors here.
  Teuchos::RefCountPtr<const Epetra_MultiVector> ey =
        (Teuchos::dyn_cast<GenSQP::YUEpetraVector>(const_cast<GenSQP::Vector&>(x))).getYVector();

  computeNy(ey);

  computeNpy(ey);

  computeAugmat();
  
}


int GLpYUEpetraDataPool::solveAugsys( const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsy,
                                       const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsu,
                                       const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsp,
                                       const Teuchos::RefCountPtr<Epetra_MultiVector> & y,
                                       const Teuchos::RefCountPtr<Epetra_MultiVector> & u,
                                       const Teuchos::RefCountPtr<Epetra_MultiVector> & p,
                                       double tol )
{
/*
  int systemChoice = 1;   // 1 for full KKT system solve, 2 for Schur complement solve
  int solverChoice = 12;  // These options are for the full KKT system solve.
                          // 11 for AztecOO with built-in Schwarz DD preconditioning and ILU on subdomains
                          // 12 for AztecOO with IFPACK Schwarz DD preconditioning and Umfpack on subdomains
                          // 13 for a direct sparse solver (Umfpack, KLU)
  
  if (systemChoice == 1) {
    // We're using the full KKT system formulation to solve the augmented system.
   
    Epetra_Map standardmap(A_->DomainMap());
    int numstates = standardmap.NumGlobalElements();
    Epetra_Map bdryctrlmap(B_->DomainMap());
    int numcontrols = bdryctrlmap.NumGlobalElements();
    Epetra_Vector rhs( (Epetra_BlockMap&)Augmat_->RangeMap() );
    Epetra_Vector soln( (Epetra_BlockMap&)Augmat_->RangeMap() );
    soln.Random();  

    std::vector<double> values(rhsy->MyLength() + rhsu->MyLength() + rhsp->MyLength());
    std::vector<int> indices(rhsy->MyLength() + rhsu->MyLength() + rhsp->MyLength());
    ((Epetra_BlockMap&)Augmat_->RangeMap()).MyGlobalElements(&indices[0]);

    for (int i=0; i<rhsy->MyLength(); i++) {
      values[i] = (*((*rhsy)(0)))[i];
    }
    for (int i=0; i<rhsu->MyLength(); i++) {
      values[i+rhsy->MyLength()] = (*((*rhsu)(0)))[i];
    }
    for (int i=0; i<rhsp->MyLength(); i++) {
      values[i+rhsy->MyLength()+rhsu->MyLength()] = (*((*rhsp)(0)))[i];
    }

    rhs.ReplaceGlobalValues(rhsy->MyLength() + rhsu->MyLength() + rhsp->MyLength(), &values[0], &indices[0]);

    if (solverChoice == 11) {
      int Overlap = 3;
      int ival = 4;

      AztecOO::AztecOO kktsolver(&(*Augmat_), &soln, &rhs);
      kktsolver.SetAztecOption( AZ_solver, AZ_gmres );
      kktsolver.SetAztecOption( AZ_precond, AZ_dom_decomp );
      kktsolver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
      //kktsolver.SetAztecOption( AZ_kspace, 2*numstates+numcontrols );
      kktsolver.SetAztecOption( AZ_kspace, 9000 );
      kktsolver.SetAztecOption(AZ_overlap,Overlap);
      kktsolver.SetAztecOption(AZ_graph_fill,ival);
      //kktsolver.SetAztecOption(AZ_poly_ord, ival);
      //kktsolver.SetAztecParam(AZ_drop, 1e-9);
      kktsolver.SetAztecParam(AZ_athresh, 1e-5);
      //kktsolver.SetAztecParam(AZ_rthresh, 0.0);
      kktsolver.SetAztecOption( AZ_reorder, 0 );
      //kktsolver.SetAztecParam44( AZ_ilut_fill, 1.5 );
      kktsolver.SetAztecOption( AZ_output, AZ_last );
      //kktsolver.Iterate(2*numstates+numcontrols,1e-12);
      kktsolver.Iterate(9000,1e-11);
      //cout << soln;
    }
    else if (solverChoice == 12) {
      // =============================================================== //
      // B E G I N N I N G   O F   I F P A C K   C O N S T R U C T I O N //
      // =============================================================== //

      Teuchos::ParameterList List;

      // allocates an IFPACK factory. No data is associated
      // to this object (only method Create()).
      Ifpack Factory;

      // create the preconditioner. For valid PrecType values,
      // please check the documentation
      string PrecType = "Amesos";
      int OverlapLevel = 2; // must be >= 0. If Comm.NumProc() == 1,
                            // it is ignored.
  
      Ifpack_Preconditioner* Prec = Factory.Create(PrecType, &(*Augmat_), OverlapLevel);
      assert(Prec != 0);

      // specify the Amesos solver to be used.
      // If the selected solver is not available,
      // IFPACK will try to use Amesos' KLU (which is usually always
      // compiled). Amesos' serial solvers are:
      // "Amesos_Klu", "Amesos_Umfpack", "Amesos_Superlu"
      List.set("amesos: solver type", "Amesos_Umfpack");

      // sets the parameters
      IFPACK_CHK_ERR(Prec->SetParameters(List));

      // initialize the preconditioner. At this point the matrix must
      // have been FillComplete()'d, but actual values are ignored.
      // At this call, Amesos will perform the symbolic factorization.
      IFPACK_CHK_ERR(Prec->Initialize());

      // Builds the preconditioners, by looking for the values of
      // the matrix. At this call, Amesos will perform the
      // numeric factorization.
      IFPACK_CHK_ERR(Prec->Compute());

      // =================================================== //
      // E N D   O F   I F P A C K   C O N S T R U C T I O N //
      // =================================================== //

      // need an Epetra_LinearProblem to define AztecOO solver
      Epetra_LinearProblem Problem;
      Problem.SetOperator(&(*Augmat_));
      Problem.SetLHS(&soln);
      Problem.SetRHS(&rhs);

      // now we can allocate the AztecOO solver
      AztecOO kktsolver(Problem);

      // specify solver
      kktsolver.SetAztecOption(AZ_solver,AZ_gmres);
      kktsolver.SetAztecOption(AZ_kspace, 300 );
      kktsolver.SetAztecOption(AZ_output,AZ_last);

      // HERE WE SET THE IFPACK PRECONDITIONER
      kktsolver.SetPrecOperator(Prec);

      // .. and here we solve
      kktsolver.Iterate(300,1e-12);

      // delete the preconditioner
      delete Prec;
    }
    else if (solverChoice == 13) {
      Epetra_LinearProblem Problem;
      Problem.SetOperator(&(*Augmat_));
      Problem.SetLHS(&soln);
      Problem.SetRHS(&rhs);
      
      EpetraExt::LinearProblem_Reindex reindex(NULL);
      Epetra_LinearProblem newProblem = reindex(Problem);
      
      Amesos_Klu kktsolver(newProblem);
   
      AMESOS_CHK_ERR(kktsolver.SymbolicFactorization());
      AMESOS_CHK_ERR(kktsolver.NumericFactorization());
      AMESOS_CHK_ERR(kktsolver.Solve());
      kktsolver.PrintTiming();
    }
    
    
    for (int i=0; i<rhsy->MyLength(); i++) {
      (*((*y)(0)))[i] = soln[i];
    }
    for (int i=0; i<rhsu->MyLength(); i++) {
      (*((*u)(0)))[i] = soln[i+rhsy->MyLength()];
    }
    for (int i=0; i<rhsp->MyLength(); i++) {
      (*((*p)(0)))[i] = soln[i+rhsy->MyLength()+rhsu->MyLength()];
    }
    
  }
  else if (systemChoice == 2) {
    // We're using the Schur complement formulation to solve the augmented system.
  
    // Form linear operator.
    GLpApp::SchurOp schurop(A_, B_, Npy_);
  
    // Form Schur complement right-hand side.
    Epetra_MultiVector ny( (Epetra_BlockMap&)Npy_->RangeMap(), 1);
    Epetra_MultiVector ay( (Epetra_BlockMap&)A_->RangeMap(), 1);
    Epetra_MultiVector schurrhs( (Epetra_BlockMap&)A_->RangeMap(), 1);
    Epetra_MultiVector bu( (Epetra_BlockMap&)B_->RangeMap(), 1);
    A_->Multiply(false, *rhsy, ay);
    Npy_->Multiply(false, *rhsy, ny);
    B_->Multiply(false, *rhsu, bu);
    schurrhs.Update(1.0, ny, 1.0, ay, 0.0);
    schurrhs.Update(-1.0, *rhsp, 1.0, bu, 1.0);
  
    p->PutScalar(0.0);
    Epetra_LinearProblem linprob(&schurop, &(*p), &schurrhs);
    AztecOO::AztecOO Solver(linprob);
    Solver.SetAztecOption( AZ_solver, AZ_cg );
    Solver.SetAztecOption( AZ_precond, AZ_none );
    Solver.SetAztecOption( AZ_output, AZ_none );
    Solver.Iterate(8000,tol);
  
    Epetra_MultiVector bp( (Epetra_BlockMap&)B_->DomainMap(), 1);
    B_->Multiply(true, *p, bp);
    u->Update(1.0, *rhsu, -1.0, bp, 0.0);

    Epetra_MultiVector ap( (Epetra_BlockMap&)A_->DomainMap(), 1);
    Epetra_MultiVector np( (Epetra_BlockMap&)A_->DomainMap(), 1);
    A_->Multiply(true, *p, ap);
    Npy_->Multiply(true, *p, np);
    y->Update(1.0, *rhsy,0.0);
    y->Update(-1.0, ap, -1.0, np, 1.0);
  }
*/
  TEST_FOR_EXCEPT(true);
  return 0;
}


Teuchos::RefCountPtr<const Epetra_Comm> GLpYUEpetraDataPool::getCommPtr()   { return commptr_; }

Teuchos::RefCountPtr<Epetra_FECrsMatrix> GLpYUEpetraDataPool::getA()  { return A_; }

Teuchos::RefCountPtr<Epetra_FECrsMatrix> GLpYUEpetraDataPool::getB()  { return B_; }

Teuchos::RefCountPtr<Epetra_FECrsMatrix> GLpYUEpetraDataPool::getH()  { return H_; }

Teuchos::RefCountPtr<Epetra_FECrsMatrix> GLpYUEpetraDataPool::getR()  { return R_; }

Teuchos::RefCountPtr<Epetra_CrsMatrix> GLpYUEpetraDataPool::getAugmat()  { return Augmat_; }

Teuchos::RefCountPtr<Epetra_FECrsMatrix> GLpYUEpetraDataPool::getNpy()  { return Npy_; }

Teuchos::RefCountPtr<Epetra_FEVector> GLpYUEpetraDataPool::getb()  { return b_; }

Teuchos::RefCountPtr<Epetra_FEVector> GLpYUEpetraDataPool::getq()  { return q_; }

Teuchos::RefCountPtr<Epetra_FEVector> GLpYUEpetraDataPool::getNy()  { return Ny_; }

double GLpYUEpetraDataPool::getbeta()  { return beta_; }

Teuchos::RefCountPtr<const Epetra_SerialDenseMatrix> GLpYUEpetraDataPool::getipcoords()  { return ipcoords_; }

Teuchos::RefCountPtr<const Epetra_IntSerialDenseVector> GLpYUEpetraDataPool::getipindx()  { return ipindx_; }

Teuchos::RefCountPtr<const Epetra_SerialDenseMatrix> GLpYUEpetraDataPool::getpcoords()  { return pcoords_; }

Teuchos::RefCountPtr<const Epetra_IntSerialDenseVector> GLpYUEpetraDataPool::getpindx()  { return pindx_; }

Teuchos::RefCountPtr<const Epetra_IntSerialDenseMatrix> GLpYUEpetraDataPool::gett()  { return t_; }

Teuchos::RefCountPtr<const Epetra_IntSerialDenseMatrix> GLpYUEpetraDataPool::gete()  { return e_; }


void GLpYUEpetraDataPool::computeNy( const Teuchos::RefCountPtr<const Epetra_MultiVector> & y )
{
  Epetra_Map overlapmap(-1, pindx_->M(), (int*)(pindx_)->A(), 1, *commptr_);
  Epetra_Map standardmap(A_->DomainMap());
  Teuchos::RefCountPtr<Epetra_MultiVector> yoverlap = Teuchos::rcp(new Epetra_MultiVector(overlapmap, 1));
  Epetra_Import Importer(overlapmap, standardmap);
  yoverlap->Import(*y, Importer, Insert);
  nonlinvec(*commptr_, *ipindx_, *ipcoords_, *pindx_, *pcoords_, *t_, yoverlap, Ny_);
}


void GLpYUEpetraDataPool::computeNpy( const Teuchos::RefCountPtr<const Epetra_MultiVector> & y )
{
  Epetra_Map overlapmap(-1, pindx_->M(), (int*)(pindx_)->A(), 1, *commptr_);
  Epetra_Map standardmap(A_->DomainMap());
  Teuchos::RefCountPtr<Epetra_MultiVector> yoverlap = Teuchos::rcp(new Epetra_MultiVector(overlapmap, 1));
  Epetra_Import Importer(overlapmap, standardmap);
  yoverlap->Import(*y, Importer, Insert);
  nonlinjac(*commptr_, *ipindx_, *ipcoords_, *pindx_, *pcoords_, *t_, yoverlap, Npy_);
}


void GLpYUEpetraDataPool::computeAugmat()
{
  Epetra_Map standardmap(A_->DomainMap());
  Epetra_Map bdryctrlmap(B_->DomainMap());

  int indexBase = 1;

  int numstates = standardmap.NumGlobalElements();
  //int numcontrols = bdryctrlmap.NumGlobalElements();
  int nummystates = standardmap.NumMyElements();
  int nummycontrols = bdryctrlmap.NumMyElements();

  Epetra_IntSerialDenseVector KKTmapindx(2*nummystates+nummycontrols);
  
  
  // Build KKT map.
  Epetra_IntSerialDenseVector states(nummystates);
  Epetra_IntSerialDenseVector controls(nummycontrols);
  standardmap.MyGlobalElements(states.Values());
  bdryctrlmap.MyGlobalElements(controls.Values());
  for (int i=0; i<nummystates; i++) {
    KKTmapindx(i) = states(i);
    KKTmapindx(nummystates+nummycontrols+i) = 2*numstates+states(i);
  }
  for (int i=0; i<nummycontrols; i++) {
    KKTmapindx(nummystates+i) = numstates+controls(i);
  }
  Epetra_Map KKTmap(-1, 2*nummystates+nummycontrols, KKTmapindx.Values(), indexBase, *(commptr_));
  
  
  // Start building the KKT matrix.
  
  Augmat_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, KKTmap, 0));

  double one[1];
  one[0]=1.0;
  for (int i=0; i<nummystates+nummycontrols; i++) {
    Augmat_->InsertGlobalValues(KKTmapindx.Values()[i], 1, one, KKTmapindx.Values()+i);
  }
  
  int checkentries=0;
  int nummyentries=0;
  Epetra_SerialDenseVector values(nummyentries);
  Epetra_IntSerialDenseVector indices(nummyentries);
  // Insert A and Npy into Augmat.
  for (int i=0; i<nummystates; i++) {
    nummyentries = A_->NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    A_->ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                             indices.Values());
    if (nummyentries > 0)
      Augmat_->InsertGlobalValues(KKTmapindx.Values()[i]+2*numstates, nummyentries, values.Values(), 
                                  indices.Values());
    nummyentries = Npy_->NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    Npy_->ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                             indices.Values());
    if (nummyentries > 0)
      Augmat_->InsertGlobalValues(KKTmapindx.Values()[i]+2*numstates, nummyentries, values.Values(), 
                                  indices.Values());
  }
  // Insert B into Augmat.
  for (int i=0; i<nummystates; i++) {
    nummyentries = B_->NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    B_->ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                             indices.Values());
    for (int j=0; j<nummyentries; j++)
      indices[j] = indices[j]+numstates;
    if (nummyentries > 0)
      Augmat_->InsertGlobalValues(KKTmapindx.Values()[i]+2*numstates, nummyentries, values.Values(), 
                                  indices.Values());
  }
  
  bool MakeDataContiguous = false;
  EpetraExt::RowMatrix_Transpose transposer( MakeDataContiguous );
  Epetra_CrsMatrix & transA = dynamic_cast<Epetra_CrsMatrix&>(transposer(*A_));
  Epetra_CrsMatrix & transB = dynamic_cast<Epetra_CrsMatrix&>(transposer(*B_));
  Epetra_CrsMatrix & transNpy = dynamic_cast<Epetra_CrsMatrix&>(transposer(*Npy_));
  // Insert transpose of A and Npy into Augmat.
  for (int i=0; i<nummystates; i++) {
    nummyentries = transA.NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    transA.ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                                indices.Values());
    for (int j=0; j<nummyentries; j++)
      indices[j] = indices[j]+2*numstates;
    if (nummyentries > 0)
      Augmat_->InsertGlobalValues(KKTmapindx.Values()[i], nummyentries, values.Values(), 
                                  indices.Values());
    nummyentries = transNpy.NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    transNpy.ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                                  indices.Values());
    for (int j=0; j<nummyentries; j++)
      indices[j] = indices[j]+2*numstates;
    if (nummyentries > 0)
      Augmat_->InsertGlobalValues(KKTmapindx.Values()[i], nummyentries, values.Values(), 
                                  indices.Values());
  }
  // Insert transpose of B into Augmat.
  for (int i=0; i<nummystates; i++) {
    nummyentries = transB.NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    transB.ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                                indices.Values());
    for (int j=0; j<nummyentries; j++)
      indices[j] = indices[j]+2*numstates;
    int err = 0;
    if (nummyentries > 0)
      err = Augmat_->InsertGlobalValues(KKTmapindx.Values()[i]+numstates, nummyentries,
                                        values.Values(), indices.Values());
    // This will give a nasty message if something goes wrong with the insertion of B transpose.
    if (err < 0) {
      cout << "Insertion of entries failed:\n";
      cout << indices;
      cout << nummyentries << endl;
      cout << "at row: " << KKTmapindx.Values()[i]+numstates << endl << endl;
    }
  }

  Augmat_->FillComplete(KKTmap, KKTmap);
  // End building the KKT matrix.

}


void GLpYUEpetraDataPool::PrintVec( const Teuchos::RefCountPtr<const Epetra_Vector> & x )
{
  Vector2MATLAB(*x, cout);
}  

} //namespace GLpApp
