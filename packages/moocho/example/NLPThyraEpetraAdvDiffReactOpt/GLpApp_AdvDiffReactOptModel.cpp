#include "GLpApp_AdvDiffReactOptModel.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"

namespace {

inline double sqr( const double& s ) { return s*s; }

} // namespace

namespace GLpApp {

AdvDiffReactOptModel::AdvDiffReactOptModel(
  Teuchos::RefCountPtr<GLpApp::GLpYUEpetraDataPool>   const& dat
  ,const double                                              x0
  ,const double                                              p0
  ,const bool                                                dumpAll
)
  :dat_(dat)
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  //
  // Get the maps
  //
  map_x_ = Teuchos::rcp(new Epetra_Map(dat_->getA()->OperatorDomainMap()));
  map_p_ = Teuchos::rcp(new Epetra_Map(dat_->getB()->OperatorDomainMap()));
  map_f_ = Teuchos::rcp(new Epetra_Map(dat_->getA()->OperatorRangeMap()));
  map_g_ = Teuchos::rcp(new Epetra_Map(1,0,Epetra_SerialComm()));
  //
  // Create vectors
  //
  x0_ = Teuchos::rcp(new Epetra_Vector(*map_x_));
  //xL_ = Teuchos::rcp(new Epetra_Vector(*map_x_));
  //xU_ = Teuchos::rcp(new Epetra_Vector(*map_x_));
  p0_ = Teuchos::rcp(new Epetra_Vector(*map_p_));
  //pL_ = Teuchos::rcp(new Epetra_Vector(*map_p_));
  //pU_ = Teuchos::rcp(new Epetra_Vector(*map_p_));
  //
  // Initialize the vectors
  //
  x0_->PutScalar(x0);
  p0_->PutScalar(p0);
  //
  // Initialize the graph for W
  //
  dat_->computeNpy(x0_);
  if(dumpAll) { *out << "\nA =\n"; { Teuchos::OSTab tab(out); dat_->getA()->Print(*out); } }
  if(dumpAll) { *out << "\nNpy =\n"; {  Teuchos::OSTab tab(out); dat_->getNpy()->Print(*out); } }
  W_graph_ = Teuchos::rcp(new Epetra_CrsGraph(dat_->getA()->Graph())); // Assume A and Npy have same graph!
  //
  isInitialized_ = true;
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RefCountPtr<const Epetra_Map>
AdvDiffReactOptModel::get_x_map() const
{
  return map_x_;
}

Teuchos::RefCountPtr<const Epetra_Map>
AdvDiffReactOptModel::get_f_map() const
{
  return map_x_;
}

Teuchos::RefCountPtr<const Epetra_Map>
AdvDiffReactOptModel::get_p_map(int l) const
{
  TEST_FOR_EXCEPT(l!=1);
  return map_p_;
}

Teuchos::RefCountPtr<const Epetra_Map>
AdvDiffReactOptModel::get_g_map(int j) const
{
  TEST_FOR_EXCEPT(j!=1);
  return map_g_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_x_init() const
{
  return x0_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_p_init(int l) const
{
  TEST_FOR_EXCEPT(l!=1);
  return p0_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_x_lower_bounds() const
{
  return xL_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_x_upper_bounds() const
{
  return xU_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_p_lower_bounds(int l) const
{
  TEST_FOR_EXCEPT(l!=1);
  return pL_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_p_upper_bounds(int l) const
{
  TEST_FOR_EXCEPT(l!=1);
  return pU_;
}

Teuchos::RefCountPtr<Epetra_Operator>
AdvDiffReactOptModel::create_W() const
{
  return Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));
}

EpetraExt::ModelEvaluator::InArgs
AdvDiffReactOptModel::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  //inArgs.set_Np(1);
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.setSupports(IN_ARG_beta,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
AdvDiffReactOptModel::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  //outArgs.set_Np_Ng(1,1);
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_W_properties(
    DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_FULL
      ,true // supportsAdjoint
      )
    );
/*
  outArgs.setSupports(OUT_ARG_DfDp,1,DERIV_MV_BY_COL);
  outArgs.set_DfDp_properties(
    1,DerivativeProperties(
      DERIV_LINEARITY_CONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
  outArgs.setSupports(OUT_ARG_DgDx,1,DERIV_TRANS_MV_BY_ROW);
  outArgs.set_DgDx_properties(
    1,DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
  outArgs.setSupports(OUT_ARG_DgDp,1,1,DERIV_TRANS_MV_BY_ROW);
  outArgs.set_DgDp_properties(
    1,1,DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
*/
  return outArgs;
}

void AdvDiffReactOptModel::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  //
  // Get the input arguments
  //
  //Teuchos::RefCountPtr<const Epetra_Vector> p_in = inArgs.get_p(1);
  //const Epetra_Vector &p = (p_in.get() ? *p_in : *p0_);
  const Epetra_Vector &x = *inArgs.get_x();
  //
  // Get the output arguments
  //
  Teuchos::RefCountPtr<Epetra_Vector>       f_inout = outArgs.get_f();
  //Teuchos::RefCountPtr<Epetra_Vector>       g_inout = outArgs.get_g(1);
  Teuchos::RefCountPtr<Epetra_Operator>     W_inout = outArgs.get_W();
  //Teuchos::RefCountPtr<Epetra_MultiVector>  DfDp_inout = get_DfDp_mv(1,outArgs);
  //Teuchos::RefCountPtr<Epetra_MultiVector>  DgDx_trans_inout = get_DgDx_mv(1,outArgs,DERIV_TRANS_MV_BY_ROW);
  //Teuchos::RefCountPtr<Epetra_MultiVector>  DgDp_trans_inout = get_DgDp_mv(1,1,outArgs,DERIV_TRANS_MV_BY_ROW);
  //
  // Compute the functions
  //
  if(f_inout.get()) {
    Epetra_Vector &f = *f_inout;
    dat_->computeNy(Teuchos::rcp(&x,false));
    Epetra_Vector Ax(*map_f_);
    dat_->getA()->Multiply(false,x,Ax);
    Epetra_Vector Bp(*map_f_);
    dat_->getB()->Multiply(false,*p0_,Bp);
    f.Update(1.0,Ax,1.0,*dat_->getNy(),0.0);
    f.Update(1.0,Bp,-1.0,*dat_->getb(),1.0);
  }
/*
  if(g_inout.get()) {
    Epetra_Vector &g = *g_inout;
    TEST_FOR_EXCEPT(true);
  }
*/
  if(W_inout.get()) {
    const double beta = inArgs.get_beta();
    dat_->computeNpy(Teuchos::rcp(&x,false));
    Epetra_CrsMatrix &DfDx = dyn_cast<Epetra_CrsMatrix>(*W_inout);
    Teuchos::RefCountPtr<Epetra_CrsMatrix>
      dat_A = dat_->getA(),
      dat_Npy = dat_->getNpy();
    const int numMyRows = dat_A->NumMyRows();
    for( int i = 0; i < numMyRows; ++i ) {
      int dat_A_num_row_entries=0; double *dat_A_row_vals=0; int *dat_A_row_inds=0;
      dat_A->ExtractMyRowView(i,dat_A_num_row_entries,dat_A_row_vals,dat_A_row_inds);
      int dat_Npy_num_row_entries=0; double *dat_Npy_row_vals=0; int *dat_Npy_row_inds=0;
      dat_Npy->ExtractMyRowView(i,dat_Npy_num_row_entries,dat_Npy_row_vals,dat_Npy_row_inds);
#ifdef _DEBUG
      TEST_FOR_EXCEPT(dat_A_num_row_entries!=dat_Npy_num_row_entries);
#endif
      int DfDx_num_row_entries=0; double *DfDx_row_vals=0; int *DfDx_row_inds=0;
      DfDx.ExtractMyRowView(i,DfDx_num_row_entries,DfDx_row_vals,DfDx_row_inds);
#ifdef _DEBUG
      TEST_FOR_EXCEPT(DfDx_num_row_entries!=dat_Npy_num_row_entries);
#endif
      for(int k = 0; k < DfDx_num_row_entries; ++k) {
#ifdef _DEBUG
        TEST_FOR_EXCEPT(dat_A_row_inds[k]!=dat_Npy_row_inds[k]||dat_A_row_inds[k]!=DfDx_row_inds[k]);
#endif
        DfDx_row_vals[k] = beta * (dat_A_row_vals[k] + dat_Npy_row_vals[k]);
      }
    }
  }
/*
  if(DfDp_inout.get()) {
    Epetra_MultiVector &DfDp = *DfDp_inout;
    TEST_FOR_EXCEPT(true);
  }
  if(DgDx_trans_inout.get()) {
    Epetra_Vector &DgDx_trans = *(*DgDx_trans_inout)(0);
    TEST_FOR_EXCEPT(true);
  }
  if(DgDp_trans_inout.get()) {
    Epetra_Vector &DgDp_trans = *(*DgDp_trans_inout)(0);
    TEST_FOR_EXCEPT(true);
  }
*/
}

} // namespace GLpApp
