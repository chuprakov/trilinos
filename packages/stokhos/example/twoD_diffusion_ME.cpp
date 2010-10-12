// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "twoD_diffusion_ME.hpp"

#include "Stokhos_KL_ExponentialRandomField.hpp"
#include "EpetraExt_MatrixMatrix.h"
#include "Teuchos_TestForException.hpp"

namespace {

  // Functor representing a diffusion function given by a KL expansion
  struct Stokhos_KL_Diffusion_Func {
    double mean;
    mutable Teuchos::Array<double> point;
    Teuchos::RCP< Stokhos::KL::ExponentialRandomField<double> > rf;
    
    Stokhos_KL_Diffusion_Func(double xyLeft, double xyRight, 
			      double mean_, double std_dev, 
			      double L, int num_KL) : mean(mean_), point(2)
    {
      Teuchos::ParameterList rfParams;
      rfParams.set("Number of KL Terms", num_KL);
      rfParams.set("Mean", mean);
      rfParams.set("Standard Deviation", std_dev);
      int ndim = 2;
      Teuchos::Array<double> domain_upper(ndim), domain_lower(ndim), 
	correlation_length(ndim);
      for (int i=0; i<ndim; i++) {
	domain_upper[i] = xyRight;
	domain_lower[i] = xyLeft;
	correlation_length[i] = L;
      }
      rfParams.set("Domain Upper Bounds", domain_upper);
      rfParams.set("Domain Lower Bounds", domain_lower);
      rfParams.set("Correlation Lengths", correlation_length);

      rf = 
	Teuchos::rcp(new Stokhos::KL::ExponentialRandomField<double>(rfParams));
    }

    double operator() (double x, double y, int k) const {
      if (k == 0)
	return mean;
      point[0] = x;
      point[1] = y;
      return rf->evaluate_eigenfunction(point, k-1);
    }
  };

  // Functor representing a diffusion function given by a KL expansion
  // where the basis is normalized
  template <typename DiffusionFunc>
  struct Normalized_KL_Diffusion_Func {
    const DiffusionFunc& func;
    int d;
    Teuchos::Array<double> psi_0, psi_1;
    
    Normalized_KL_Diffusion_Func(
      const DiffusionFunc& func_,
      const Stokhos::OrthogPolyBasis<int,double>& basis) : 
      func(func_),
      d(basis.dimension()), 
      psi_0(basis.size()), 
      psi_1(basis.size())    
    {
      Teuchos::Array<double> zero(d), one(d);
      for(int k=0; k<d; k++) {
	zero[k] = 0.0;
	one[k] = 1.0;
      }
      basis.evaluateBases(zero, psi_0);
      basis.evaluateBases(one, psi_1);
    }

    double operator() (double x, double y, int k) const {
      if (k == 0) {
	double val = func(x, y, 0);
	for (int i=1; i<=d; i++)
	  val -= psi_0[i]/(psi_1[i]-psi_0[i])*func(x, y, i);
	val /= psi_0[0];
	return val;
      }
      else
	return 1.0/(psi_1[k]-psi_0[k])*func(x, y, k);
    }
  };

  // Functor representing a diffusion function given by a log-normal PC expansion
  template <typename DiffusionFunc>
  struct LogNormal_Diffusion_Func {
    double mean;
    const DiffusionFunc& func;
    Teuchos::RCP<const Stokhos::ProductBasis<int, double> > prodbasis;

    LogNormal_Diffusion_Func(
      double mean_,
      const DiffusionFunc& func_,
      const Teuchos::RCP<const Stokhos::ProductBasis<int, double> > prodbasis_) 
      : mean(mean_), func(func_), prodbasis(prodbasis_) {}
    
    double operator() (double x, double y, int k) const {
      int d = prodbasis->dimension();
      const Teuchos::Array<double>& norms = prodbasis->norm_squared();
      Teuchos::Array<int> multiIndex = prodbasis->getTerm(k);
      double sum_g = 0.0, efval;
      for (int l=0; l<d; l++) {
	sum_g += std::pow(func(x,y,l+1),2);
      }
      efval = std::exp(mean + 0.5*sum_g)/norms[k];
      for (int l=0; l<d; l++) {
	efval *= std::pow(func(x,y,l+1),multiIndex[l]); 
      }
      return efval;
    }
  };

}

twoD_diffusion_ME::
twoD_diffusion_ME(
  const Teuchos::RCP<Epetra_Comm>& comm, int n, int d, 
  double s, double mu, 
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis_,
  bool log_normal_) :
  mesh(n),
  basis(basis_),
  log_normal(log_normal_)
{
  //////////////////////////////////////////////////////////////////////////////
  // Construct the mesh.  The mesh is just the tensor of the below array with 
  // itself.
  // The mesh is uniform and the nodes are numbered
  // LEFT to RIGHT, DOWN to UP.
  //
  // 5-6-7-8-9
  // | | | | |
  // 0-1-2-3-4
  /////////////////////////////////////////////////////////////////////////////
  double xyLeft = -.5;
  double xyRight = .5;
  h = (xyRight - xyLeft)/((double)(n-1));
  for(int idx = 0; idx < n; idx++){
    mesh[idx] = xyLeft + (idx)*h;
  }
  h = mesh[1]-mesh[0];
  
  // Solution vector map
  x_map = Teuchos::rcp(new Epetra_Map(n*n, 0, *comm));

  // Overlapped solution vector map
  x_overlapped_map = Teuchos::rcp(new Epetra_Map(n*n, 0, *comm));

  // Importer
  importer = Teuchos::rcp(new Epetra_Import(*x_overlapped_map, *x_map));

  // Initial guess, initialized to 0.0
  x_init = Teuchos::rcp(new Epetra_Vector(*x_map));
  x_init->PutScalar(0.0);

  // Overlapped solution vector
  x_overlapped = Teuchos::rcp(new Epetra_Vector(*x_overlapped_map));

  // Parameter vector map
  p_map = Teuchos::rcp(new Epetra_LocalMap(d, 0, *comm));

  // Response vector map
  g_map = Teuchos::rcp(new Epetra_LocalMap(1, 0, *comm));

  // Initial parameters
  p_init = Teuchos::rcp(new Epetra_Vector(*p_map));
  p_init->PutScalar(0.0);

  // Parameter names
  p_names = Teuchos::rcp(new Teuchos::Array<std::string>(d));
  for (int i=0;i<d;i++) {
    std::stringstream ss;
    ss << "KL Random Variable " << i+1;
    (*p_names)[i] = ss.str(); 
  }

  // Build Jacobian graph
  int NumMyElements = x_map->NumMyElements();
  int *MyGlobalElements = x_map->MyGlobalElements();
  int *NumNz = new int[NumMyElements];
  int Indices[4];
  int NumEntries;
  bcIndices.resize(NumMyElements);
  for (int i = 0; i<NumMyElements; i++) {
    // MyGlobalElements[i]<n ==> Boundary node on bottom edge.
    // MyGlobalElements[i]%n == 0 ==> Boundary node on left edge.
    // MyGlobalElements[i]+1%n == 0 ==> right edge.
    // MyGlobalElements[i] >= n - n ==> top edge.
    if (MyGlobalElements[i] < n || MyGlobalElements[i]%n == 0 ||
	(MyGlobalElements[i]+1)%n == 0 || MyGlobalElements[i] >= n*n - n) {
      NumNz[i] = 1;
      bcIndices[i] = 1;
    }
    else {
      NumNz[i] = 5;
      bcIndices[i] = 0;
    }
  }
  graph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *x_map, NumNz));
  for(int i=0; i<NumMyElements; ++i ) {
    if (bcIndices[i] == 0) {
      Indices[0] = MyGlobalElements[i]-n; //Down
      Indices[1] = MyGlobalElements[i]-1; //left
      Indices[2] = MyGlobalElements[i]+1; //right
      Indices[3] = MyGlobalElements[i]+n; //up
      NumEntries = 4;
    }
    if(bcIndices[i] == 0) 
      graph->InsertGlobalIndices(MyGlobalElements[i], NumEntries, Indices);
    graph->InsertGlobalIndices(MyGlobalElements[i], 1, MyGlobalElements+i);
  }
  graph->FillComplete();
  graph->OptimizeStorage();

  //KL_Diffusion_Func klFunc;
  Stokhos_KL_Diffusion_Func klFunc(xyLeft, xyRight, mu, s, 1.0, d);
  fillMatrices(klFunc, d+1);
  if (!log_normal) {
    // Fill coefficients of KL expansion of operator
    if (basis == Teuchos::null) { 
      fillMatrices(klFunc, d+1);
    }
    else {
      Normalized_KL_Diffusion_Func<Stokhos_KL_Diffusion_Func> nklFunc(klFunc, *basis);
      fillMatrices(nklFunc, d+1);
    }
  }
  else {
    // Fill coefficients of PC expansion of operator
    int sz = basis->size();
    Teuchos::RCP<const Stokhos::ProductBasis<int, double> > prodbasis =
      Teuchos::rcp_dynamic_cast<const Stokhos::ProductBasis<int, double> >(basis, true);
    LogNormal_Diffusion_Func<Stokhos_KL_Diffusion_Func> lnFunc(mu, klFunc, prodbasis);
    fillMatrices(lnFunc, sz);
  }

  // Construct deterministic operator
  A = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));
 
  // Construct the RHS vector.
  b = Teuchos::rcp(new Epetra_Vector(*x_map));
  for( int i=0 ; i<NumMyElements; ++i ) {
    if (bcIndices[i] == 1 )
      (*b)[i] = 0;
    else 
      (*b)[i] = 1;
  }

  delete [] NumNz;
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
twoD_diffusion_ME::
get_x_map() const
{
  return x_map;
}

Teuchos::RCP<const Epetra_Map>
twoD_diffusion_ME::
get_f_map() const
{
  return x_map;
}

Teuchos::RCP<const Epetra_Map>
twoD_diffusion_ME::
get_p_map(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  twoD_diffusion_ME::get_p_map():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return p_map;
}

Teuchos::RCP<const Epetra_Map>
twoD_diffusion_ME::
get_p_sg_map(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  twoD_diffusion_ME::get_p_sg_map():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return p_map;
}

Teuchos::RCP<const Epetra_Map>
twoD_diffusion_ME::
get_g_map(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  twoD_diffusion_ME::get_g_map():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return g_map;
}

Teuchos::RCP<const Epetra_Map>
twoD_diffusion_ME::
get_g_sg_map(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  twoD_diffusion_ME::get_g_sg_map():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return g_map;
}

Teuchos::RCP<const Teuchos::Array<std::string> >
twoD_diffusion_ME::
get_p_names(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  twoD_diffusion_ME::get_p_names():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return p_names;
}

Teuchos::RCP<const Teuchos::Array<std::string> >
twoD_diffusion_ME::
get_p_sg_names(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  twoD_diffusion_ME::get_p_sg_names():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return p_names;
}

Teuchos::RCP<const Epetra_Vector>
twoD_diffusion_ME::
get_x_init() const
{
  return x_init;
}

Teuchos::RCP<const Epetra_Vector>
twoD_diffusion_ME::
get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  twoD_diffusion_ME::get_p_init():  " <<
                     "Invalid parameter index l = " << l << std::endl);
  
  return p_init;
}

Teuchos::RCP<Epetra_Operator>
twoD_diffusion_ME::
create_W() const
{
  Teuchos::RCP<Epetra_CrsMatrix> AA = 
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));
  AA->FillComplete();
  AA->OptimizeStorage();
  return AA;
}

EpetraExt::ModelEvaluator::InArgs
twoD_diffusion_ME::
createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription("TwoD Diffusion Model Evaluator");

  // Deterministic InArgs
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.set_Np(1);    // 1 parameter vector

  // Stochastic InArgs
  inArgs.setSupports(IN_ARG_x_sg,true);
  inArgs.set_Np_sg(1); // 1 SG parameter vector
  inArgs.setSupports(IN_ARG_sg_basis,true);
  inArgs.setSupports(IN_ARG_sg_quadrature,true);
  inArgs.setSupports(IN_ARG_sg_expansion,true);
  
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
twoD_diffusion_ME::
createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription("TwoD Diffusion Model Evaluator");

  // Deterministic OutArgs
  outArgs.set_Np_Ng(1, 1);
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  
  // Stochastic OutArgs
  outArgs.set_Np_Ng_sg(1, 1);
  outArgs.setSupports(OUT_ARG_f_sg,true);
  outArgs.setSupports(OUT_ARG_W_sg,true);

  return outArgs;
}

void 
twoD_diffusion_ME::
evalModel(const InArgs& inArgs, const OutArgs& outArgs) const
{

  //
  // Determinisic calculation
  //

  // Solution vector
  Teuchos::RCP<const Epetra_Vector> det_x = inArgs.get_x();

  // Parameters
  Teuchos::RCP<const Epetra_Vector> p = inArgs.get_p(0);
  if (p == Teuchos::null)
    p = p_init;

  Teuchos::RCP<Epetra_Vector> f = outArgs.get_f();
  Teuchos::RCP<Epetra_Operator> W = outArgs.get_W();
  if (f != Teuchos::null || W != Teuchos::null) {
    *A = *(A_k[0]);
    for (int k=1;k<A_k.size();k++) {
      EpetraExt::MatrixMatrix::Add((*A_k[k]), false, (*p)[k-1], *A, 1.0);
    }
  }

  // Residual  
  if (f != Teuchos::null) {
    Teuchos::RCP<Epetra_Vector> kx = Teuchos::rcp(new Epetra_Vector(*x_map));
    A->Apply(*det_x,*kx);
    f->Update(1.0,*kx,-1.0, *b, 0.0);
  }

  // Jacobian
  if (W != Teuchos::null) {
    Teuchos::RCP<Epetra_CrsMatrix> jac =
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W, true);
    *jac = *A;
  }

  // Responses (mean value)
  Teuchos::RCP<Epetra_Vector> g = outArgs.get_g(0);
  if (g != Teuchos::null)
    det_x->MeanValue(&(*g)[0]);

  //
  // Stochastic Galerkin calculation
  //

  // Stochastic solution vector
  InArgs::sg_const_vector_t x_sg = inArgs.get_x_sg();

  // Stochastic parameters
  InArgs::sg_const_vector_t p_sg = inArgs.get_p_sg(0);

  // Stochastic residual
  OutArgs::sg_vector_t f_sg = outArgs.get_f_sg();
  if (f_sg != Teuchos::null) {

    // Get stochastic expansion data
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expn = 
      inArgs.get_sg_expansion();
    typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;
    Teuchos::RCP<const Cijk_type> Cijk = expn->getTripleProduct(); 
    const Teuchos::Array<double>& norms = basis->norm_squared();

    if (sg_kx_vec_all.size() != basis->size()) {
      sg_kx_vec_all.resize(basis->size());
      for (int i=0;i<basis->size();i++) {
	sg_kx_vec_all[i] = Teuchos::rcp(new Epetra_Vector(*x_map));
      }
    }
    f_sg->init(0.0);

    double pc_size;
    if(!log_normal) 
      pc_size = basis->dimension()+1;
    else
      pc_size = basis->size(); 
    for (int k=0; k<pc_size; k++) {
      for (Cijk_type::kj_iterator j_it = Cijk->j_begin(k); 
	   j_it != Cijk->j_end(k); ++j_it) {
	int j = Stokhos::index(j_it);
	A_k[k]->Apply((*x_sg)[j],*(sg_kx_vec_all[j]));
      }
      for (Cijk_type::kj_iterator j_it = Cijk->j_begin(k); 
	   j_it != Cijk->j_end(k); ++j_it) {
	int j = Stokhos::index(j_it);
	for (Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
	     i_it != Cijk->i_end(j_it); ++i_it) {
	  int i = Stokhos::index(i_it);
	  double c = Stokhos::value(i_it);  // C(i,j,k)
	  (*f_sg)[i].Update(1.0*c/norms[i],*(sg_kx_vec_all[j]),1.0);
	}
      }
    } //End 
    (*f_sg)[0].Update(-1.0,*b,1.0);
  }

  // Stochastic Jacobian
  OutArgs::sg_operator_t W_sg = outArgs.get_W_sg();
  if (W_sg != Teuchos::null) {
    for (int i=0; i<W_sg->size(); i++) {
      Teuchos::RCP<Epetra_CrsMatrix> jac = 
	Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_sg->getCoeffPtr(i), true);
      *jac = *A_k[i];
    }
  }

  // Stochastic responses
  Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly > g_sg = 
    outArgs.get_g_sg(0);
  if (g_sg != Teuchos::null) {
    int sz = x_sg->size();
    for (int i=0; i<sz; i++)
      (*x_sg)[i].MeanValue(&(*g_sg)[i][0]);
  }

}

