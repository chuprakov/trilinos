#include "usr_par.h"
#include "includes.h"


class Epetra_BLAS;
int compproduct(Epetra_SerialDenseVector &, double *, double *);
int compproduct(Epetra_SerialDenseVector &, double *, double *, double *);


Usr_Par::Usr_Par() {
  Epetra_BLAS blas;
  Epetra_SerialDenseVector product(4);

  // get nodes and weights
  quadrature(2,3,Nodes,Weights);
    
  // Evaluate nodal basis functions and their derivatives at quadrature
  // pts N(i,j) = value of the j-th basis function at quadrature node i.
  N.Shape(Nodes.M(),3);
  for (int i=0; i<Nodes.M(); i++) {
    N(i,0) = 1.0 - Nodes(i,0) - Nodes(i,1);
    N(i,1) = Nodes(i,0);
    N(i,2) = Nodes(i,1);
  }
  // Nx1(i,j) partial derrivative of basis function j wrt x1 at node i
  Nx1.Shape(Nodes.M(),3);
  for (int i=0; i<Nodes.M(); i++) {
    Nx1(i,0) = -1.0;
    Nx1(i,1) = 1.0;
    Nx1(i,2) = 0.0;
  }
  // Nx2(i,j) partial derrivative of basis function j wrt x2 at node i
  Nx2.Shape(Nodes.M(),3);
  for (int i=0; i<Nodes.M(); i++) {
    Nx2(i,0) = -1.0;
    Nx2(i,1) = 0.0;
    Nx2(i,2) = 1.0;
  }

  // Evaluate integrals of various combinations of partial derivatives
  // of the local basis functions (they're constant).
  S1.Shape(3,3);
  S1(0,0)= 1.0; S1(0,1)=-1.0; S1(0,2)= 0.0;
  S1(1,0)=-1.0; S1(1,1)= 1.0; S1(1,2)= 0.0;
  S1(2,0)= 0.0; S1(2,1)= 0.0; S1(2,2)= 0.0;
  S2.Shape(3,3);
  S2(0,0)= 2.0; S2(0,1)=-1.0; S2(0,2)=-1.0;
  S2(1,0)=-1.0; S2(1,1)= 0.0; S2(1,2)= 1.0;
  S2(2,0)=-1.0; S2(2,1)= 1.0; S2(2,2)= 0.0;
  S3.Shape(3,3);
  S3(0,0)= 1.0; S3(0,1)= 0.0; S3(0,2)=-1.0;
  S3(1,0)= 0.0; S3(1,1)= 0.0; S3(1,2)= 0.0;
  S3(2,0)=-1.0; S3(2,1)= 0.0; S3(2,2)= 1.0;
    
  // Evaluate integrals of basis functions N(i).
  Nw.Size(3);
  Nw.Multiply('T', 'N', 1.0, N, Weights, 0.0);

  // Evaluate integrals of 9 products N(i)*N(j).
  NNw.Shape(3,3);
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      compproduct(product, N[i], N[j]);
      NNw(i,j) = blas.DOT(Weights.M(), Weights.A(), product.A());
    }
  }

  // Evaluate integrals of 27 products N(i)*N(j)*N(k).
  NNNw = new Epetra_SerialDenseMatrix[3];
  NNNw[0].Shape(3,3); NNNw[1].Shape(3,3); NNNw[2].Shape(3,3); 
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        compproduct(product, N[i], N[j], N[k]);
        NNNw[k](i,j) = blas.DOT(Weights.M(), Weights.A(), product.A());
      }
    }
  }

  // Evaluate integrals of 27 products N(i)*(dN(j)/dx1)*N(k).
  NdNdx1Nw = new Epetra_SerialDenseMatrix[3];
  NdNdx1Nw[0].Shape(3,3); NdNdx1Nw[1].Shape(3,3); NdNdx1Nw[2].Shape(3,3); 
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        compproduct(product, N[i], Nx1[j], N[k]);
        NdNdx1Nw[k](i,j) = blas.DOT(Weights.M(), Weights.A(), product.A());
      }
    }
  }

  // Evaluate integrals of 27 products N(i)*(dN(j)/dx2)*N(k).
  NdNdx2Nw = new Epetra_SerialDenseMatrix[3];
  NdNdx2Nw[0].Shape(3,3); NdNdx2Nw[1].Shape(3,3); NdNdx2Nw[2].Shape(3,3); 
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        compproduct(product, N[i], Nx2[j], N[k]);
        NdNdx2Nw[k](i,j) = blas.DOT(Weights.M(), Weights.A(), product.A());
      }
    }
  }

}

void Usr_Par::Print(ostream& os) const {
  os << endl;
  os << "\n\nQuadrature nodes:";
  os << "\n-----------------";
  Nodes.Print(os);
  os << "\n\nQuadrature weights:";
  os << "\n-------------------\n";
  Weights.Print(os);
  os << "\n\nNodal basis functions:";
  os << "\n----------------------";
  N.Print(os);
  os << "\n\nIntegrals of combinations of partial derivatives:";
  os << "\n-------------------------------------------------";
  S1.Print(os);
  S2.Print(os);
  S3.Print(os);
  os << "\n\nIntegrals of basis functions:";
  os << "\n-----------------------------\n";
  Nw.Print(os);
  os << "\n\nIntegrals of products N(i)*N(j):";
  os << "\n--------------------------------\n";
  NNw.Print(os);
  os << "\n\nIntegrals of products N(i)*N(j)*N(k):";
  os << "\n-------------------------------------\n";
  NNNw[0].Print(os); NNNw[1].Print(os); NNNw[2].Print(os);
  os << "\n\nIntegrals of products N(i)*(dN(j)/dx1)*N(k):";
  os << "\n--------------------------------------------\n";
  NdNdx1Nw[0].Print(os); NdNdx1Nw[1].Print(os); NdNdx1Nw[2].Print(os);
  os << "\n\nIntegrals of products N(i)*(dN(j)/dx2)*N(k):";
  os << "\n--------------------------------------------\n";
  NdNdx2Nw[0].Print(os); NdNdx2Nw[1].Print(os); NdNdx2Nw[2].Print(os);
}


ostream& operator <<(ostream & out, const Usr_Par & usr_par)
{
  usr_par.Print(out);
  return out;
}


int compproduct(Epetra_SerialDenseVector & product,
                double *first, double *second)
{
  for (int i=0; i<product.M(); i++) {
    product[i] = first[i]*second[i];
  }
  return(0);
}


int compproduct(Epetra_SerialDenseVector & product,
                double *first, double *second, double *third)
{
  for (int i=0; i<product.M(); i++) {
    product[i] = first[i]*second[i]*third[i];
  }
  return(0);
}

