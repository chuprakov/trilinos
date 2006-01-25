#include "Epetra_LAPACK.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"

#include "includes.h"

#include "usr_par.h"

double determinant(const Epetra_SerialDenseMatrix &);

int inverse(const Epetra_SerialDenseMatrix &, Epetra_SerialDenseMatrix &);

/** \brief Computes local stiffness matrix and local RHS vector for simplex
           (triangles in two dimensions).
                            
  \param  vertices   [in]  - Matrix containing the global coordinates of the current simplex.
  \param  k          [in]  - Vector containing the value of the diffusion \f$k\f$ at each vertex
                             (\f$k\f$ is assumed to be piecewise linear), where \f$k\f$ is the coeff in the
                             term \f$ \nabla \cdot (k \nabla(u)) \f$ in the advection-diffusion equation.
  \param  c          [in]  - Matrix containing the value of the advection field \f$ \mathbf{c} \f$ at each
                             vertex (\f$ \mathbf{c} \f$ is assumed to be piecewise linear), where
                             \f$ \mathbf{c} \f$ is the 2-vector of coeffs in the term
                             \f$ \mathbf{c}(x) \cdot \nabla y(x) \f$ in the advection-diffusion equation.
  \param  r          [in]  - Vector containing the value of \f$ r \f$ at each vertex (\f$ r \f$ is assumed
                             to be piecewise linear), where \f$ r \f$ is the coefficient in the term
                             \f$ ru \f$ in the advection-diffusion equation.
  \param  f          [in]  - Vector containing the value of \f$ f \f$ at each vertex (\f$ f \f$ is assumed to be
                             piecewise linear), where \f$ f \f$ is the right hand side in the adv-diff eq
  \param  usr_par    [in]  - class containing: \n
                              - S1, S2, S3 (3x3) the integrals of various combinations of partials
                                of local basis functions
                              - N (1x3) integrals of local basis functions
                              - NNN[3] (3x3), integrals of products of local basis functions N(i)*N(j)*N(k)
                              - etc.
  \param  At         [out] - stiffness matrix contribution for the simplex
  \param  bt         [out] - right-hand side contribution for the simplex

  \return 0                 if successful.

  \par Detailed Description:

     Computes the local (per simplex) contributions to the FE volume stiffness matrix \e A and the
     right-hand side \e b for the advection-diffusion equation
    \f{align*}
       - \nabla \cdot \left( K(x) \nabla y(x) \right) + \mathbf{c}(x) \cdot \nabla y(x) + r(x) y(x)
       &= f(x), & x &\in \Omega,\;  \\
       y(x) &= d(x),    & x &\in {\partial \Omega}_d, \\
       K(x) \frac{\partial}{\partial \mathbf{n}}  y(x) &= g(x), & x &\in {\partial \Omega}_n, \\
       \sigma_0 K(x) \frac{\partial}{\partial \mathbf{n}} y(x)
       + \sigma_1 y(x) &= g(x), & x &\in {\partial \Omega}_r.
    \f}
*/
int lassembly(const Epetra_SerialDenseMatrix & vertices,
              const Epetra_SerialDenseVector & k,
              const Epetra_SerialDenseMatrix & c,
              const Epetra_SerialDenseVector & r,
              const Epetra_SerialDenseVector & f,
              const Usr_Par & usr_par,
              Epetra_SerialDenseMatrix & At,
              Epetra_SerialDenseVector & bt)
{
  // Note that the constructors below initialize all entries to 0.
  Epetra_SerialDenseMatrix B(2,2);
  Epetra_SerialDenseMatrix b(2,2);
  Epetra_SerialDenseMatrix BtB(2,2);  
  Epetra_SerialDenseMatrix C(2,2);  
  Epetra_SerialDenseMatrix M1(3,3);
  Epetra_SerialDenseMatrix M2(3,3);
  Epetra_SerialDenseMatrix M3(3,3);
  Epetra_SerialDenseMatrix tmp(3,3);
  double dB, adB;
  At.Shape(3,3);
  bt.Size(3);

  // Construct affine transformation matrix.
  for(int i=0; i<2; i++) {
    B(i,0) = vertices(1,i)-vertices(0,i);
    B(i,1) = vertices(2,i)-vertices(0,i);
  }
  dB  = determinant(B);
  adB = abs(dB);

  // Compute matrix C = inv(B'*B).
  BtB.Multiply('T', 'N', 1.0, B, B, 0.0);
  inverse(BtB, C);

  inverse(B, b); b.Scale(dB);
  
  // Compute the part corresponding to div(K*grad(u)).
  tmp = usr_par.S1; tmp.Scale(C(0,0));
  M1 += tmp;
  tmp = usr_par.S2; tmp.Scale(C(0,1));
  M1 += tmp;
  tmp = usr_par.S3; tmp.Scale(C(1,1));
  M1 += tmp;
  M1.Scale( (k(0)*usr_par.Nw(0) + k(1)*usr_par.Nw(1) +
             k(2)*usr_par.Nw(2)) * adB );

  // Compute the part corresponding to c'*grad(u).
  tmp = usr_par.NdNdx1Nw[0]; tmp.Scale(b(0,0)*c(0,0)+b(0,1)*c(0,1));
  M2 += tmp;
  tmp = usr_par.NdNdx1Nw[1]; tmp.Scale(b(0,0)*c(1,0)+b(0,1)*c(1,1));
  M2 += tmp;
  tmp = usr_par.NdNdx1Nw[2]; tmp.Scale(b(0,0)*c(2,0)+b(0,1)*c(2,1));
  M2 += tmp;
  tmp = usr_par.NdNdx2Nw[0]; tmp.Scale(b(1,0)*c(0,0)+b(1,1)*c(0,1));
  M2 += tmp;
  tmp = usr_par.NdNdx2Nw[1]; tmp.Scale(b(1,0)*c(1,0)+b(1,1)*c(1,1));
  M2 += tmp;
  tmp = usr_par.NdNdx2Nw[2]; tmp.Scale(b(1,0)*c(2,0)+b(1,1)*c(2,1));
  M2 += tmp;
  M2.Scale(adB/dB);

  // Compute the part corresponding to r*u.
  tmp = usr_par.NNNw[0]; tmp.Scale(r(0));
  M3 += tmp;
  tmp = usr_par.NNNw[1]; tmp.Scale(r(1));
  M3 += tmp;
  tmp = usr_par.NNNw[2]; tmp.Scale(r(2));
  M3 += tmp;
  M3.Scale(adB);

  At += M1;
  At += M2;
  At += M3;

  bt.Multiply('T', 'N', adB, usr_par.NNw, f, 0.0);  
  
  return(0);
}

/**  \brief Computes the inverse of a dense matrix.

  \param  mat  [in]  - the matrix that is to be inverted
  \param  inv  [in]  - the inverse

  \return 0            if successful
*/
int inverse(const Epetra_SerialDenseMatrix & mat,
            Epetra_SerialDenseMatrix & inv)
{
  Epetra_LAPACK lapack;
  int dim = mat.M();
  int info;
  Epetra_IntSerialDenseVector ipiv(dim);
  Epetra_SerialDenseVector work(2*dim);

  inv.Shape(dim,dim);
  inv = mat;

  lapack.GETRF(dim, dim, inv.A(), dim, ipiv.A(), &info);
  lapack.GETRI(dim, inv.A(), dim, ipiv.A(), work.A(), &dim, &info);
  
  return(0);
}


/**  \brief Computes the determinant of a dense matrix.

  \param  mat  [in]  - the matrix

  \return the determinant
*/
double determinant(const Epetra_SerialDenseMatrix & mat)
{
  //Teuchos::LAPACK<int,double> lapack;
  Epetra_LAPACK lapack;
  double det;
  int swaps; 
  int dim = mat.M();
  int info;
  Epetra_IntSerialDenseVector ipiv(dim);
 
  Epetra_SerialDenseMatrix mymat(mat);
  lapack.GETRF(dim, dim, mymat.A(), dim, ipiv.A(), &info);

  det = 1.0;
  for (int i=0; i<dim; i++) {
    det *= mymat(i,i);
  }

  swaps = 0;
  for (int i=0; i<dim; i++) {
    if ((ipiv[i]-1) != i)
      swaps++;
  }

  det *= pow((double)(-1.0),swaps);

  return(det);
}

