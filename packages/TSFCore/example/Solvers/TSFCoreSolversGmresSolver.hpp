// //////////////////////////////////////////////////////////////////////////
// TSFCoreSolversGmresSolver.hpp
//

#ifndef TSFCORE_SOLVERS_GMRES_SOLVER_HPP
#define TSFCORE_SOLVERS_GMRES_SOLVER_HPP

#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreTypes.hpp"
#include "Teuchos_DenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_BLAS.hpp"

namespace TSFCore {
namespace Solvers {

template<class Scalar>
class GMRESSolver {

public:
	GMRESSolver( 	const LinearOp<Scalar> 	&Op_In, 
			const Vector<Scalar> 	&b_In,
			Vector<Scalar> 		*curr_soln, 
			const ETransp		default_trans = NOTRANS,
			int			default_max_iter = 1000,
			Scalar			default_tol = 1e-10			
		  	);
				
	int currIteration() const { return curr_iter; };

	Scalar currEstRelResidualNorm() const { return (curr_res/r0); };

	void solve();	
	
private:

	void doIteration();

	int				max_iter, curr_iter;
	bool				isConverged;
	Scalar				tol, curr_res, actual_res, r0; 
	Scalar				*ptr_H, *z;
	const LinearOp<Scalar>		&Op;
	const Vector<Scalar>		&b;
	Vector<Scalar>			*x;
	ETransp				Op_trans;
	MemMngPack::ref_count_ptr<MultiVector<Scalar> >		
					V;
	MemMngPack::ref_count_ptr<Teuchos::DenseMatrix<int,Scalar> >
					H;
	MemMngPack::ref_count_ptr<Vector<Scalar> >
					w, r;
	std::vector<Scalar>		cs, sn;	

}; // end class GMRESSolver


template<class Scalar>
GMRESSolver<Scalar>::GMRESSolver(	const LinearOp<Scalar>	&Op_In,
					const Vector<Scalar>	&b_In,
					Vector<Scalar>		*curr_soln,
					const ETransp		default_trans,
					int			default_max_iter,
					Scalar			default_tol
				)
    : Op(Op_In), x(curr_soln), b(b_In), Op_trans(default_trans), max_iter(default_max_iter), tol(default_tol),
    isConverged(false), curr_iter(0), ptr_H(0), z(0)
{	
    // 
    //  Check compatability of linear operator with rhs and solution vector.
    //
    const VectorSpace<Scalar> &Op_domain = *Op.domain();
    const VectorSpace<Scalar> &Op_range = *Op.range();
    const bool 	domain_compatable = Op_domain.isCompatable( x->range() ),
		range_compatable = Op_range.isCompatable( b.range() );
    if (!(domain_compatable && range_compatable) ) { cout<<"Op is not compatable with x or b"<<endl; }
    //
    //  Initialize internal data structures
    //		
    V = Op_domain.createMembers(max_iter+1);
    w = Op_domain.createMember();
    r = Op_domain.createMember();
    H->shape(max_iter+1, _max_iter);
    ptr_H = H->values;
    z = new Scalar[max_iter+1];
    cs.resize(max_iter); sn.resize(max_iter); 
    //
    // Determine the residual from the current solution.
    //
    assign( r, b );
    Op.apply( Op_trans, x, r, 1.0, -1.0 );
    curr_res = norm_2( r ); 	        
    if (curr_res < tol) { isConverged = true; }
    //
    // Set up initial vector.
    //
    r0 = curr_res;
    z[0] = curr_res;
    w = V.col(1);		// get a mutable view of the first column of V.
    assign( w, r );		// copy r to the first column of V.
    Vt_S( w, 1.0/r0 );		// v_1 = r_0 / ||r_0||
}	
	
template<class Scalar>
void GMRESSolver<Scalar>::solve()
{
//
// Calls doIteration() using the current linear system parameters.
//
   while( !isConverged || curr_iter < max_iter ) { doIteration(); }
//
// Solve least squares problem.
//
   Teuchos::LAPACK<int,Scalar> lapack;
   lapack.TRSM("L", "U", "N", "N", curr_iter, 1, 1.0, ptr_H, 
		max_iter+1, z, max_iter+1 );
//
// Compute the new solution.
//
   for( int i = 1; i < curr_iter+1; i++ )
	Vp_StV( x, z[i-1], V.col(i) );
}

template<class Scalar>
void GMRESSolver<Scalar>::doIteration()
{
    int i;
    Scalar temp;
    Teuchos::BLAS<int, Scalar> blas;
    // 
    w = V.col(curr_iter+2);				// w = v_{j+1}
    Op.apply( Op_trans, V.col(curr_iter+1), w );	// w = Op * v_{j}	
    //
    // Perform MGS to orthogonalize new Krylov vector.
    //
    for( i=0; i<curr_iter+1; i++ ) {	
	H( i, curr_iter ) = dot( w, V.col(i+1) );	// h_{i,j} = ( w, v_{i} )
	Vt_StV( w, H(i, curr_iter), V.col(i+1) );	// w = w - h_{i,j} * v_{i}
    }
    H( curr_iter+1, curr_iter ) = norm_2( w );		// h_{j+1,j} = || w ||
    Vt_S( w, 1.0 / H( curr_iter+1, curr_iter ) ); 	// v_{j+1} = w / h_{j+1,j}			
    //
    // Apply previous Givens rotations
    //
    for( i=0; i<curr_iter; i++ ) {
	temp = cs[i]*H( i, curr_iter ) + sn[i]*H( i+1, curr_iter );
	H( i+1, curr_iter ) = -sn[i]*H( i, curr_iter ) + cs[i]*H( i+1, curr_iter );
	H( i, curr_iter ) = temp;
    }
    //
    // Calculate new Givens rotation
    //
    blas.ROTG( H( curr_iter, curr_iter ), H( curr_iter+1, curr_iter ), 
		cs[curr_iter], sn[curr_iter] );
    //
    // Update RHS and residual w/ new transform and compute residual.
    //
    z[curr_iter+1] = -sn[curr_iter]*z[curr_iter];
    z[curr_iter] *= cs[curr_iter];
    curr_res = abs( z[curr_iter+1] ) / r0; 
    if (curr_res < tol) { isConverged = true; }
    //    
    // Increment the iteration counter.
    //
    curr_iter++;
}

} // namespace Solvers
} // namespace TSFCore

#endif // TSFCORE_SOLVERS_GMRES_SOLVER_HPP
