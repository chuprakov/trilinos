// //////////////////////////////////////////////////////
// TSFCoreNonlinEpetraNPFO.hpp

#ifndef TSFCORE_NONLIN_EPETRA_NPFO_HPP
#define TSFCORE_NONLIN_EPETRA_NPFO_HPP

#include "TSFCoreNonlinNonlinearProblemFirstOrder.hpp"
#include "Epetra_NonlinearProblemFirstOrder.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraVector.hpp"

namespace TSFCore {
namespace Nonlin {

class EpetraNPFO : public NonlinearProblemFirstOrder<double> {
public:

	/** @name Constructors / Initializers / accessors */
	//@{

  /// Construct to uninitialized
  EpetraNPFO();

  /// Calls <tt>initialize()</tt>
  EpetraNPFO(
    const Teuchos::RefCountPtr<Epetra::NonlinearProblemFirstOrder>   &epetra_np
    );

  ///
  void initialize(
    const Teuchos::RefCountPtr<Epetra::NonlinearProblemFirstOrder>   &epetra_np
    );

	//@}

	/** @name Overridden from NonlinearProblem */
	//@{

	///
	void initialize( bool testSetup );
	///
	bool isInitialized() const;
	///
	int Nu() const;
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> >  space_y() const;
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> >  space_u(int l) const;
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> >  space_c() const;
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> >  space_g() const;
	///
	const Vector<Scalar>& yL() const;
	///
	const Vector<Scalar>& yU() const;
	///
	const Vector<Scalar>& uL(int l) const;
	///
	const Vector<Scalar>& uU(int l) const;
	///
	const Vector<Scalar>& gL() const;
	///
	const Vector<Scalar>& gU() const;
	///
	const Vector<Scalar>& y0() const;
	///
	const Vector<Scalar>& u0(int l) const;
	///
	void set_c(Vector<Scalar>* c);
	///
	Vector<Scalar>* get_c();
	///
	void set_g(Vector<Scalar>* g);
	///
	Vector<Scalar>* get_g();
	///
	void unsetQuantities();
	///
	void calc_c(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
		) const;
	///
	void calc_g(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
		) const;

	//@}

	/** @name Overridden from NonlinearProblemFirstOrder */
	//@{

	///
	Teuchos::RefCountPtr< const MemMngPack::AbstractFactory<LinearOpWithSolve<Scalar> > > factory_DcDy() const;
	///
	Teuchos::RefCountPtr< const MemMngPack::AbstractFactory<LinearOp<Scalar > > > factory_DcDu(int l) const;
	///
	ETransp opDcDy() const;
	///
	ETransp opDcDu(int l) const;
	///
	void set_DcDy(LinearOpWithSolve<Scalar>* DcDy);
	///
	LinearOpWithSolve<Scalar>* get_DcDy();
	///
	void set_DcDu(int l, LinearOp<Scalar>* DcDu_l);
	///
	LinearOp<Scalar>* get_DcDu(int l);
	///
	void set_DgDy(MultiVector<Scalar>* DgDy);
	///
	MultiVector<Scalar>* get_DgDy();
	///
	void set_DgDu(int l, MultiVector<Scalar>* DgDu_l);
	///
	MultiVector<Scalar>* get_DgDu(int l);
	///
	void calc_DcDy(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
		) const;
	///
	void calc_DcDu(
		int                      l
		,const Vector<Scalar>    &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
		) const;
	///
	void calc_DgDy(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
		) const;
	///
	void calc_DgDu(
		int                      l
		,const Vector<Scalar>    &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
		) const;

	//@}

private:

	// //////////////////////////////////////
	// Private data members

	bool isInitialized_;

  Teuchos::RefCountPtr<Epetra::NonlinearProblemFirstOrder>        epetra_np_;

	Teuchos::RefCountPtr<const EpetraVectorSpace >                  space_y_;
	std::vector<Teuchos::RefCountPtr<const EpetraVectorSpace > >    space_u_;
	Teuchos::RefCountPtr<const EpetraVectorSpace >                  space_c_;
	Teuchos::RefCountPtr<const EpetraVectorSpace >                  space_g_;

	Teuchos::RefCountPtr<Vector<Scalar> >                           yL_;
	Teuchos::RefCountPtr<Vector<Scalar> >                           yU_;
	Teuchos::RefCountPtr<Vector<Scalar> >                           y0_;
	std::vector<Teuchos::RefCountPtr<Vector<Scalar> > >             uL_;
 	std::vector<Teuchos::RefCountPtr<Vector<Scalar> > >             uU_;
	std::vector<Teuchos::RefCountPtr<Vector<Scalar> > >             u0_;
	Teuchos::RefCountPtr<Vector<Scalar> >                           gL_;
 	Teuchos::RefCountPtr<Vector<Scalar> >                           gU_;

	Teuchos::RefCountPtr<const MemMngPack::AbstractFactory<LinearOpWithSolve<Scalar> > >       factory_DcDy_;
  std::vector<Teuchos::RefCountPtr<const MemMngPack::AbstractFactory<LinearOp<Scalar> > > >  factory_DcDu_;

  mutable bool c_updated_, g_updated_, DcDy_updated_, DgDy_updated_;
  mutable std::vector<bool> DcDu_updated_, DgDu_updated_;

  EpetraVector                                    *c_;
	EpetraVector                                    *g_;
//	LinearOpWithSolveIter<Scalar>                   *DcDy_;
  std::vector<EpetraLinearOp*>                    DcDu_;
  EpetraMultiVector                               *DgDy_;
  std::vector<EpetraMultiVector*>                 DgDu_;

	// //////////////////////////////////////
	// Private member functions

  ///
  static const Epetra_Vector& get_epetra_vec( const Vector<Scalar> &v );

  //
	void set_u( const Vector<Scalar>* u[], bool newPoint ) const;
  
  //
  void calc_Dc(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
    ,bool                    computeGradients
    ) const;

  //
  void calc_Dg(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]
		,bool                    newPoint
    ,bool                    computeGradients
    ) const;

	// Not defined and not to be called
	EpetraNPFO(const EpetraNPFO&);
	EpetraNPFO& operator=(const EpetraNPFO&);

}; // class EpetraNPFO

// ///////////////////////////////////
// Inline members

inline
const Epetra_Vector& EpetraNPFO::get_epetra_vec( const Vector<Scalar> &v )
{
  using DynamicCastHelperPack::dyn_cast;
  return *dyn_cast<const TSFCore::EpetraVector>(v).epetra_vec();
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_EPETRA_NPFO_HPP
