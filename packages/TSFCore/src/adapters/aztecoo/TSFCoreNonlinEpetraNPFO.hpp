// //////////////////////////////////////////////////////
// TSFCoreNonlinEpetraNPFO.hpp

#ifndef TSFCORE_NONLIN_EPETRA_NPFO_HPP
#define TSFCORE_NONLIN_EPETRA_NPFO_HPP

#include "TSFCoreNonlinNonlinearProblemFirstOrder.hpp"
#include "TSFCoreNonlinLinearOpWithSolveAztecOO.hpp"
#include "Epetra_NonlinearProblemFirstOrder.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraVector.hpp"
#include "Ifpack_PrecGenerator.hpp"

namespace TSFCore {
namespace Nonlin {

///
/** Adds first-derivatives to <tt>NonlinearProblem</tt> as Epetra objects.
 *
 * Note: this class is setup for the default case where only a square
 * set if equations are defined and if optimization is performed then
 * DcDu are all <tt>Epetra_MultiVector</tt> objects.
 *
 * ToDo: Finish documentation!
 */
class EpetraNPFO : public NonlinearProblemFirstOrder<double> {
public:

  ///
  /** Give mutable access to the object used to generate
   * preconditioners.
   *
   * The purpose of this function is to allow clients to change the
   * options that affect how preconditioners are generated.
   */
  Ifpack::PrecGenerator& precGenerator();

  ///
  /** Same as above except is constant.
   */
  const Ifpack::PrecGenerator& precGenerator() const;

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

  // ///////////////////////////////////////
  // Private types

  ///
  class DcDu_Allocator {
  public:
    
    /** @name Constructors */
    //@{
    
    ///
    DcDu_Allocator(
      const bool  useEO
      );
      
    //@}
      
    /** @name AbstractFactoryStd compliant interface */
    //@{
    
    ///
    typedef Teuchos::RefCountPtr<TSFCore::LinearOp<double> >  ptr_t;
    ///
    const ptr_t allocate() const;

    //@}
    
  private:
    
    DcDu_Allocator();  // Not derfined and not to be called
    
    const bool  useEO_;
    
  }; // class DcDu_Allocator

	// //////////////////////////////////////
	// Private data members

	bool isInitialized_;

  Ifpack::PrecGenerator   precGenerator_;

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

	mutable std::vector<const Epetra_Vector*>  u_in_;

  mutable bool c_updated_, g_updated_, DcDy_updated_, DgDy_updated_;
  mutable std::vector<bool> DcDu_updated_, DgDu_updated_;

  EpetraVector                                    *c_;
	EpetraVector                                    *g_;
	LinearOpWithSolveAztecOO                        *DcDy_;
  std::vector<EpetraLinearOp*>                    DcDu_op_;
  std::vector<EpetraMultiVector*>                 DcDu_mv_;
  EpetraMultiVector                               *DgDy_;
  std::vector<EpetraMultiVector*>                 DgDu_;

  mutable std::vector<Teuchos::RefCountPtr<Epetra_Operator> >     epetra_DcDu_op_;
  mutable std::vector<Teuchos::RefCountPtr<Epetra_MultiVector> >  epetra_DcDu_mv_;
  mutable std::vector<Epetra::EpetraOp_or_EpetraMV>               epetra_DcDu_args_;

  mutable std::vector<Teuchos::RefCountPtr<Epetra_MultiVector> >  epetra_DgDu_;
  mutable std::vector<Epetra_MultiVector*>                        epetra_DgDu_args_;

	// //////////////////////////////////////
	// Private member functions

  ///
  static const Epetra_Vector& get_epetra_vec( const Vector<Scalar> &v );

  //
	const Epetra_Vector** set_u( const Vector<Scalar>* u[], bool newPoint ) const;

  //
	void updateNewPoint( bool newPoint ) const;
  
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

// public

inline
Ifpack::PrecGenerator& EpetraNPFO::precGenerator()
{
  return precGenerator_;
}

inline
const Ifpack::PrecGenerator& EpetraNPFO::precGenerator() const
{
  return precGenerator_;
}

// private

inline
const Epetra_Vector& EpetraNPFO::get_epetra_vec( const Vector<Scalar> &v )
{
  using DynamicCastHelperPack::dyn_cast;
  return *dyn_cast<const TSFCore::EpetraVector>(v).epetra_vec();
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_EPETRA_NPFO_HPP
