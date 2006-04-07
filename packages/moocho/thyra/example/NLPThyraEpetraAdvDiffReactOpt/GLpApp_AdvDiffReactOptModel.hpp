#ifndef GLP_APP_ADV_DIFF_REACT_OPT_MODEL_HPP
#define GLP_APP_ADV_DIFF_REACT_OPT_MODEL_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "GLpApp_GLpYUEpetraDataPool.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"
#include "Teuchos_VerboseObject.hpp"

namespace GLpApp {

/** \brief 
 *
 * ToDo: Finish Documentation!
 */
class AdvDiffReactOptModel
  : public EpetraExt::ModelEvaluator
  , public Teuchos::VerboseObject<AdvDiffReactOptModel>
{
public:

  /** \brief Constructor. */
  AdvDiffReactOptModel(
    Teuchos::RefCountPtr<GLpApp::GLpYUEpetraDataPool>   const& dat
    ,const int                                                 np           = -1
    ,const double                                              x0           = 0.0
    ,const double                                              p0           = 1.0
    ,const double                                              reactionRate = 1.0
    );

  /** \brief . */
  void set_q( Teuchos::RefCountPtr<Epetra_Vector> const& q );

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Map> get_f_map() const;
  /** \breif . */
  Teuchos::RefCountPtr<const Epetra_Map> get_p_map(int l) const;
  /** \breif . */
  Teuchos::RefCountPtr<const Epetra_Map> get_g_map(int j) const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_p_init(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_x_lower_bounds() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_x_upper_bounds() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_p_lower_bounds(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_p_upper_bounds(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<Epetra_Operator> create_W() const;
  /** \brief . */
  Teuchos::RefCountPtr<Epetra_Operator> create_DfDp_op(int l) const;
  /** \brief . */
  InArgs createInArgs() const;
  /** \brief . */
  OutArgs createOutArgs() const;
  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  //@}

private:

  // /////////////////////////////////////
  // Private member data

	bool      isInitialized_;

  Teuchos::RefCountPtr<GLpApp::GLpYUEpetraDataPool>   dat_;
  int                                                 np_;
  double                                              reactionRate_;
	Teuchos::RefCountPtr<Epetra_Vector>                 q_;

	Teuchos::RefCountPtr<const Epetra_Map>              map_p_bar_;
	Teuchos::RefCountPtr<Epetra_MultiVector>            B_bar_;

  Teuchos::RefCountPtr<const Epetra_Comm>  epetra_comm_;
	Teuchos::RefCountPtr<const Epetra_Map>   map_x_;
	Teuchos::RefCountPtr<const Epetra_Map>   map_p_;
	Teuchos::RefCountPtr<const Epetra_Map>   map_f_;
	Teuchos::RefCountPtr<const Epetra_Map>   map_g_;

	Teuchos::RefCountPtr<Epetra_Vector> xL_;
	Teuchos::RefCountPtr<Epetra_Vector> xU_;
	Teuchos::RefCountPtr<Epetra_Vector> pL_;
	Teuchos::RefCountPtr<Epetra_Vector> pU_;
	Teuchos::RefCountPtr<Epetra_Vector> gL_;
	Teuchos::RefCountPtr<Epetra_Vector> gU_;
	Teuchos::RefCountPtr<Epetra_Vector> x0_;
	Teuchos::RefCountPtr<Epetra_Vector> p0_;

  Teuchos::RefCountPtr<Epetra_CrsGraph>  W_graph_;

};

} // namespace GLpApp

#endif // GLP_APP_ADV_DIFF_REACT_OPT_MODEL_HPP
