// ////////////////////////////////////////////////////
// AbstractLinAlgPack_ThyraAccessors.hpp

#ifndef ALAP_THYRA_ACCESSORS_HPP
#define ALAP_THYRA_ACCESSORS_HPP

#include "AbstractLinAlgPack_VectorSpaceThyra.hpp"

namespace AbstractLinAlgPack {

///
void get_thyra_vector(
  const VectorSpaceThyra                                         &thyra_vec_spc
  ,const Vector                                                  &vec
  ,Teuchos::RefCountPtr<const Thyra::VectorBase<value_type> >    *thyra_vec
  );

///
void free_thyra_vector(
  const VectorSpaceThyra                                         &thyra_vec_spc
  ,const Vector                                                  &vec
  ,Teuchos::RefCountPtr<const Thyra::VectorBase<value_type> >    *thyra_vec
  );

///
void get_thyra_vector(
  const VectorSpaceThyra                                         &thyra_vec_spc
  ,VectorMutable                                                 *vec
  ,Teuchos::RefCountPtr<Thyra::VectorBase<value_type> >          *thyra_vec
  );

///
void commit_thyra_vector(
  const VectorSpaceThyra                                         &thyra_vec_spc
  ,VectorMutable                                                 *vec
  ,Teuchos::RefCountPtr<Thyra::VectorBase<value_type> >          *thyra_vec
  );

}

#endif // ALAP_THYRA_ACCESSORS_HPP
