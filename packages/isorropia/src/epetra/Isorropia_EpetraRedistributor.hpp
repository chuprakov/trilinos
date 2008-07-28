//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc, 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Alan Williams (william@sandia.gov)
                or Erik Boman    (egboman@sandia.gov)

************************************************************************
*/
//@HEADER

#ifndef _Isorropia_EpetraRedistributor_hpp_
#define _Isorropia_EpetraRedistributor_hpp_

#include <Isorropia_Redistributor.hpp>
#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_EPETRA
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_LinearProblem;
class Epetra_SrcDistObject;
class Epetra_DistObject;

namespace Isorropia {

namespace Epetra {
  class Partitioner;

/** Class which is constructed with a Partitioner instance, and
     provides several methods for redistributing Epetra objects
     given the partitioning computed by the Partitioner object.
*/

class Redistributor : public Isorropia::Redistributor {
public:
  /** This constructor is part of the Isorropia API.

      This constructor calls partitioner.compute_partitioning() if it
      has not already been called.
   */
  Redistributor(Teuchos::RefCountPtr<Isorropia::Epetra::Partitioner> partitioner);

  /** Destructor
   */
  virtual ~Redistributor();

  /** This method, part of the Isorropia API, is used to redistribute an Epetra object.

      Method to redistribute a Epetra_SrcDistObject into a
      Epetra_DistObject. The caller is required to have constructed
      the target object using the correct target map.
  */
  void redistribute(const Epetra_SrcDistObject& src,
		    Epetra_DistObject& target);

  /** This method is part of the Isorropia API, used to redistribute an Epetra object.

      Method to accept a Epetra_CrsGraph object, and
      return a redistributed Epetra_CrsGraph object.

      Note that the 'input_graph' argument may be a
      different object than the one which was used to
      construct the partitioner.
  */
  Teuchos::RefCountPtr<Epetra_CrsGraph>
     redistribute(const Epetra_CrsGraph& input_graph, bool callFillComplete= true);

  /** This method is part of the Isorropia API, used
      to redistribute an Epetra object.

      Method to accept a Epetra_CrsMatrix object, and
      return a redistributed Epetra_CrsMatrix object.

      Note that the 'input_matrix' argument may be a
      different object than the one which was used to
      construct the partitioner.
  */
  Teuchos::RefCountPtr<Epetra_CrsMatrix>
     redistribute(const Epetra_CrsMatrix& input_matrix, bool callFillComplete= true);

  /** This method is part of the Isorropia API, used
      to redistribute an Epetra object.

      Method to accept a Epetra_RowMatrix object, and
      return a redistributed Epetra_CrsMatrix object.
  */
  Teuchos::RefCountPtr<Epetra_CrsMatrix>
     redistribute(const Epetra_RowMatrix& input_matrix, bool callFillComplete= true);

  /** This method is part of the Isorropia API, used
      to redistribute an Epetra object.

      Method to accept a Epetra_Vector object, and
      return a redistributed Epetra_Vector object.
  */
  Teuchos::RefCountPtr<Epetra_Vector>
     redistribute(const Epetra_Vector& input_vector);

  /** This method is part of the Isorropia API, used
      to redistribute an Epetra object.

      Method to accept a Epetra_MultiVector object, and
      return a redistributed Epetra_MultiVector object.
  */
  Teuchos::RefCountPtr<Epetra_MultiVector>  
     redistribute(const Epetra_MultiVector& input_vector);

  /** Reverse redistribute an Epetra_Vector.
  */
  void
     redistribute_reverse(const Epetra_Vector& input_vector, Epetra_Vector& output_vector);

  /** Reverse redistribute an Epetra_MultiVector.
   */
  void
     redistribute_reverse(const Epetra_MultiVector& input_vector, Epetra_MultiVector& output_vector);
private:
  /** This is an internal method of the Redistributor class.
   */
  void create_importer(const Epetra_BlockMap& src_map);

  Teuchos::RefCountPtr<Isorropia::Epetra::Partitioner> partitioner_;
  Teuchos::RefCountPtr<Epetra_Import> importer_;
  Teuchos::RefCountPtr<Epetra_Map> target_map_;

  bool created_importer_;

}; //class Redistributor

}//namespace Epetra

}//namespace Isorropia

#endif //HAVE_EPETRA

#endif

