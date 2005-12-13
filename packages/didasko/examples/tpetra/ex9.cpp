//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_Version.hpp"
#ifdef TPETRA_MPI
#include "Tpetra_MpiPlatform.hpp"
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"
#endif
#include "Teuchos_ScalarTraits.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_Import.hpp"
#include <map>

// \author Marzio Sala, ETHZ/COLAB
//
// \date Last updated on 12-Dec-05

namespace Tpetra 
{
  template<class OrdinalType, class ScalarType>
  class CrsMatrix 
  {
    public:
      CrsMatrix(VectorSpace<OrdinalType, ScalarType>& VectorRowSpace) :
        VectorRowSpace_(VectorRowSpace),
        RowSpace_(VectorRowSpace.elementSpace()),
        NumMyRows_(RowSpace_.getNumMyElements()),
        NumGlobalRows_(RowSpace_.getNumGlobalElements()),
        MyGlobalElements_(RowSpace_.getMyGlobalElements())
      {
        Indices_.resize(NumMyRows_);
        Values_.resize(NumMyRows_);
      }

      void fillComplete()
      {
#ifdef HAVE_MPI
        // =============================== //
        // Part 0: send off-image elements //
        // =============================== //

        int NumImages = RowSpace_.comm().getNumImages();
        map<OrdinalType, OrdinalType> containter_map;

        // this is a list of non-locally owned rows, in a map (should become a
        // hash some day for faster access)
        for (typename map<OrdinalType, vector<pair<OrdinalType, ScalarType> > >::iterator iter = nonlocals_.begin() ; 
             iter != nonlocals_.end() ; ++iter)
        {
          containter_map[iter->first] = 1;
        }

        // convert the map to a vector so that I can use get getRemoteIDList()
        vector<OrdinalType> container_vector;

        for (typename map<OrdinalType, OrdinalType>::iterator iter = containter_map.begin() ;
             iter != containter_map.end() ; ++iter)
        {
          container_vector.push_back(iter->first);
        }

        vector<OrdinalType> image_vector(container_vector.size());

        RowSpace_.getRemoteIDList (container_vector, image_vector);

        map<OrdinalType, OrdinalType> image_map;

        for (OrdinalType i = 0 ; i < image_vector.size() ; ++i)
        {
          image_map[container_vector[i]] = image_vector[i];
        }

        vector<OrdinalType> local_neighbors(RowSpace_.comm().getNumImages());
        for (OrdinalType i = 0 ; i < local_neighbors.size() ; ++i) local_neighbors[i] = 0;

        for (OrdinalType i = 0 ; i < image_vector.size() ; ++i)
        {
          local_neighbors[image_vector[i]] = 1;
        }

        vector<OrdinalType> global_neighbors(NumImages * NumImages);

        RowSpace_.comm().gatherAll(&local_neighbors[0], &global_neighbors[0], NumImages);

        // `global_neighbors' at this point contains (on all images) the
        // connectivity between the images. On the row `i', a nonzero on col `j' means
        // that image i will send something to image j. On the column `j', a
        // nonzero on row `i' means that image j will receive something from
        // image i.
        
        // now I loop over all columns to know which image is supposed to send
        // me something
        vector<OrdinalType> recvImages;

        for (int j = 0 ; j < NumImages ; ++j)
        {
          OrdinalType what = global_neighbors[j * NumImages + RowSpace_.comm().getMyImageID()];
          if (what > 0)
          {
            recvImages.push_back(j);
          }
        }

        // do the same but with send
        vector<OrdinalType> sendImages;

        // now I pack what has to be sent to the other images
        map<OrdinalType, vector<OrdinalType> > sendRows;
        map<OrdinalType, vector<OrdinalType> > sendCols;
        map<OrdinalType, vector<ScalarType> >  sendVals;

        for (typename map<OrdinalType, vector<pair<OrdinalType, ScalarType> > >::iterator iter = nonlocals_.begin() ; 
             iter != nonlocals_.end() ; ++iter)
        {
          OrdinalType row   = iter->first;
          OrdinalType image = image_map[row];

          for (OrdinalType i = 0 ; i < iter->second.size() ; ++i)
          {
            OrdinalType col   = iter->second[i].first;
            ScalarType  val   = iter->second[i].second;

            sendRows[image].push_back(row);
            sendCols[image].push_back(col);
            sendVals[image].push_back(val);
          }
        }

        int MyImageID = RowSpace_.comm().getMyImageID();

        vector<MPI_Request> send_requests(NumImages * 3);

        OrdinalType send_count = 0;
        for (int j = 0 ; j < NumImages ; ++j)
        {
          OrdinalType what = global_neighbors[j + NumImages * RowSpace_.comm().getMyImageID()];
          if (what > 0)
          {
            sendImages.push_back(j);
            OrdinalType size = sendRows[j].size();

            MPI_Isend(&size, 1, MPI_INT, j, 23, MPI_COMM_WORLD, &(send_requests[send_count]));
            ++send_count;
          }
        }

        // Now receive the actual sizes
        vector<MPI_Request> recv_requests(NumImages * 3);
        vector<OrdinalType> recv_sizes(NumImages);
        vector<OrdinalType> recv_images(NumImages);

        OrdinalType recv_count = 0;
        for (int j = 0 ; j < NumImages ; ++j)
        {
          OrdinalType what = global_neighbors[j * NumImages + RowSpace_.comm().getMyImageID()];
          if (what > 0)
          {
            recv_images[recv_count] = j;
            MPI_Irecv(&(recv_sizes[recv_count]), 1, MPI_INT, j, 23, MPI_COMM_WORLD, &(recv_requests[recv_count]));
            ++recv_count;
          }
        }

        for (int i = 0 ; i < send_count ; ++i)
          MPI_Wait(&(send_requests[i]), MPI_STATUS_IGNORE);

        for (int i = 0 ; i < recv_count ; ++i)
          MPI_Wait(&(recv_requests[i]), MPI_STATUS_IGNORE);

        MPI_Barrier(MPI_COMM_WORLD);

        map<OrdinalType, vector<OrdinalType> > recvRows;
        map<OrdinalType, vector<OrdinalType> > recvCols;
        map<OrdinalType, vector<ScalarType> >  recvVals;

        map<OrdinalType, OrdinalType> xxx;

        for (int i = 0 ; i < recv_count ; ++i)
        {
          OrdinalType image = recv_images[i];
          recvRows[image].resize(recv_sizes[i]);
          recvCols[image].resize(recv_sizes[i]);
          recvVals[image].resize(recv_sizes[i]);

          xxx[image] = recv_sizes[i];
        }

        // At this point I know:
        // - I have to receive from `recv_count' images;
        // - image `i' will send recv_count[i] things, split in
        //   two vectors of OrdinalType and a vector of ScalarType.
        // First I start sending, then receiving

        send_count = 0;
        for (int j = 0 ; j < NumImages ; ++j)
        {
          OrdinalType what = global_neighbors[j + NumImages * RowSpace_.comm().getMyImageID()];
          if (what > 0)
          {
            // want to send to image `j', first Rows, then Cols, then Vals
            OrdinalType size = sendRows[j].size();
            
            MPI_Isend(&(sendRows[j][0]), size, MPI_INT, j, 32, MPI_COMM_WORLD, &(send_requests[send_count]));
            ++send_count;

            MPI_Isend(&(sendCols[j][0]), size, MPI_INT, j, 33, MPI_COMM_WORLD, &(send_requests[send_count]));
            ++send_count;

            MPI_Isend(&(sendVals[j][0]), size, MPI_DOUBLE, j, 34, MPI_COMM_WORLD, &(send_requests[send_count]));
            ++send_count;
          }
        }

        recv_count = 0;
        for (int j = 0 ; j < NumImages ; ++j)
        {
          OrdinalType what = global_neighbors[j * NumImages + RowSpace_.comm().getMyImageID()];
          if (what > 0)
          {
            OrdinalType size = xxx[j];

            // I want to receive from image `j', first Rows, then Cols, then Vals.
            MPI_Irecv(&(recvRows[j][0]), size, MPI_INT, j, 32, MPI_COMM_WORLD, &(recv_requests[recv_count]));
            ++recv_count;

            MPI_Irecv(&(recvCols[j][0]), size, MPI_INT, j, 33, MPI_COMM_WORLD, &(recv_requests[recv_count]));
            ++recv_count;

            MPI_Irecv(&(recvVals[j][0]), size, MPI_DOUBLE, j, 34, MPI_COMM_WORLD, &(recv_requests[recv_count]));
            ++recv_count;
          }
        }

        for (int i = 0 ; i < send_count ; ++i)
          MPI_Wait(&(send_requests[i]), MPI_STATUS_IGNORE);

        for (int i = 0 ; i < recv_count ; ++i)
          MPI_Wait(&(recv_requests[i]), MPI_STATUS_IGNORE);

        MPI_Barrier(MPI_COMM_WORLD);

        // now I add all the received elements to the list of local elements.

        for (int i = 0 ; i < recv_images.size() ; ++i)
        {
          OrdinalType image = recv_images[i];
          for (int j = 0 ; j < recv_sizes[i] ; ++j)
          {
            setGlobalElement(recvRows[image][j], recvCols[image][j], recvVals[image][j]);
          }
        }
#endif

        // =============================== //
        // Part I: remove repeated indices //
        // =============================== //
        
        // I load all matrix entries in a hash table, then I re-fill
        // the row with the last inserted value.
        for (OrdinalType i = 0 ; i < NumMyRows_ ; ++i)
        {
          std::map<OrdinalType, ScalarType> singleRow;

          for (OrdinalType j = 0 ; j < Indices_[i].size() ; ++j)
          {
            singleRow[Indices_[i][j]] = Values_[i][j];
          }

          OrdinalType count = 0;
          for (typename std::map<OrdinalType,ScalarType>::iterator iter = singleRow.begin() ; 
               iter != singleRow.end() ; ++iter)
          {
            Indices_[i][count] = iter->first;
            Values_[i][count] = iter->second;
            ++count;
          }

          Indices_[i].resize(count);
          Values_[i].resize(count);
        }

        // =============================== //
        // Part II: build the column space //
        // =============================== //
        
        // I have to find the list of non-locally owned columns

        map<OrdinalType, bool> container; // replace with a hash table

        for (OrdinalType i = 0 ; i < NumMyRows_ ; ++i)
        {
          for (OrdinalType j = 0 ; j < Indices_[i].size() ; ++j)
          {
            OrdinalType what = Indices_[i][j];
            if (RowSpace_.isMyGID(what)) continue;
            else
              container[what] = true;
          }
        }

        vector<OrdinalType> MyPaddedGlobalElements(MyGlobalElements_);

        for (typename std::map<OrdinalType, bool>::iterator iter = container.begin() ; 
             iter != container.end() ; ++iter)
        {
          MyPaddedGlobalElements.push_back(iter->first);
        }

        // now I can build the column space

        ColSpace_ = new ElementSpace<OrdinalType>(-1, MyPaddedGlobalElements.size(),
                                                  MyPaddedGlobalElements, RowSpace_.getIndexBase(), RowSpace_.platform());
        VectorColSpace_ = new VectorSpace<OrdinalType, ScalarType>(*ColSpace_, VectorRowSpace_.platform());

        PaddedVector_ = new Vector<OrdinalType, ScalarType>(*VectorColSpace_);
        Importer_ = new Import<OrdinalType>(RowSpace_, *ColSpace_);

        // =============================== //
        // Part III: move to local indices //
        // =============================== //
        
        for (OrdinalType i = 0 ; i < NumMyRows_ ; ++i)
        {
          for (OrdinalType j = 0 ; j < Indices_[i].size() ; ++j)
          {
            Indices_[i][j] = ColSpace_->getLID(Indices_[i][j]);
          }
        }
      }

      void setGlobalElement(const OrdinalType& GlobalRow, const OrdinalType& GlobalCol,
                            const ScalarType& Value)
      {
        if (RowSpace_.isMyGID(GlobalRow))
        {
          OrdinalType LocalRow = RowSpace_.getLID(GlobalRow);
          Indices_[LocalRow].push_back(GlobalCol);
          Values_[LocalRow].push_back(Value);
        }
        else
        {
          nonlocals_[GlobalRow].push_back(pair<OrdinalType, ScalarType>(GlobalCol, Value));
        }
      }

      void setGlobalElements(std::vector<OrdinalType>& GlobalRows, 
                             std::vector<OrdinalType>& GlobalCols,
                             std::vector<ScalarType>& Values)
      {
        if (GlobalRows.size() != GlobalCols.size())
          throw(-1);

        for (OrdinalType i = 0 ; i < GlobalRows.size() ; ++i)
        {
          OrdinalType LocalRow = RowSpace_.getLID(GlobalRows[i]);
          Indices_[LocalRow].push_back(GlobalCols[i]);
          Values_[LocalRow].push_back(Values[i]);
        }
      }

      void getMyRowCopy(const OrdinalType MyRow, vector<OrdinalType> Indices,
                        vector<ScalarType> Values)
      {
        OrdinalType length = Indices_[MyRow].size();

        if (Indices.size() < length || Values.size() < length)
          throw(-1);

        for (OrdinalType i = 0 ; i < length ; ++i)
        {
          Indices[i] = Indices_[MyRow][i];
          Values[i] = Values_[MyRow][i];
        }
      }

      void apply(const Vector<OrdinalType, ScalarType>& x,
                 Vector<OrdinalType, ScalarType> y) const
      {
        y.setAllToScalar(0.0);

        PaddedVector_->doImport(x, *Importer_, Insert);

        for (OrdinalType i = 0 ; i < NumMyRows_ ; ++i)
        {
          for (OrdinalType j = 0 ; j < Indices_[i].size() ; ++j)
          {
            OrdinalType col = Indices_[i][j];
            ScalarType val = Values_[i][j];
            y[i] += val * (*PaddedVector_)[col];
          }
        }
      }

    private:
      const vector<vector<OrdinalType> >& getIndices() const
      {
        return(Indices_);
      }

      const vector<vector<ScalarType> >& getValues() const
      {
        return(Values_);
      }

      const VectorSpace<OrdinalType, ScalarType>& VectorRowSpace_;
      const ElementSpace<OrdinalType>& RowSpace_;

      OrdinalType NumMyRows_, NumGlobalRows_;
      const vector<OrdinalType>& MyGlobalElements_;

      std::vector<vector<OrdinalType> > Indices_;
      std::vector<vector<ScalarType> > Values_;

      ElementSpace<OrdinalType>* ColSpace_;
      VectorSpace<OrdinalType, ScalarType>* VectorColSpace_;

      Vector<OrdinalType, ScalarType>* PaddedVector_;
      Import<OrdinalType>* Importer_;

      map<OrdinalType, vector<pair<OrdinalType, ScalarType> > > nonlocals_;
  }; // class CrsMatrix
} // namespace Tpetra



typedef int OrdinalType;
typedef double ScalarType;

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Tpetra::MpiComm<OrdinalType, ScalarType> Comm(MPI_COMM_WORLD);
#else
  Tpetra::SerialComm<OrdinalType, ScalarType> Comm;
#endif

  // Get zero and one for the OrdinalType
  
  OrdinalType const OrdinalZero = Teuchos::ScalarTraits<OrdinalType>::zero();
  OrdinalType const OrdinalOne  = Teuchos::ScalarTraits<OrdinalType>::one();

  // Get zero and one for the ScalarType
  
  ScalarType const ScalarZero = Teuchos::ScalarTraits<ScalarType>::zero();
  ScalarType const ScalarOne  = Teuchos::ScalarTraits<ScalarType>::one();

  // Creates a vector of size `length', then set the elements values.
  
  OrdinalType length    = OrdinalOne * 4 * Comm.getNumImages();
  OrdinalType indexBase = OrdinalZero;

  // 1) Creation of a platform
  
#ifdef HAVE_MPI
  const Tpetra::MpiPlatform <OrdinalType, OrdinalType> platformE(MPI_COMM_WORLD);
  const Tpetra::MpiPlatform <OrdinalType, ScalarType> platformV(MPI_COMM_WORLD);
#else
  const Tpetra::SerialPlatform <OrdinalType, OrdinalType> platformE;
  const Tpetra::SerialPlatform <OrdinalType, ScalarType> platformV;
#endif

  // 2) We can now create a space:

  Tpetra::ElementSpace<OrdinalType> elementSpace(length, indexBase, platformE);

  Tpetra::VectorSpace<OrdinalType, ScalarType> vectorSpace(elementSpace, platformV);

  // 3) and the vector, which has type int for the OrdinalType
  //    and double for the ScalarType

  Tpetra::CrsMatrix<OrdinalType, ScalarType> Matrix(vectorSpace);

  if (Comm.getMyImageID() == 0)
  {
    for (int i = 0 ; i < elementSpace.getNumGlobalElements() ; ++i)
    {
      int GRID = i; // elementSpace.getGID(i);

      if (GRID != 0)
        Matrix.setGlobalElement(GRID, GRID - 1, -1.0);
      if (GRID != elementSpace.getNumGlobalElements() - 1)
        Matrix.setGlobalElement(GRID, GRID + 1, -1.0);
      Matrix.setGlobalElement(GRID, GRID, 2.0);
    }
  }

  Tpetra::Vector<OrdinalType, ScalarType> x(vectorSpace), y(vectorSpace);

  x.setAllToScalar(1.0);

  Matrix.fillComplete();
  Matrix.apply(x, y);

  cout << y;

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);
}
