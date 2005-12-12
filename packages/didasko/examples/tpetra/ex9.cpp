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
        // =============================== //
        // Part 0: send off-image elements //
        // =============================== //

        int NumImages = RowSpace_.comm().getNumImages();
        map<OrdinalType, OrdinalType> containter;

        for (typename map<OrdinalType, vector<OrdinalType> >::iterator iter = nonlocalIndices_.begin() ; 
             iter != nonlocalIndices_.end() ; ++iter)
        {
          containter[iter->first] = 1;
        }

        // convert the map to a vector so that I can use get getRemoteIDList()
        vector<OrdinalType> myrowlist;

        for (typename map<OrdinalType, OrdinalType>::iterator iter = containter.begin() ;
             iter != containter.end() ; ++iter)
        {
          myrowlist.push_back(iter->first);
        }

        vector<OrdinalType> imageIDList(myrowlist.size());

        RowSpace_.getRemoteIDList (myrowlist, imageIDList);

        vector<OrdinalType> mylist(RowSpace_.comm().getNumImages());
        for (OrdinalType i = 0 ; i < mylist.size() ; ++i) mylist[i] = 0;

        for (OrdinalType i = 0 ; i < myrowlist.size() ; ++i)
        {
          mylist[imageIDList[i]] = 1;
        }

        vector<OrdinalType> globalList(NumImages * NumImages);

        RowSpace_.comm().gatherAll(&mylist[0], &globalList[0], NumImages);

        // now I loop over all columns to know which image is supposed to send
        // me something
        vector<OrdinalType> rcvImages;

        for (int j = 0 ; j < NumImages ; ++j)
        {
          OrdinalType what = globalList[j * NumImages + RowSpace_.comm().getMyImageID()];
          if (what > 0)
          {
            rcvImages.push_back(j);
          }
        }

        for (int i = 0 ; i < rcvImages.size() ; ++i)
          cout << rcvImages[i] << " " << endl;

#if 0
        // At this point I know how many processors will send a message to
        // this image. I have to compute what to send, size and content.
        // I make another loop over all the elements, by separing all data
        // image by image.
        
        map<OrdinalType, OrdinalType> sendSizes;
        map<OrdinalType, vector<OrdinalType> > sendIndices;
        map<OrdinalType, vector<ScalarType> > sendValues;

        for (typename map<OrdinalType, vector<OrdinalType> >::iterator iter = nonlocalIndices_.begin() ; 
             iter != nonlocalIndices_.end() ; ++iter)
        {
          OrdinalType Row = iter->first;
          for (int j = 0 ; j < (iter->second).size() ; ++j)
          {
            OrdinalType Col = (iter->second)[j];
            ScalarType Val = 
          containter[iter->first] = 1;
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

      void setGlobalElement(OrdinalType GlobalRow, OrdinalType GlobalCol,
                            ScalarType Value)
      {
        if (RowSpace_.isMyGID(GlobalRow))
        {
          OrdinalType LocalRow = RowSpace_.getLID(GlobalRow);
          Indices_[LocalRow].push_back(GlobalCol);
          Values_[LocalRow].push_back(Value);
        }
        else
        {
          nonlocalIndices_[GlobalRow].push_back(GlobalCol);
          nonlocalValues_[GlobalRow].push_back(Value);
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

      map<OrdinalType, vector<OrdinalType> > nonlocalIndices_;
      map<OrdinalType, vector<ScalarType> > nonlocalValues_;
  };
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

#if 0
  Tpetra::Vector<OrdinalType, ScalarType> x(vectorSpace), y(vectorSpace);

  x.setAllToScalar(1.0);

#endif
  Matrix.fillComplete();
#if 0
  Matrix.apply(x, y);

  cout << y;
#endif



#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);
}
