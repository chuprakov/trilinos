//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
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

#include "checkmap.h"
int checkmap(Epetra_BlockMap & Map, int NumGlobalElements, int NumMyElements, 
	     int *MyGlobalElements, int ElementSize, int * ElementSizeList,
	     int NumGlobalPoints, int NumMyPoints,
	     int IndexBase, Epetra_Comm& Comm,
	     bool DistributedGlobal,
	     bool IsOneToOne)
{

  int i, ierr=0, forierr=0;// forierr is used in for loops, then is tested
  // after for loop completes to see if it is non zero - potentially prevents
  // thousands of error messages

  if (ElementSizeList==0)
    {
      EPETRA_TEST_ERR(!Map.ConstantElementSize(),ierr);
    }
  else
    EPETRA_TEST_ERR(Map.ConstantElementSize(),ierr);
  
  EPETRA_TEST_ERR(DistributedGlobal!=Map.DistributedGlobal(),ierr);

  EPETRA_TEST_ERR(IsOneToOne!=Map.IsOneToOne(),ierr);

  int *MyElementSizeList;

  if (ElementSizeList==0)
    {
      EPETRA_TEST_ERR(Map.ElementSize()!=ElementSize,ierr);
      
      MyElementSizeList = new int[NumMyElements];
      
      EPETRA_TEST_ERR(Map.ElementSizeList(MyElementSizeList)!=0,ierr);
      forierr = 0;
      for (i=0; i<NumMyElements; i++) 
        forierr += MyElementSizeList[i]!=ElementSize;
      EPETRA_TEST_ERR(forierr,ierr);

      EPETRA_TEST_ERR(Map.MaxMyElementSize() != ElementSize,ierr);
      EPETRA_TEST_ERR(Map.MinMyElementSize() != ElementSize,ierr);
    }
  else
    {
      MyElementSizeList = new int[NumMyElements];
      EPETRA_TEST_ERR(Map.ElementSizeList(MyElementSizeList)!=0,ierr);
      int MaxSize = MyElementSizeList[0];
      int MinSize = MyElementSizeList[0];
      forierr=0;
      for (i=0; i<NumMyElements; i++) {
        forierr += MyElementSizeList[i]!=ElementSizeList[i];
	if (MyElementSizeList[i] > MaxSize)
	  MaxSize = MyElementSizeList[i];
	if (MyElementSizeList[i] < MinSize)
	  MinSize = MyElementSizeList[i];

	// Test ElementSize(int LID) method	

	forierr += Map.ElementSize(Map.LID(MyGlobalElements[i])) != ElementSizeList[i];
      }
      EPETRA_TEST_ERR(forierr,ierr);
   
      EPETRA_TEST_ERR(MaxSize !=Map.MaxMyElementSize(),ierr);
      EPETRA_TEST_ERR(MinSize !=Map.MinMyElementSize(),ierr);
    }

  const Epetra_Comm & Comm1 = Map.Comm();

  EPETRA_TEST_ERR(Comm1.NumProc()!=Comm.NumProc(),ierr);

  EPETRA_TEST_ERR(Comm1.MyPID()!=Comm.MyPID(),ierr);

  EPETRA_TEST_ERR(Map.IndexBase()!=IndexBase,ierr);

  EPETRA_TEST_ERR(!Map.LinearMap() && MyGlobalElements==0,ierr);

  EPETRA_TEST_ERR(Map.LinearMap() && MyGlobalElements!=0,ierr);

  EPETRA_TEST_ERR(Map.MaxAllGID()!=NumGlobalElements-1+IndexBase,ierr);

  EPETRA_TEST_ERR(Map.MaxElementSize()!=ElementSize,ierr);

  int MaxLID = Map.MaxLID();
  EPETRA_TEST_ERR(MaxLID!=NumMyElements-1,ierr);

  int MaxMyGID = (Comm.MyPID()+1)*NumMyElements-1+IndexBase;
  if (Comm.MyPID()>2) MaxMyGID+=3;
  if (!DistributedGlobal) MaxMyGID = NumMyElements-1+IndexBase;
  EPETRA_TEST_ERR(Map.MaxMyGID()!=MaxMyGID,ierr);

  EPETRA_TEST_ERR(Map.MinAllGID()!=IndexBase,ierr);

  if (ElementSizeList==0)
    {
      EPETRA_TEST_ERR(Map.MinElementSize()!=ElementSize,ierr);
    }
  else EPETRA_TEST_ERR(Map.MinElementSize()!=2,ierr);

  int MinLID = Map.MinLID();
  EPETRA_TEST_ERR(MinLID!=0,ierr);

  int MinMyGID = Comm.MyPID()*NumMyElements+IndexBase;
  if (Comm.MyPID()>2) MinMyGID+=3;
  if (!DistributedGlobal) MinMyGID = IndexBase; // Not really needed
  EPETRA_TEST_ERR(Map.MinMyGID()!=MinMyGID,ierr);

  int * MyGlobalElements1 = new int[NumMyElements];
  EPETRA_TEST_ERR(Map.MyGlobalElements(MyGlobalElements1)!=0,ierr);
  
  forierr = 0;
  if (MyGlobalElements==0) {
    for (i=0; i<NumMyElements; i++) 
      forierr += MyGlobalElements1[i]!=MinMyGID+i;
    EPETRA_TEST_ERR(forierr,ierr);
  }
  else {
    for (i=0; i<NumMyElements; i++)
      forierr += MyGlobalElements[i]!=MyGlobalElements1[i];
    EPETRA_TEST_ERR(forierr,ierr);
  }
  EPETRA_TEST_ERR(Map.NumGlobalElements()!=NumGlobalElements,ierr);
  
  EPETRA_TEST_ERR(Map.NumGlobalPoints()!=NumGlobalPoints,ierr);
  
  EPETRA_TEST_ERR(Map.NumMyElements()!=NumMyElements,ierr);  

  EPETRA_TEST_ERR(Map.NumMyPoints()!=NumMyPoints,ierr);

  int MaxMyGID2 = Map.GID(Map.LID(MaxMyGID));
  EPETRA_TEST_ERR(MaxMyGID2 != MaxMyGID,ierr);
  int MaxLID2 = Map.LID(Map.GID(MaxLID));
  EPETRA_TEST_ERR(MaxLID2 != MaxLID,ierr);

  EPETRA_TEST_ERR(Map.GID(MaxLID+1) != IndexBase-1,ierr);// MaxLID+1 doesn't exist
  EPETRA_TEST_ERR(Map.LID(MaxMyGID+1) != -1,ierr);// MaxMyGID+1 doesn't exist or is on a different processor

  EPETRA_TEST_ERR(!Map.MyGID(MaxMyGID),ierr);
  EPETRA_TEST_ERR(Map.MyGID(MaxMyGID+1),ierr);

  EPETRA_TEST_ERR(!Map.MyLID(MaxLID),ierr);
  EPETRA_TEST_ERR(Map.MyLID(MaxLID+1),ierr);

  EPETRA_TEST_ERR(!Map.MyGID(Map.GID(MaxLID)),ierr);
  EPETRA_TEST_ERR(Map.MyGID(Map.GID(MaxLID+1)),ierr);

  EPETRA_TEST_ERR(!Map.MyLID(Map.LID(MaxMyGID)),ierr);
  EPETRA_TEST_ERR(Map.MyLID(Map.LID(MaxMyGID+1)),ierr);

  // Test the FirstPointInElementList methods, begin by testing that they produce identical results
  int * FirstPointInElementList = new int[NumMyElements+1];
  Map.FirstPointInElementList(FirstPointInElementList);
  int * FirstPointInElementList1 = Map.FirstPointInElementList();
  forierr = 0;
  for (i=0; i<=NumMyElements; i++)
    forierr += FirstPointInElementList[i]!=FirstPointInElementList1[i];
  EPETRA_TEST_ERR(forierr,ierr);
  // Now make sure values are correct
  forierr = 0;
  if (Map.ConstantElementSize()) {
    for (i=0; i<=NumMyElements; i++)
      forierr += FirstPointInElementList1[i]!=(i*ElementSize);// NOTE:FirstPointInElement[NumMyElements] is not the first point of an element
    EPETRA_TEST_ERR(forierr,ierr);
  }
  else {
    int FirstPoint = 0;
    for (i=0; i<NumMyElements; i++) {
      forierr += FirstPointInElementList1[i]!=FirstPoint;
      FirstPoint += ElementSizeList[i];
    }
    EPETRA_TEST_ERR(forierr,ierr);
    EPETRA_TEST_ERR(FirstPointInElementList[NumMyElements] != NumMyPoints,ierr);// The last entry in the array = the total number of Points on the proc
  }
  delete [] FirstPointInElementList;

  // Declare some variables for the FindLocalElementID test
  int ElementID, Offset;
  // Test the PointToElementList methods, begin by testing that they produce identical results
  int * PointToElementList = new int[NumMyPoints];
  Map.PointToElementList(PointToElementList);
  int * PointToElementList1 = Map.PointToElementList();
  forierr = 0;
  for (i=0; i<NumMyPoints; i++)
    forierr += PointToElementList1[i] != PointToElementList[i];
  EPETRA_TEST_ERR(forierr,ierr);
  //Now make sure values are correct
  forierr=0;
  if (Map.ConstantElementSize()) {
    for (i=0; i<NumMyElements; i++)
      for (int j=0; j<ElementSize; j++) {
	forierr += PointToElementList[i*ElementSize+j] != i;
	// Test FindLocalElementID method
	Map.FindLocalElementID(i*ElementSize+j,ElementID,Offset);
	forierr += ElementID != i || Offset != j;
      }
    EPETRA_TEST_ERR(forierr,ierr);
  }
  else {
    int MyPointTot = 0; // Keep track of total number of points in all previously completely checked elements
    for (i=0; i<NumMyElements; i++) {
      for (int j=0; j<ElementSizeList[i]; j++) {
	forierr += PointToElementList[MyPointTot+j] != i;
	// Test FindLocalElementID method
	Map.FindLocalElementID(MyPointTot+j,ElementID,Offset);
	forierr += ElementID != i || Offset != j;
      }
      MyPointTot += ElementSizeList[i];
    }
    EPETRA_TEST_ERR(forierr,ierr);
  }
  delete [] PointToElementList;

  // Check RemoteIDList function that includes a parameter for size
  // Get some GIDs off of each processor to test
  int TotalNumEle, NumElePerProc, NumProc = Comm.NumProc();
  int MinNumEleOnProc;
  int NumMyEle = Map.NumMyElements();
  Comm.MinAll(&NumMyEle,&MinNumEleOnProc,1);
  if (MinNumEleOnProc > 5) NumElePerProc = 6;
  else NumElePerProc = MinNumEleOnProc;
  if (NumElePerProc > 0) {
    TotalNumEle = NumElePerProc*NumProc;
    int * MyGIDlist = new int[NumElePerProc];
    int * GIDlist = new int[TotalNumEle];
    int * PIDlist = new int[TotalNumEle];
    int * LIDlist = new int[TotalNumEle];
    int * SizeList = new int[TotalNumEle];
    for (i=0; i<NumElePerProc; i++)
	  MyGIDlist[i] = MyGlobalElements1[i];
    Comm.GatherAll(MyGIDlist,GIDlist,NumElePerProc);// Get a few values from each proc
    Map.RemoteIDList(TotalNumEle, GIDlist, PIDlist, LIDlist, SizeList);
    int MyPID= Comm.MyPID();
    forierr = 0;
    for (i=0; i<TotalNumEle; i++) {
      if (Map.MyGID(GIDlist[i])) {
	forierr += PIDlist[i] != MyPID;
	forierr += !Map.MyLID(Map.LID(GIDlist[i])) || Map.LID(GIDlist[i]) != LIDlist[i] || Map.GID(LIDlist[i]) != GIDlist[i];
	forierr += SizeList[i] != Map.ElementSize(LIDlist[i]);
      }
      else {
	forierr += PIDlist[i] == MyPID; // If MyGID comes back false, the PID listed should be that of another proc
      }
    }
    EPETRA_TEST_ERR(forierr,ierr);

    delete [] MyGIDlist;
    delete [] GIDlist;
    delete [] PIDlist;
    delete [] LIDlist;
    delete [] SizeList;
  }

  delete [] MyGlobalElements1;
  delete [] MyElementSizeList;

  // Check RemoteIDList function (assumes all maps are linear, even if not stored that way)

  if (Map.LinearMap()) {

    int * GIDList = new int[3];
    int * PIDList = new int[3];
    int * LIDList = new int[3];
    int MyPID = Map.Comm().MyPID();
  
    int NumIDs = 0;
    //GIDList[NumIDs++] = Map.MaxAllGID()+1; // Should return -1 for both PID and LID
    if (Map.MinMyGID()-1>=Map.MinAllGID()) GIDList[NumIDs++] = Map.MinMyGID()-1;
    if (Map.MaxMyGID()+1<=Map.MaxAllGID()) GIDList[NumIDs++] = Map.MaxMyGID()+1;

    Map.RemoteIDList(NumIDs, GIDList, PIDList, LIDList);

    NumIDs = 0;

    //EPETRA_TEST_ERR(!(PIDList[NumIDs]==-1),ierr);
    //EPETRA_TEST_ERR(!(LIDList[NumIDs++]==-1),ierr);

    if (Map.MinMyGID()-1>=Map.MinAllGID()) EPETRA_TEST_ERR(!(PIDList[NumIDs++]==MyPID-1),ierr);
    if (Map.MaxMyGID()+1<=Map.MaxAllGID()) EPETRA_TEST_ERR(!(PIDList[NumIDs]==MyPID+1),ierr);
    if (Map.MaxMyGID()+1<=Map.MaxAllGID()) EPETRA_TEST_ERR(!(LIDList[NumIDs++]==0),ierr);

    delete [] GIDList;
    delete [] PIDList;
    delete [] LIDList;


  }
  return (ierr);
}
