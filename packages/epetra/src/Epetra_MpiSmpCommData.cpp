
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

#include "Epetra_MpiSmpCommData.h"
//#include "Epetra_ConfigDefs.h" //DATA_DEBUG

//=============================================================================
Epetra_MpiSmpCommData::Epetra_MpiSmpCommData(MPI_Comm& Comm) 
	: Comm_(Comm),
		curTag_(minTag_),
		ThreadID_(0),
		NumThreads_(1)
{
	MPI_Comm_size(Comm, &size_);
	MPI_Comm_rank(Comm, &rank_);
	minTag_ = 24050;
	maxTag_ = 24099;
	NodeID_ = rank_;
	//cout << "--MSCD created (dc), addr: " << this << endl; //DATA_DEBUG
}

//=============================================================================
Epetra_MpiSmpCommData::~Epetra_MpiSmpCommData() {
	//cout << "--MSCD destroyed, addr: " << this << endl; //DATA_DEBUG
}
