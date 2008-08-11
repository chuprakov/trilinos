// @HEADER
// ************************************************************************
//
//            Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

#ifndef GALERI_VIZ_MEDIT_H
#define GALERI_VIZ_MEDIT_H

#include "Galeri_grid_Element.h"

#include "Epetra_Import.h"

namespace Galeri {
namespace viz {

class MEDIT
{
public:

  static
  void write(Galeri::grid::Loadable& patch,
             const string& BaseName,
             const Epetra_MultiVector& vector)
  {
    const Epetra_Comm& comm = patch.getComm();

    std::ofstream medit;
    string FileName = BaseName + ".mesh";

    const Epetra_MultiVector& linearCoord = patch.getLinearCoordinates();

    if (comm.MyPID() == 0)
    {
      medit.open(FileName.c_str());
      medit << "MeshVersionFormatted 1" << endl;
      medit << "Dimension 3" << endl;
      medit << "# mesh from phoenix" << endl << endl;
      medit << "Vertices " << patch.getNumGlobalVertices() << endl;
      medit.close();
    }

    for (int iproc = 0; iproc < comm.NumProc(); ++iproc)
    {
      if (iproc == comm.MyPID())
      {
        medit.open(FileName.c_str(),ios::app);
        if (Galeri::core::Workspace::getNumDimensions() == 2)
        {
          for (int i = 0; i < linearCoord.MyLength(); ++i)
          {
            medit << setw(12) << setiosflags(ios::showpoint) 
              << setw(12) << linearCoord[0][i] << " "
              << setw(12) << linearCoord[1][i] << " "
              << setw(12) << "0.0 1" << endl;
          }
        }
        else if (Galeri::core::Workspace::getNumDimensions() == 3)
        {
          for (int i = 0; i < linearCoord.MyLength(); ++i)
          {
            medit << setw(12) << setiosflags(ios::showpoint) 
              << setw(12) << linearCoord[0][i] << " "
              << setw(12) << linearCoord[1][i] << " "
              << setw(12) << linearCoord[2][i] << " 1" << endl;
          }
        }
          
        medit.close();
      }
      comm.Barrier();
    }

    for (int iproc = 0; iproc < comm.NumProc(); ++iproc)
    {
      if (comm.MyPID() == iproc) 
      {
        medit.open(FileName.c_str(),ios::app);
        if (iproc == 0)
        {
          if (patch.getElement().getLabel() == "Galeri::grid::Triangle")
            medit << "Triangles " << patch.getNumGlobalElements() << endl;
          else if (patch.getElement().getLabel() == "Galeri::grid::Quad")
            medit << "Quadrilaterals " << patch.getNumGlobalElements() << endl;
          else if (patch.getElement().getLabel() == "Galeri::grid::Hex")
            medit << "Hexahedra " << patch.getNumGlobalElements() << endl;
          else
          {
            cout << "NOT SUPPORTED YET" << endl;
            exit(0);
          }
        }

        for (int i = 0 ; i < patch.getNumMyElements(); ++i) 
        {
          for (int j = 0; j < patch.getNumVerticesPerElement(); ++j)
            medit << patch.getMyConnectivity(i, j) + 1 << " ";

          medit << comm.MyPID() << endl;
        }

        if (iproc == comm.NumProc() - 1) 
          medit << endl << "End" << endl;

        medit.close();
      }

      comm.Barrier();
    }

    Epetra_Import importer(linearCoord.Map(), vector.Map());
    Epetra_MultiVector linearVector(linearCoord.Map(), vector.NumVectors());
    linearVector.Import(vector, importer, Insert);

    // ======== //
    // .bb file //
    // ======== //

    FileName = BaseName + ".bb";    

    for (int iproc = 0; iproc < comm.NumProc(); ++iproc)
    {
      if (comm.MyPID() == iproc) 
      {
        if (iproc == 0)
        {
          medit.open(FileName.c_str());
          medit << "3 1 " << patch.getNumGlobalVertices() << " 2" << endl;
        }
        else
          medit.open(FileName.c_str(),ios::app);

        for (int i = 0; i < linearVector.MyLength(); ++i)
          medit << setiosflags(ios::showpoint) << linearVector[0][i] << endl;

        medit.close();
      }
      comm.Barrier();
    }
  }
};

} // namespace viz
} // namespace Galeri

#endif
