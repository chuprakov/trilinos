/* True C header file! */


#ifndef FEPETRA_MAP_H
#define FEPETRA_MAP_H


#ifdef __cplusplus
extern "C" {
#endif


typedef int MapID;

MapID Epetra_Map_Create( int numGlobalElements );

void Epetra_Map_Destroy( MapID mapID );

int Epetra_Map_NumGlobalElements( MapID mapID );

int Epetra_Map_NumMyElements( MapID mapID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* FEPETRA_MAP_H */
