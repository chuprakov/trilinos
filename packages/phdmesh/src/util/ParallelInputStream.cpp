/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards
 */

#include <cstdio>
#include <stdexcept>
#include <util/ParallelInputStream.hpp>

/*--------------------------------------------------------------------*/

namespace phdmesh {
namespace {

#if defined( PHDMESH_HAS_MPI )

void broadcast( ParallelMachine comm , void * buf , int n )
{ MPI_Bcast( buf , n , MPI_BYTE , 0 , comm ); }

#else

void broadcast( ParallelMachine , void * , int ) {}

#endif

}
}

/*--------------------------------------------------------------------*/

namespace phdmesh {
namespace {

//----------------------------------------------------------------------

class ParInBuf : public std::streambuf {
public:
  enum { BUFFER_LENGTH  = 0x010000 /* 64k bytes */ };
  enum { BUFFER_PUTBACK = 0x000010 /*  16 bytes */ };
  enum { MAX_READ       = BUFFER_LENGTH - BUFFER_PUTBACK };
  ParInBuf( ParallelMachine , const char * const );
  virtual ~ParInBuf();

protected:
   virtual int underflow(); // refill the input buffer
   virtual int overflow( int c = EOF ); // Not called
   virtual int sync(); // No-op
   virtual std::streambuf * setbuf( char * , std::streamsize ); // No-op

private:
  void close();

  ParallelMachine m_comm ;
  std::FILE     * m_root_fp ;
  char            m_buffer[ BUFFER_LENGTH ];
};

ParInBuf::ParInBuf( ParallelMachine comm , const char * const file_name )
  : m_comm( comm ), m_root_fp( NULL )
{
  int result = 1 ;

  if ( 0 == parallel_machine_rank( comm ) && NULL != file_name ) {
    result = NULL != ( m_root_fp = std::fopen( file_name , "r" ) );
  }

  broadcast( m_comm , & result , sizeof(int) );

  if ( ! result ) {
    std::string msg;
    msg.append("phdmesh::ParallelInputStream( " );
    if ( 0 == parallel_machine_rank( comm ) && NULL != file_name ) {
      msg.append( file_name );
    }
    else {
      msg.append( "<NULL>" );
    }
    msg.append( " ) FAILED" );
    throw std::runtime_error(msg);
  }
}

void ParInBuf::close()
{
  if ( NULL != m_root_fp ) { std::fclose( m_root_fp ); m_root_fp = NULL ; }
  setg(NULL,NULL,NULL);
}

ParInBuf::~ParInBuf()
{ close(); }

int ParInBuf::underflow()
{
  char * const buf = m_buffer + BUFFER_PUTBACK ;
  int nread = 0 ;

  if ( gptr() == NULL || egptr() <= gptr() ) {
    if ( NULL != m_root_fp ) { nread = std::fread(buf,1,MAX_READ,m_root_fp); }
    broadcast( m_comm , & nread , sizeof(int) );
  }

  if ( 0 < nread ) {
    broadcast( m_comm , buf , nread );
    setg( m_buffer , buf , buf + nread );
  }
  else {
    close();
  }

  return 0 < nread ? *buf : EOF ;
}

namespace {

void throw_overflow()
{
  std::string msg ;
  msg.append("phdmesh::ParallelInputStream::overflow CALL IS ERRONEOUS" );
  throw std::runtime_error(msg);
}

}

int ParInBuf::overflow( int )
{ throw_overflow(); return EOF ; }

int ParInBuf::sync()
{ return 0 ; }

std::streambuf * ParInBuf::setbuf( char * , std::streamsize )
{
  return this ;
}

//----------------------------------------------------------------------

} // namespace


ParallelInputStream::ParallelInputStream(
  ParallelMachine comm ,
  const char * const file_name )
  : std::istream( new ParInBuf( comm , file_name ) )
{}

ParallelInputStream::~ParallelInputStream()
{ delete rdbuf(); }

} // namespace phdmesh


