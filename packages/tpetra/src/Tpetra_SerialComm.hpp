// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#ifndef TPETRA_SERIALCOMM_HPP
#define TPETRA_SERIALCOMM_HPP

#include <Teuchos_OrdinalTraits.hpp>
#include "Tpetra_Comm.hpp"
#include "Tpetra_Distributor.hpp"
#include <Teuchos_Object.hpp>

namespace Tpetra {

	//! Tpetra::SerialComm:  The Tpetra Serial Communication Class.
	/*! The SerialComm class is an implementation of Tpetra::Comm, providing the general
	    information and services needed for other Tpetra classes to run on a serial computer.
	*/

	template<typename OrdinalType, typename PacketType>
	class SerialComm : public Teuchos::Object, public virtual Comm<OrdinalType, PacketType> {
	public:

		//@{ \name Constructor/Destructor Methods

		//! Constructor
		/*! Builds an instance of a serial communicator.  Even
		    if the application is running in parallel via MPI, this communicator
			will execute in serial.  The access functions return the number of
			memory images to be 1 and the image ID to be 0.
		*/
		SerialComm() : Teuchos::Object("Tpetra::SerialComm") {}
  
		//! Copy constructor
		/*! Makes an exact copy of an existing SerialComm instance.
		 */
		SerialComm(SerialComm<OrdinalType, PacketType> const& comm) : Teuchos::Object(comm.label()) {}
  
		//! Destructor.
		/*! Completely deletes a SerialComm object.  
		  \warning Note:  All objects that depend on a SerialComm instance 
		           should be destroyed prior to calling this function.
		*/
		~SerialComm() {}

		//@}

		//@{ \name Barrier Methods

		//! Barrier function. 
		/*! A no-op for a serial communicator.
		 */
		void barrier() const {}

		//@}

		//@{ \name Broadcast Methods

		//! SerialComm Broadcast function.
		/*! A no-op for a serial communicator.
		  \param myVals InOut
		         On entry, the root image contains the list of values.  On exit,
				 all images will have the same list of values.  Note that values must be
				 allocated on all images before the broadcast.
		  \param count In
		         On entry, contains the length of myVals.
		  \param root In
		         On entry, contains the imageID from which all images will receive a copy of myVals.
		*/
		void broadcast(PacketType* myVals, OrdinalType const count, int const root) const {}

		//@}

		//@{ \name Gather Methods

		//! SerialComm All Gather function.
		/*! A copy for a serial communicator.
		  \param myVals In
		         On entry, contains the list of values, to be sent to all images.
		  \param allVals Out
		         On exit, contains the list of values from all images. Must by of size numImages*count.
		  \param count In
		         On entry, contains the length of myVals.
		*/
		void gatherAll(PacketType* myVals, PacketType* allVals, OrdinalType const count) const {
			copy(myVals, allVals, count);
		}

		//@}

		//@{ \name Sum Methods

		//! SerialComm Global Sum function.
		/*! A copy for a serial communicator.
		  \param partialSums In
		         On entry, contains the list of values, usually partial sums computed locally,
				 to be summed across all images.
		  \param globalSums Out
		         On exit, contains the list of values summed across all images.
		  \param count In
		         On entry, contains the length of partialSums.
		*/
		void sumAll(PacketType* partialSums, PacketType* globalSums, OrdinalType const count) const {
			copy(partialSums, globalSums, count);
		}

		//! Scattered Global Sum function
		/*! Take a list of input values from each image, compute the global sums, and distribute 
		    the list of sums across all images.
		  \param sendVals In
		         On entry, contains the list of values to sum from this image.
		  \param recvVals Out
		         On exit, contains the list of sums distributed to this image.
		  \param recvCounts In
		         On entry, contains a list of sizes. On exit, recvVals on image i will contain
				 recvCounts[i] entries.
		*/
		void sumAllAndScatter(PacketType* sendVals, PacketType* recvVals, OrdinalType* recvCounts) const {
			// In serial, we're on Image 0, so we should copy recvCounts[0] elements.
			copy(sendVals, recvVals, recvCounts[0]);
		}

		//@}
	
		//@{ \name Max/Min Methods

		//! SerialComm Global Max function.
		/*! A copy for a serial communicator.
		  \param partialMaxs In
		         On entry, contains the list of values, usually partial maxs computed locally;
				 using these Partial Maxs, the max across all images will be computed.
		  \param globalMaxs Out
		         On exit, contains the list of maxs computed across all images.
		  \param count In
		         On entry, contains the length of partialMaxs.
		*/
		void maxAll(PacketType* partialMaxs, PacketType* globalMaxs, OrdinalType const count) const {
			copy(partialMaxs, globalMaxs, count);
		}

		//! SerialComm Global Min function.
		/*! A copy for a serial communicator.
		  \param partialMins In
		         On entry, contains the list of values, usually partial mins computed locally;
				 using these Partial Mins, the min across all images will be computed.
		  \param globalMins Out
		         On exit, contains the list of mins computed across all images.
		  \param count In
		         On entry, contains the length of partialMins.
		*/
		void minAll(PacketType* partialMins, PacketType* globalMins, OrdinalType const count) const {
			copy(partialMins, globalMins, count);
		}

		//@}

		//@{ \name Parallel Prefix Methods

		//! SerialComm Scan Sum function.
		/*! A copy for a serial communicator.
		  \param myVals In
		         On entry, contains the list of values to be summed across all images.
		  \param scanSums Out
		         On exit, contains the list of values summed across images 0 through i.
		  \param count In
		         On entry, contains the length of myVals.
		*/
		void scanSum(PacketType* myVals, PacketType* scanSums, OrdinalType const count) const {
			copy(myVals, scanSums, count);
		}

		//@}

		//@{ \name Point-to-Point Methods

		//! SerialComm Blocking Send function
		/*! Sends a list of values to another image. Control will not exit this function until
		    the data have been acknowledged by the receiving image.
		  \param myVals In
		         On entry, contains the list of values to be sent.
		  \param count In
		         On entry, contains the length of myVals.
		  \param destImageID In
		         On entry, contains the ImageID of the image to send the values to.
		*/
		void send(PacketType* myVals, OrdinalType const count, int destImageID) const {
			throw reportError("This method should never be called.", -1);
		}

		//! SerialComm Blocking Receive function
		/*! Receive a list of values from another image. Control will not exit this function until
		    data has been received from the source image specified.
		  \param myVals Out
			     On exit, contains the list of values received.
		  \param count In
		         On entry, contains the length (capacity) of myVals.
		  \param sourceImageID In
		         On entry, contains the ImageID of the image to receive the values from.
				 (A sourceImageID of -1 means receive values from any image.)
		  \return If successful, the ImageID of the sender (>= 0).
		          If not successful, aminteger return code (< 0).
		*/
		int receive(PacketType* myVals, OrdinalType const count, int sourceImageID) const {
			throw reportError("This method should never be called.", -1);
		}

		//@}

		//@{ \name Execute Distributor Plan Methods

		//! doPostsAndWaits
		/*! Execute a plan specified by the distributor object passed in. 
		  \param distributor In
			     Contains the specifications of the plan we're executing.
		  \param exports In
		         On entry, contains the values we're exporting.
		  \param packetSize In
		         On entry, the number of PacketType variables that make up an element.
		  \param imports Out
		         On exit, contains the values exported to us. (imports will be resized
				 if necessary, and any existing values will be overwritten.)
		*/
		void doPostsAndWaits(Distributor<OrdinalType> const& distributor,
							 std::vector<PacketType>& exports,
							 OrdinalType packetSize,
							 std::vector<PacketType>& imports) {
			throw reportError("This method should never be called.", -1);
		}

		//! doPosts
		void doPosts(Distributor<OrdinalType> const& distributor,
					 std::vector<PacketType>& exports,
					 OrdinalType packetSize,
					 std::vector<PacketType>& imports) {
			throw reportError("This method should never be called.", -1);
		}

		//! doWaits
		void doWaits(Distributor<OrdinalType> const& distributor) {
			throw reportError("This method should never be called.", -1);
		}

		//! doReversePostsAndWaits
		/*! Execute a reverse plan specified by the distributor object passed in. 
		  \param distributor In
			     Contains the specifications of the plan we're reverse-executing.
		  \param exports In
		         On entry, contains the values we're exporting.
		  \param packetSize In
		         On entry, the number of PacketType variables that make up an element.
		  \param imports Out
		         On exit, contains the values exported to us. (imports will be resized
				 if necessary, and any existing values will be overwritten.)
		*/
		void doReversePostsAndWaits(Distributor<OrdinalType> const& distributor,
									std::vector<PacketType>& exports,
									OrdinalType packetSize,
									std::vector<PacketType>& imports) 
		{
			throw reportError("This method should never be called.", -1);
		}

		//! doReversePosts
		void doReversePosts(Distributor<OrdinalType> const& distributor,
							std::vector<PacketType>& exports,
							OrdinalType packetSize,
							std::vector<PacketType>& imports)
		{
			throw reportError("This method should never be called.", -1);
		}

		//! doReverseWaits
		void doReverseWaits(Distributor<OrdinalType> const& distributor)
		{
			throw reportError("This method should never be called.", -1);
		}

		//@}

		//@{ \name Image Info Methods

		//! getMyImageID - In serial mode, always returns 0.
		int getMyImageID() const {return(0);}

		//! getNumImages - In serial mode, always returns 1.
		int getNumImages() const {return(1);}

		//@}

		//@{ \name I/O Methods

		//! Print methods
		void print(ostream& os) const {}
		void printInfo(ostream& os) const {os << *this;}

		//@}

	private:
		// convenience function for copying
		void copy(PacketType* source, PacketType* dest, OrdinalType const count) const {
			for(OrdinalType i = Teuchos::OrdinalTraits<OrdinalType>::zero(); i < count; i++)
				dest[i] = source[i];
		}

	}; // class SerialComm

} // namespace Tpetra

#endif // TPETRA_SERIALCOMM_HPP
