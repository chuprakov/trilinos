#ifndef XMLINTERFACE_H
#define XMLINTERFACE_H

#include "TSFConfig.h"

#include "Socket.h"
#include "XMLObject.h"
#include "Mutex.h"


namespace TSF
{

	/** \ingroup XML 
	 * Interface for using XML for socket-based communication of objects between
	 * processes.
	 */

	class XMLInterface
		{
		public:
			/** empty ctor */
			XMLInterface() : socket_() {;}
			/** create an XMLInterface on a specified socket */
			XMLInterface(const Socket& socket);

			/** Receive an XML-encapsulated object from the interface socket */
			XMLObject recvObject() const ;
			/** Send an XML-encapsulated object to the interface socket */
			void sendObject(const XMLObject& obj) ;

			/** nonblocking probe for incoming data */
			bool hasDataToRead();

			/** Turn on display of all messages */
			static void traceOn() {trace_ = true;}
			/** Turn off display of all messages */
			static void traceOff() {trace_ = false;}

			void sendEOT();

		private:
			Socket socket_;
			Mutex lock_;
			static bool trace_;
		};



}
#endif

