#include "XMLInterface.h"
#include "InterruptibleInputSource.h"
#include "InterruptibleInputStream.h"
#include "Synch.h"
#include "TSFOut.h"



using namespace TSF;


bool XMLInterface::trace_ = false;

XMLInterface::XMLInterface(const Socket& socket)
	: socket_(socket), lock_()
{}

XMLObject XMLInterface::recvObject() const 
{
	Synch synch(lock_);
	XMLObject rtn;

	try
		{
			InterruptibleInputSource source(socket_);
			rtn = source.getObject();
			if (trace_) TSFOut::println("XMLInterface recvd: " + rtn.toString());
			//			source.resume();
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in XMLInterface::getObject()");
		}
	return rtn;
}

void XMLInterface::sendObject(const XMLObject& obj) 
{
	Synch synch(lock_);
	socket_.write(obj.toString());
	if (trace_) TSFOut::println("XMLInterface send: " + obj.toString());
	sendEOT();
}

bool XMLInterface::hasDataToRead()
{
	Synch synch(lock_);

	return socket_.hasDataToRead();
}

void XMLInterface::sendEOT()
{
	socket_.write("\f");
}





