#include "TSFSystem.h"


#include <unistd.h>
#include <time.h>


using namespace TSF;


string TSFSystem::getcwd()
{
	char cwd[201];

	::getcwd(cwd, 200);

	string rtn(cwd);
	
	return rtn;
}
 
int TSFSystem::system(const string& cmd)
{

	int rtn = ::system(cmd.c_str());
	if (rtn != 0)
		{
			TSFError::raise("system call " + cmd + " failed");
		}
	return rtn;
}
	
	
int TSFSystem::chdir(const string& newDir)
{
	int rtn = ::chdir(newDir.c_str());
	if (rtn != 0) 
		{
			TSFError::raise("could not chdir to " + newDir);
		}
	return rtn;
}

string TSFSystem::getenv(const string& varName)
{
	char* rtn = ::getenv(varName.c_str());
	if (rtn==0) return string();
	else return string(rtn);
}

string TSFSystem::date()
{
	time_t t;
	t = time(0);
	struct tm* lt = localtime((const time_t*) &t);
	char str[101];
	strftime(str, 100, "%a %b %d %T %Z %Y", lt);
	return str;
}

void TSFSystem::nsleep(int sec, int nsec)
{
	struct timespec request ;
	struct timespec remaining ;

	request.tv_sec = (time_t) sec;
	request.tv_nsec = nsec;

	nanosleep(&request, &remaining);
}

void TSFSystem::usleep(int usec)
{
	int sec = usec/1000000L;
	int uleft = usec - 1000000L*sec;
	TSFSystem::nsleep(sec, uleft*1000);
}

void TSFSystem::msleep(int msec)
{
	int sec = msec/1000L;
	int mleft = msec - 1000L*sec;
	TSFSystem::nsleep(sec, mleft*1000000L);
}

void TSFSystem::sleep(int sec)
{
	TSFSystem::msleep(sec*1000);
}

