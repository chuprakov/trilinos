#ifndef TSFSYSTEM_H
#define TSFSYSTEM_H

#include "TSFDefs.h"

#include <cstdlib>
#include "TSFArray.h"

namespace TSF
{

  /**
   * \ingroup General
   * Operating system calls
   *
   * @author Kevin Long
   */

  class TSFSystem
    {
    public:
      /**
       * get current working directory
       */
      static string getcwd();

      /**
       * Change working directory
       */
      static int chdir(const string& newDir);

      /**
       * get current date
       */
      static string date();

      /**
       * look up an environment variable
       */
      static string getenv(const string& varName);

      /**
       * make a system call to execute a command
       */
      static int system(const string& cmd);

      /**
       * sleep for a specified number of seconds
       */
      static void sleep(int sec);

      /**
       * sleep for a specified number of milliseconds
       */
      static void msleep(int millis);

      /**
       * sleep for a specified number of microseconds
       */
      static void usleep(int usec);

      /**
       * sleep for a specified number of nanoseconds
       */
      static void nsleep(int sec, int nanos);

      static void hold();
    };



}
#endif




