#ifndef TSFSYSTEM_H
#define TSFSYSTEM_H

#include "TSFDefs.h"

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

      static void hold();
    };



}
#endif




