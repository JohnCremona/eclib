/**
 * logger.h
 *
 * Declarations for logger class. 
 * Verbosity levels are integers:
 *    -1 : Error. Always print
 *     0 : Print if verbose > 0
 *     1 : Print if verbose > 1
 *    ... etc ...
 */

#ifndef LOGGER_H
#define LOGGER_H

// Include headers
#include <iostream>
#include <iomanip>
#include <sstream>

#ifdef MULTITHREAD
#include <boost/thread/thread.hpp>

// Uncomment if multithreading debug messages required
//#define ECLIB_MULTITHREAD_DEBUG
#endif

// Logging configurations
#ifndef ECLOG1
#define ECLOG1(v) if(eclogger::level() <= (v)) ; else eclogger().stream()
#endif
#ifndef ECLOG2
#define ECLOG2(v) if(eclogger::level() <= (v)) ; else eclogger().stream(__FILE__,__LINE__)
#endif

// Expose logging macro
#ifndef ECLOG
#define ECLOG ECLOG1
#endif

class eclogger {
  public:
    eclogger();
    ~eclogger();

    std::ostringstream& stream();
    std::ostringstream& stream( const char *file, const unsigned long line );

    static int  level();
    static void setLevel( const int verbose );

  private:
    static int level_;
    std::ostringstream s;
};

#endif
