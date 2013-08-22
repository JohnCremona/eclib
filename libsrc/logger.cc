/**
 * logger.cc
 *
 * Simple thread-safe logger class. Prints to stdout.
 */

#include "eclib/logger.h"

int eclogger::level_ = 0;

/**
 * eclogger()
 *
 * Main constructor. Defaults verbosity to 0.
 */
eclogger::eclogger() {}

/**
 * ~eclogger()
 *
 * Destructor. Appends newline character.
 * Uses thread-safe fprintf, and flushes.
 */
eclogger::~eclogger() {
  //s << std::endl;
  fprintf( stdout, "%s", s.str().c_str() );  
  fflush( stdout );
}

/**
 * level()
 *
 * Returns logger verbosity level.
 */
int eclogger::level() {
  return level_;
}

/**
 * setLevel()
 *
 * Allows for global change of verbosity level.
 */
void eclogger::setLevel( const int verbose ) {
  level_ = verbose;
}

/**
 * stream()
 *
 * Returns stream object.
 */
std::ostringstream& eclogger::stream() {
#if defined(ECLIB_MULTITHREAD) && defined(ECLIB_MULTITHREAD_DEBUG)
  s << "Thread " << boost::this_thread::get_id() << "\t"; 
#endif

  return s;
}

/**
 * stream()
 *
 * Returns stream object with FILE and LINE details.
 */
std::ostringstream& eclogger::stream( const char *file, const unsigned long line ) {
  std::string filename = file;
  s << std::setw(20) << filename << std::setw(5) << line << " ";
  return stream();  
}
