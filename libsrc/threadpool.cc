/**
 * threadpool.cc
 *
 * Simple, no-frills threadpool for C++11 threads or Boost threads.
 * Some features have been made available specifically for ECLIB.
 *
 * These include use of verbosity level in constructor, to be 
 * used for logging. Tasks must also have an overloaded
 * operator() which describes what each thread must perform on 
 * that object. The return value of the task must be
 * void, i.e. nothing is returned.
 */

// Compile this only if Boost is installed
#ifdef MULTITHREAD

// Include header files
#include "eclib/threadpool.h"

/**
 * Threadpool()
 *
 * Main constructor. Creates required number of
 * threads, and sets them to wait for jobs.
 * Additional counters are initialised and
 * cached for ease of access.
 */
//template< class Task >
threadpool::threadpool( unsigned int numThreads, int verbose )
  : maxThreads( boost::thread::hardware_concurrency() ),
    work_( new boost::asio::io_service::work( io_service_ ) ),
    verbose( verbose ) {
 
  // Store actual number of threads to be used.
  // If not specified, we use the system limit.
  threadCount = ( numThreads > 0 ) ? numThreads : maxThreads;

  // We limit the number of threads to the system limit.
  // Note: it is best to not use all available threads
  if( threadCount > maxThreads ) {
    // Reset limit
    threadCount = maxThreads;

    // Notify
    if( verbose ) std::cout << "Requested more threads than available." 
                            << std::endl; 
  }

  // Declare the final number of threads to be used
  if( verbose ) std::cout << "Threadpool will be using " << threadCount 
                          << " threads." << std::endl;

  // Create threads and add to threadpool
  for( unsigned int i = 0; i < threadCount; i++  ) {
    threads_.create_thread( boost::bind( &boost::asio::io_service::run, &io_service_ ) );
  }

}

/**
 * ~Threadpool()
 *
 * Desctructor. Simply calls close() 
 */
//template< class Task >
threadpool::~threadpool() {
  close();
}


/**
 * close()
 *
 * Closes io_service_ to prevent further
 * jobs added to job queue. Joins all threads in 
 * threadpool; currently running jobs are completed
 * before returning control to calling thread.
 *
 * This blocking method exists in case we wish 
 * to close the threadpool before end-of-scope,
 * or to detect when all previously posted 
 * jobs have finished.
 */
void threadpool::close() {
  // We destroy the work class on the io_service object
  // so that we can exit once all jobs posted have finished
  work_.reset();

  // run() blocks until all posted jobs have finished
  io_service_.run();

  // We close the threadpool and join all threads
  io_service_.stop();
  threads_.join_all();
}

/**
 * getThreadCount()
 *
 * Returns number of threads initiated 
 * in threadpool.
 */
//template< class Task >
unsigned int threadpool::getThreadCount() {
  return threadCount;
}

/**
 * getMaxThreads()
 * 
 * Returns maximum number of threads available.
 * Provides easier access.
 */
//template< class Task >
unsigned int threadpool::getMaxThreads() {
  return maxThreads;
}

#endif // MULTITHREAD
