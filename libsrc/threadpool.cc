/**
 * threadpool.cc
 *
 * Simple, no-frills threadpool for Boost threads.
 * Some features have been made available specifically for ECLIB.
 *
 * These include use of verbosity level in constructor, to be 
 * used for logging. Tasks must also have an overloaded
 * operator() which describes what each thread must perform on 
 * that object. The return value of the task must be
 * void, i.e. nothing is returned.
 */

// Compile this only if Boost is installed
#ifdef ECLIB_MULTITHREAD

// Include header files
#include "eclib/threadpool.h"

/**
 * Threadpool()
 *
 * Default constructor. Must call start() before using threadpool.
 */
threadpool::threadpool() 
  : maxThreads_( 0 ), threadCount_( 0 ), verbose_( -1 ),
    work_( new boost::asio::io_service::work( io_service_ ) )
{}

/**
 * threadpool()
 *
 * Main constructor.
 */
threadpool::threadpool( unsigned int numThreads, int verbose ) 
  : work_( new boost::asio::io_service::work( io_service_ ) ) {
  start( numThreads, verbose );
}

/**
 * ~Threadpool()
 *
 * Desctructor. Simply calls close() 
 */
threadpool::~threadpool() {
  close();
}

/**
 * start();
 *
 * Must be called after constructor, and allows for threadpool 
 * to be restarted after a call to close().
 */
void threadpool::start( unsigned int numThreads, int verbose ) {
  // Store verbosity
  verbose_ = verbose;

  // Store maximum number of threads system can support
  maxThreads_ = boost::thread::hardware_concurrency();

  // Store actual number of threads to be used.
  // If not specified, we use the system limit.
  threadCount_ = ( numThreads > 0 ) ? numThreads : maxThreads_;

  // We limit the number of threads to the system limit.
  // Note: it is best to not use all available threads
  if( threadCount_ > maxThreads_ ) {
    // Reset limit
    threadCount_ = maxThreads_;

    // Notify
    if( verbose_ ) std::cout << "Requested more threads than available." 
                            << std::endl; 
  }

  // Declare the final number of threads to be used
  if( verbose_ ) std::cout << "Threadpool will be using " << threadCount_ 
                           << " threads from a total of " 
                           << maxThreads_ << " threads." << std::endl;

  // Create threads and add to threadpool
  for( unsigned int i = 0; i < threadCount_-1; i++  ) {
    threads_.create_thread( boost::bind( &boost::asio::io_service::run, &io_service_ ) );
  }
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
unsigned int threadpool::getThreadCount() {
  return threadCount_;
}

/**
 * getMaxThreads()
 * 
 * Returns maximum number of threads available.
 * Provides easier access.
 */
unsigned int threadpool::getMaxThreads() {
  return maxThreads_;
}

#endif // ECLIB_MULTITHREAD
