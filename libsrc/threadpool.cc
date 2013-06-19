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
Threadpool::Threadpool( unsigned int numThreads, int verbose )
  : maxThreads( boost::thread::hardware_concurrency() ),
    work( io_service ),
    verbose( verbose ) {
  
  // Store actual number of threads to be used.
  // If not specified, we use the system limit.
  threadCount = ( !numThreads ) ? numThreads : maxThreads;

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
    threads.create_thread( boost::bind( &boost::asio::io_service::run, &io_service ) );
  }

  // Declare threadpool operational
  if( verbose ) std::cout << "Threadpool is now operational." << std::endl;

}

/**
 * ~Threadpool()
 *
 * Desctructor. Closes io_service to prevent further
 * jobs added to job queue. Joins all threads in 
 * threadpool; currently running jobs are completed
 * before returning control to calling thread.
 */
Threadpool::~Threadpool() {
  io_service.stop();
  threads.join_all();
}

/**
 * post()
 *
 * Add a job to the job queue so that the next
 * idle thread in the threadpool may execute it.
 */
void Threadpool::post( Task task ) {
  // Add new task to job queue
  io_service.post( boost::bind( task ) );
}

/**
 * getThreadCount()
 *
 * Returns number of threads initiated 
 * in threadpool.
 */
unsigned int Threadpool::getThreadCount() {
  return threadCount;
}

/**
 * getMaxThreads()
 * 
 * Returns maximum number of threads available.
 * Provides easier access.
 */
unsigned int Threadpool::getMaxThreads() {
  return maxThreads;
}
