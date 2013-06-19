/**
 * xsplit_data.h
 *
 * Declarations for data structure used
 * in xsplit.cc (form_finder class).
 */

#ifndef XSPLIT_DATA_H
#define XSPLIT_DATA_H

// Include headers
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>

#include "eclib/xsplit.h"
#include "eclib/threadpool.h"

// Forward declaration of classes
class form_finder;

/**
 * ff_data (form_finder_data)
 *
 * Class name follows convention of that
 * contained in xsplit.cc for ease of use.
 */
class ff_data {
  public:
    // Constructor, destructor.
    ff_data();
    ~ff_data();

    // Overload operators
    void operator()();                      // Executed upon submission to job queue

    // Getters
    ssubspace*     submats( long depth );   // Return parent subspace of current depth
    long           depth();                 // Return current depth
    long           subdim();                // Return subdimension
    long           eig();                   // Return associated eigenvalue
    vector< long > eiglist();               // Return sequence of eigenvalues

    // Modifiers (available if needed i.e. not a friend class)
    void increaseDepth( long delta = 1 );   // Increase current depth
    void decreaseDepth( long delta = 1 );   // Decrease current depth
    void setEiglist( long idx, long eig );  // Store eigenvalue in local eiglist
    void storeBplus( vec bp );              // Store vector bplus
    void storeBminus( vec bm );             // Store vector bminus
    void addChild( ff_data *child );        // Store new child of current node 

    // Make form_finder class a friend to gain access to protected/private methods
    friend form_finder;

  private:
    long               depth_;              // Indicator of current depth
    long               subdim_;             // Dimension of current subspace
    long               eigenvalue_;         // Corresponding eigenvalue
    vector< long >     eiglist_;            // Sequence of eigenvalues leading to current
    vec                bplus_, bminus_;     //
    ssubspace*         nest_;               // Pointer to current subspace (to be created)
    smat               conjmat_;            // Used only when plus==0 and bigmats==1
    smat               the_opmat_;          //
    smat*              submat_;             // Pointer to 

    ff_data*           parent_;             // Pointer to parent data node
    vector< ff_data* > children_;           // Pointers to corresponding data nodes

    // Multithreading
    boost::mutex child_lock_;               // Lock for adding a new child node
}

#endif

// End of XSPLIT_DATA.H
