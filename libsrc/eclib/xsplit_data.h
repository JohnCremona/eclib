/**
 * xsplit_data.h
 *
 * Declarations for data structure used in xsplit.cc (form_finder class).
 */

#ifndef XSPLIT_DATA_H
#define XSPLIT_DATA_H

// Include headers
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <vector>

#ifdef MULTITHREAD 
#include <boost/thread/mutex.hpp>
#endif

// Header files for macros and custom data structures
#include "method.h"
#include "smatrix_elim.h"

// Forward declaration of classes ... prevents circular includes.
// form_finder only required as pointer
class form_finder;

// Ennumerate node and child flags, exposed to global scope for outside access
enum nodestatus  { INTERNAL, ALL_OLD, FOUND_NEW, COMMON_EIGENSPACE };
enum childstatus { NOT_COMPLETE, COMPLETE, DESTROYED };

/**
 * ff_data (form_finder_data)
 *
 * Class name follows convention of that contained in xsplit.cc for ease of use.
 */
class ff_data {
  public:
    // Constructor, destructor.
    ff_data( form_finder* ff );
    ~ff_data();

#ifdef MULTITHREAD
    void operator()();                              // Executed upon submission to job queue
#endif

    // Getters (to maintain consistency)
    nodestatus     status();                        // Return status of current node
    ssubspace*     submats();                       // Return parent subspace
    long           depth();                         // Return current depth
    long           subdim();                        // Return subdimension
    long           eig();                           // Return associated eigenvalue
    vector< long > eiglist();                       // Return sequence of eigenvalues
    ff_data*       child( long eig );               // Return pointer to child
    int            numCompleteChildren();           // Return number of complete children
    bool           complete();                      // Return true if all children complete
  
    // Modifiers (available if needed i.e. not a friend class)
    void setStatus( nodestatus s );                 // Store status of current node

    void increaseDepth( long delta = 1 );           // Increase current depth
    void decreaseDepth( long delta = 1 );           // Decrease current depth
    
    void increaseSubmatUsage();                     // Locked counter incrementer

    void storeBplus( vec bp );                      // Store vector bplus
    void storeBminus( vec bm );                     // Store vector bminus
    
    void addChild( long eig, ff_data &child );      // Store new child of current node 
    void eraseChild( long eig );                    // Destroy child given eigenvalue
    void eraseChild( int idx );                     // Destroy child given index
    
    void setParent( ff_data *parent );              // Store parent node pointer
    void setEigenvalue( long eig );                 // Store eigenvalue

    void numChildren( int size );                   // Store number of children
    void childStatus( long eig, childstatus flag ); // Sets completed children
    
    void eraseCompletedChildren();                  // Destroys completed children

    // Make form_finder class a friend to gain access to protected/private methods
    friend class form_finder;

  private:
    form_finder*        ff_;                 // Back-pointer to form_finder class
                                             // Allows access to form_finder methods      

    nodestatus          status_;             // Status of current node
    long                depth_;              // Indicator of current depth
    long                subdim_;             // Dimension of current subspace
    long                eigenvalue_;         // Corresponding eigenvalue
    vector< long >      eiglist_;            // Sequence of eigenvalues leading to current
    vec                 bplus_, bminus_;
    ssubspace*          nest_;               // Current subspace (dynamically created)
    smat                conjmat_;            // Used only when plus==0 and bigmats==1
    smat                the_opmat_;
    smat                submat_;

    ff_data*            parent_;             // Pointer to parent data node
    vector<ff_data*>    children_;           // Pointers to corresponding data nodes
    int                 numChildren_;        // Store number of children
    int                 childCount_;         // Class counter for children
    vector<childstatus> completedChildren_;  // Flags for child completion
    int                 submatUsage_;        // Counter for submat

#ifdef MULTITHREAD
    boost::mutex        childComplete_lock_; // Lock for completed children
    boost::mutex        submat_lock_;        // Lock for submat usage
    boost::mutex        go_up_lock_;         // Lock for go_up() method    
#endif

    // Helper methods
    int map( long eig );                    // Map eigenvalue to index
};

#endif

// end of XSPLIT_DATA.H
