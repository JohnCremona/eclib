/**
 * xsplit_data.cc
 *
 * Implementation of data structure used in xsplit.cc (form_finder class).
 */

// Include headers
#include "eclib/logger.h"
#include "eclib/xsplit_data.h"
#include "eclib/xsplit.h"

/**
 * ff_data()
 *
 * Main constructor.
 * Simply initiates private variables.
 */
ff_data::ff_data( form_finder* ff )
  : ff_( ff ),
    status_( INTERNAL ),
    depth_( 0 ),
    subdim_( 0 ),
    eigenvalue_( 0 ),
    eiglist_(),
    nest_( NULL ),
    conjmat_( NULL ),
    the_opmat_( NULL ),
    parent_( NULL ) {}

/**
 * ~ff_data()
 *
 * Destructor.
 * Only delete the objects which this instance created. 
 * Prevents destruction of objects which may still be required
 * by others, especially when running concurrently.
 */
ff_data::~ff_data() {
  // Delete dynamically created objects
  delete nest_;
  nest_ = NULL;
}

#ifdef MULTITHREAD
/**
 * operator()
 *
 * Overloaded operator(). Required for object to be passed to 
 * job queue. Task is to scan test eignvalues to identify 
 * newforms, or  whether further recursion is necessary.
 */
void ff_data::operator()() {
  // Call go_down() on current node, passing through its parent
  // to keep interface consistent with original.
  ff_ -> go_down( *(this->parent_), eigenvalue_, 0 );

#ifdef ECLIB_MULTITHREAD_DEBUG
  ECLOG(1) << "Executing node (eig=" << eigenvalue_ << " depth=" <<depth_ << ")" << endl;
#endif

  // Call find() on current node
  if( subdim_ > 0 ) ff_ -> find( *this );

  // Call go_up() only if this branch has ended
  ECLOG(0) << "child->status=" << child->status_ << " child->subdim=" << child->subdim_ << endl;
  if( status_ != INTERNAL || subdim_ == 0 ) ff_ -> go_up( *this );

#ifdef ECLIB_MULTITHREAD_DEBUG
  ECLOG(1) << "Completed node (eig=" << eigenvalue_ << " depth=" <<depth_ << " status=" << status_ << ")" << endl;
#endif
}
#endif

/**
 * status()
 *
 * Return status of current node.
 */
nodestatus ff_data::status() {
  return status_;
}

/**
 * submats()
 *
 * Returns relevant subspace of current depth.
 */
ssubspace* ff_data::submats() {
  return nest_;
}

/**
 * depth()
 *
 * Returns current depth. 
 */
long ff_data::depth() {
  return depth_;
}

/**
 * subdim()
 *
 * Return subdimension
 */
long ff_data::subdim() {
  return subdim_;
}

/** 
 * eig()
 *
 * Returns eigenvalue corresponding to current instance of class.
 */
long ff_data::eig() {
  return eigenvalue_;
}

/**
 * eiglist()
 *
 * Build sequence of eigenvalues for current node by traversing up the tree.
 * If an eiglist has already been computed, we simply return it to 
 * avoid accessing more nodes than necessary.
 */
vector< long > ff_data::eiglist() {
  // Return precomputed eiglist if available
  if( !eiglist_.empty() ) return eiglist_;

  // Root node (depth==0). Return empty vector.
  if( parent_ == NULL ) return vector< long >();

  // Else, we concatenate lists from further up the tree.
  eiglist_ = parent_ -> eiglist();
  eiglist_.push_back( eigenvalue_ );

  return eiglist_;
}

/**
 * child()
 *
 * Returns pointer to child w.r.t. given eigenvalue.
 */
ff_data* ff_data::child( long eig ) {
  return children_[ map(eig) ];
}

/**
 * numCompleteChildren()
 *
 * Returns number of completed children for current node.
 */
int ff_data::numCompleteChildren() {
  int completeCount = 0;

  vector< childstatus >::iterator it;
  for( it = completedChildren_.begin(); it != completedChildren_.end(); it++ ) {
    if( *it != NOT_COMPLETE ) completeCount++; 
  }

  return completeCount;
}

/**
 * complete()
 *
 * Return true if all children complete.
 */
bool ff_data::complete() {
  return ( numCompleteChildren() == numChildren_ ) ? true : false;
}

/** 
 * setStatus()
 *
 * Store status of current node.
 */
void ff_data::setStatus( nodestatus s ) {
  status_ = s;
}

/**
 * increaseDepth()
 *
 * Increases current depth. First check delta is positive.
 */
void ff_data::increaseDepth( long delta ) {
  assert( delta > 0 );
  depth_ += delta;
}

/**
 * decreaseDepth()
 *
 * Decrease current depth. First check delta is positive.
 */
void ff_data::decreaseDepth( long delta ) {
  assert( delta > 0 );
  depth_ -= delta;
}

/**
 * increaseSubmatUsage()
 *
 * Locked counter increment method.
 */
void ff_data::increaseSubmatUsage() {
#ifdef MULTITHREAD
  boost::mutex::scoped_lock lock( submat_lock_ );
#endif

#ifdef ECLIB_MULTITHREAD_DEBUG
  ECLOG(2) << "Increasing submat usage from " << submatUsage_ << " to " << submatUsage_+1 
           << " for node eig=" << eigenvalue_ << " depth=" << depth_ << endl;
#endif

  ++submatUsage_;
}

/**
 * storeBplus()
 *
 * Copies given vector into object storage.
 * Use vec class overloaded operator=.
 */
void ff_data::storeBplus( vec bp ) {
  bplus_ = bp;
}

/**
 * storeBminus()
 *
 * Copies given vector into object storage.
 * Use vec class overloaded operator=.
 */
void ff_data::storeBminus( vec bm ) {
  bminus_ = bm;
}

/** 
 * addChild()
 *
 * Adds a new data node to the children vector.
 */
void ff_data::addChild( long eig, ff_data &child ) {
  child.setParent( this ); 
  child.setEigenvalue( eig );
  children_[map(eig)] = &child;
}

/**
 * eraseChild()
 * 
 * Calls the destructor for the data node corresponding to given eigenvalue.
 */
void ff_data::eraseChild( long eig ) {
  eraseChild( map(eig) );
}

/**
 * eraseChild()
 *
 * Overloaded method. Main method for destroying children.
 */
void ff_data::eraseChild( int idx ) {
#ifdef ECLIB_MULTITHREAD_DEBUG
  ECLOG(2) << "Deleting node (eig=" << children_[idx]->eigenvalue_ 
           << " depth=" << depth_+1 << " status=" << children_[idx]->status_ << ")" << endl;
#endif

  delete children_[ idx ];
  children_[ idx ] = NULL;
  completedChildren_[ idx ] = DESTROYED;
}

/**
 * setParent()
 *
 * Stores pointer to parent data node.
 */
void ff_data::setParent( ff_data *parent ) {
  parent_ = parent;
}

/**
 * setEigenvalue()
 *
 * Stores eigenvalue.
 */
void ff_data::setEigenvalue( long eig ) {
  eigenvalue_ = eig;
}

/**
 * numChildren()
 *
 * Store number of children and resize vectors to correct size.
 */
void ff_data::numChildren( int size ) {
  numChildren_ = size;

  children_.resize( size, NULL );
  completedChildren_.resize( size, NOT_COMPLETE );
}

/**
 * childStatus()
 *
 * Sets value in map to 'flag', given a child node.
 * Monitors how many of a nodes children have completed.
 */
void ff_data::childStatus( long eig, childstatus flag ) {
#ifdef MULTITHREAD
  boost::mutex::scoped_lock lock( childComplete_lock_ );
#endif

  completedChildren_[map(eig)] = flag;
}

/**
 * eraseCompletedChldren()
 *
 * Loops through containers and destroying completed children.
 */
void ff_data::eraseCompletedChildren() {
  for( int i = 0; i < numChildren_; i++ ) {
    if( completedChildren_[i] == COMPLETE ) eraseChild( i );
  }
}

/**
 * map()
 *
 * Hash function to map given eigenvalue
 * to an index value. Removes dependancy on unordered_map.
 */
int ff_data::map( long eig ) {
  if( numChildren_ == 2 ) {
    return ( eig == 1 ) ? eig : 0;
  } else {
    return eig + ( numChildren_ - 1 ) * 0.5;
  }
}

// end of XSPLIT_DATA.CC
