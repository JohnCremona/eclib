/**
 * xsplit_data.cc
 *
 * Implementation of data structure used
 * in xsplit.cc (form_finder class).
 */

// Include headers
#include "eclib/xsplit_data.h"

/**
 * ff_data()
 *
 * Main constructor.
 * Simply initiates private variables.
 */
ff_data::ff_data()
  : depth_( 0 ),
    subdim_( 0 ),
    eigenvalue_( 0 ),
    nest_( NULL ),
    conjmat_( NULL ),
    the_opmat_( NULL ),
    parent_( NULL ) {
}

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
  delete submat_

  // Delete children data nodes
  for( int i = 0; i < children_.size(); i++ ) {
    if( children_[i] != NULL ) delete children_[i];
  }
}

/**
 * operator()
 *
 * Overloaded operator(). Required for object to be passed to 
 * job queue. Task is to scan test eignvalues to identify 
 * newforms, or  whether further recursion is necessary.
 */
void ff_data::operator()() {
 // Call find() on current eigenvalue
 if( subdim_ > 0 ) ff_.find( this );
 //   In find(), create a new instance of ff_data for each test
 //   eigenvalue. Pass on depth, nest, and eiglist.
 
 // Upon completion, delete current instance of ff_data 
 // i.e. commit suicide
 // PROVIDED we do not reference anything of this instance 
 // anywhere in the remaining code, and otherwise.
 delete this;
}

/**
 * submats()
 *
 * Returns relevant subspace of current depth.
 * Performs checks for valid access.
 */
ssubspace* ff_data::submats( long depth ) {
  assert( depth == depth_ || depth == depth_ - 1  );
  return ( depth == depth_ ) ? nest_current : nest_parent;
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
 * Returns eigenvalue corresponding to
 * current instance of class
 */
long ff_data::eig() {
  return eig_;
}

/**
 * eiglist()
 *
 * Return sequence of eigenvalues
 */
vector< long > ff_data::eiglist() {
  return eiglist_;
}
/**
 * eiglist()
 *
 * Build sequence of eigenvalues for current 
 * node by traversing up the tree.
 * If an eiglist has already been computed,
 * we simply return it to avoid accessing
 * more nodes than necessary.
 */
vector< long > ff_data::eiglist() {
  // Return precomputed eiglist if available
  if( eiglist_ != NULL ) {
    return eiglist_;
  }

  // Root node
  if( parent_ == NULL ) {
    return;
  }

  // Else, we concatenate lists from further
  // up the tree.
  eiglist_ = parent_.eiglist();
  eiglist_.push_back( eigenvalue_ );

  return eiglist_;
}

/**
 * increaseDepth()
 *
 * Increases current depth.
 * First check delta is positive.
 */
void ff_data::increaseDepth( long delta ) {
  assert( delta > 0 );
  depth_ += delta;
}

/**
 * decreaseDepth()
 *
 * Decrease current depth.
 * First check delta is positive.
 */
void ff_data::decreaseDepth( long delta ) {
  assert( delta > 0 );
  depth_ -= delta;
}

/**
 * setEiglist()
 *
 * Store given eigenvalue in local
 * eiglist. 
 * Index must equal current depth.
 */
void ff_data::setEiglist( long idx, long eig ) {
  assert( idx == depth_ );
  eiglist_[depth] = eig;
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
 * Provides multithreaded support to adding a 
 * new data node to the children vector.
 */
void ff_data::addChild( ff_data *child ) {
  // Lock vector with scoped lock 
  boost::mutex::scoped_lock lock( child_lock_ );

  // Add to vector
  children_.push_back( child );
}

// End of XSPLIT_DATA.CC
