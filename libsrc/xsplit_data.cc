// FILE xsplit_data.cc : Implementation of member functions for class ff_data
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2023 Marcus Mo
// 
// This file is part of the eclib package.
// 
// eclib is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// eclib is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with eclib; if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
// 
//////////////////////////////////////////////////////////////////////////

// Include headers
#include "eclib/logger.h"
#include "eclib/xsplit.h"

template class ff_data<int>;
template class ff_data<long>;
template class ff_data<ZZ>;

/**
 * ff_data()
 *
 * Main constructor.
 * Simply initiates private variables.
 */
template<class T>
ff_data<T>::ff_data( form_finderT<T>* ff )
  : ff_( ff ),
    status_( INTERNAL ),
    depth_( 0 ),
    subdim_( 0 ),
    eigenvalue_( 0 ),
    eigrange_( 0 ),
    eiglist_(),
    abs_space_( NULL ),
    rel_space_( NULL ),
    conjmat_(),
    the_opmat_(),
    parent_( NULL ),
    numChildren_( 0 ) {}

/**
 * ~ff_data()
 *
 * Destructor.
 * Only delete the objects which this instance created. 
 * Prevents destruction of objects which may still be required
 * by others, especially when running concurrently.
 */
template<class T>
ff_data<T>::~ff_data() {
  // Delete dynamically created objects
  delete abs_space_;
  delete rel_space_;
}

#ifdef ECLIB_MULTITHREAD
/**
 * operator()
 *
 * Overloaded operator(). Required for object to be passed to 
 * job queue. Task is to scan test eigenvalues to identify 
 * newforms, or  whether further recursion is necessary.
 */
template<class T>
void ff_data<T>::operator()() {
  // Call go_down() on current node, passing through its parent
  // to keep interface consistent with original.
  ff_ -> go_down( *(this->parent_), eigenvalue_, 0 );

#ifdef ECLIB_MULTITHREAD_DEBUG
  ECLOG(1) << "Executing node (eig=" << eigenvalue_ << " depth=" <<depth_ << ")" << endl;
#endif

  // Call find() on current node
  if( subdim_ > 0 ) ff_ -> find( *this );

  // Call go_up() only if this branch has ended
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
template<class T>
nodestatus ff_data<T>::status() {
  return status_;
}

/**
 * abs_space()
 *
 * Returns absolute subspace of current depth.
 */
template<class T>
ssubZspace<T>* ff_data<T>::abs_space() {
  return abs_space_;
}

/**
 * rel_space()
 *
 * Returns relative subspace of current depth.
 */
template<class T>
ssubZspace<T>* ff_data<T>::rel_space() {
  return rel_space_;
}

/**
 * depth()
 *
 * Returns current depth. 
 */
template<class T>
long ff_data<T>::depth() {
  return depth_;
}

/**
 * subdim()
 *
 * Return subdimension
 */
template<class T>
long ff_data<T>::subdim() {
  return subdim_;
}

/** 
 * eig()
 *
 * Returns eigenvalue corresponding to current instance of class.
 */
template<class T>
long ff_data<T>::eig() {
  return eigenvalue_;
}

/**
 * eiglist()
 *
 * Build sequence of eigenvalues for current node by traversing up the tree.
 * If an eiglist has already been computed, we simply return it to 
 * avoid accessing more nodes than necessary.
 */
template<class T>
vector< long > ff_data<T>::eiglist() {
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
template<class T>
ff_data<T>* ff_data<T>::child( long eig ) {
  return children_[ map(eig) ];
}

/**
 * numCompleteChildren()
 *
 * Returns number of completed children for current node.
 */
template<class T>
int ff_data<T>::numCompleteChildren() {
  return std::count_if(completedChildren_.begin(), completedChildren_.end(),
                       [](childstatus s) {return s!=NOT_COMPLETE;});
}

/**
 * complete()
 *
 * Return true if all children complete.
 */
template<class T>
bool ff_data<T>::complete() {
  return ( numCompleteChildren() == numChildren_ ) ? true : false;
}

/** 
 * setStatus()
 *
 * Store status of current node.
 */
template<class T>
void ff_data<T>::setStatus( nodestatus s ) {
  status_ = s;
}

/**
 * increaseDepth()
 *
 * Increases current depth. First check delta is positive.
 */
template<class T>
void ff_data<T>::increaseDepth( long delta ) {
  assert( delta > 0 );
  depth_ += delta;
}

/**
 * decreaseDepth()
 *
 * Decrease current depth. First check delta is positive.
 */
template<class T>
void ff_data<T>::decreaseDepth( long delta ) {
  assert( delta > 0 );
  depth_ -= delta;
}

/**
 * increaseSubmatUsage()
 *
 * Locked counter increment method.
 */
template<class T>
void ff_data<T>::increaseSubmatUsage() {
#ifdef ECLIB_MULTITHREAD
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
template<class T>
void ff_data<T>::storeBplus( Zvec<T> bp ) {
  bplus_ = bp;
}

/**
 * storeBminus()
 *
 * Copies given vector into object storage.
 * Use vec class overloaded operator=.
 */
template<class T>
void ff_data<T>::storeBminus( Zvec<T> bm ) {
  bminus_ = bm;
}

/** 
 * addChild()
 *
 * Adds a new data node to the children vector.
 */
template<class T>
void ff_data<T>::addChild( long eig, ff_data<T> &child ) {
  child.setParent( this ); 
  child.setEigenvalue( eig );
  children_[map(eig)] = &child;
}

/**
 * eraseChild()
 * 
 * Calls the destructor for the data node corresponding to given eigenvalue.
 */
template<class T>
void ff_data<T>::eraseChild( long eig ) {
  eraseChild( map(eig) );
}

/**
 * eraseChild()
 *
 * Overloaded method. Main method for destroying children.
 */
template<class T>
void ff_data<T>::eraseChild( int idx ) {
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
template<class T>
void ff_data<T>::setParent( ff_data<T> *parent ) {
  parent_ = parent;
}

/**
 * setEigenvalue()
 *
 * Stores eigenvalue.
 */
template<class T>
void ff_data<T>::setEigenvalue( long eig ) {
  eigenvalue_ = eig;
}

/**
 * numChildren()
 *
 * Stores number of children and eigrange, and resize vectors to correct size.
 */
template<class T>
void ff_data<T>::setChildren( vector<long> eigs ) {
  numChildren_ = eigs.size();
  eigrange_ = eigs;

  children_.resize( numChildren_, NULL );
  completedChildren_.resize( numChildren_, NOT_COMPLETE );
}

/**
 * childStatus()
 *
 * Sets value in map to 'flag', given a child node.
 * Monitors how many of a nodes children have completed.
 */
template<class T>
void ff_data<T>::childStatus( long eig, childstatus flag ) {
#ifdef ECLIB_MULTITHREAD
  boost::mutex::scoped_lock lock( childComplete_lock_ );
#endif

  completedChildren_[map(eig)] = flag;
}

/**
 * eraseChldren()
 *
 * Loops through containers and destroys all children.
 */
template<class T>
void ff_data<T>::eraseChildren( ) {
  if( numChildren_ > 0 ) {
    for( int i = 0; i < numChildren_; i++ ) {
      if ( children_[i] != NULL) {
        children_[i] -> eraseChildren();
        eraseChild( i );
      }
    }
  }
}

/**
 * map()
 *
 * Hash function to map given eigenvalue
 * to an index value. Removes dependancy on unordered_map.
 *
 * N.B. This function no longer makes assumptions on the eigenvalues,
 * the number of which is numChildren_: in current practice, if
 * numChildren_==2, they are [-1,+1] while otherwise
 * numChildren_==2n+1 and they are [-n,...,-2,-1,0,1,2,...,n], but
 * this specific choice is no longer relied on.
 */
template<class T>
int ff_data<T>::map( long eig ) {
  int i = (int)(find(eigrange_.begin(),eigrange_.end(),eig)-eigrange_.begin());
  return i;
}

// end of XSPLIT_DATA.CC
