// FILE MSPLITTER.H: (REDUNDANT) Declaration of class splitter 
//                               version using Integers
//
// (also dsplitter (dual splitter))

#include "marith.h"
#include "mvector.h"
#include "mmatrix.h"
#include "msubspace.h"

int startswith(const longlist& a, const longlist& b, int l, int from=0);

class splitter {
private:
  homspace* h1;
  int maxdepth, depth, use_ws;
  msubspace** nest;     // array of pointers to subspaces
  longlist aplist, plist;
  mvector basisvector, bplus, bminus;
  int *havemat;
  matrix conjmat;
  matrix *opmats;
  int plusflag, dual, method, verbose;
//method=0 for usual elimination using Integers
//method=3 works mod P throughout, lifting only at end.
public:
  splitter(homspace* h, int plus, int d, int dualflag, int method, int v=0, int uw=1);
  ~splitter(void); 
  void splitoff(const longlist& apl);
  vector getbasis() const {return basisvector.shorten();}
  vector getbasisplus() const {return bplus.shorten();}
  vector getbasisminus() const {return bminus.shorten();}
};

class dsplitter {
private:
  homspace* h1;
  int maxdepth, depth;
  msubspace **lnest, **rnest;
  longlist aplist, plist;
  vector lbasisvector, lbplus, lbminus, rbasisvector, rbplus, rbminus;
  int *havemat;
  matrix conjmat, dualconjmat;
  matrix *opmats, *dualmats;
  int plusflag, dualflag, method, verbose;
//method=0 for usual elimination using Integers
//method=3 works mod P throughout, lifting only at end.
public:
  dsplitter(homspace* h, int plus, int d, int dflag, int method, int v=0);
  ~dsplitter(void);
// {if(dualflag) {delete rnest; delete opmats;} 
//		    delete lnest; delete havemat; delete dualmats;
//		  }
  void splitoff(const longlist& apl);
  vector getbasis() const {return rbasisvector;}
  vector getbasisplus() const {return rbplus;}
  vector getbasisminus() const {return rbminus;}
  vector getdbasis() const {return lbasisvector;}
  vector getdbasisplus() const {return lbplus;}
  vector getdbasisminus() const {return lbminus;}
  long getindex();
};

