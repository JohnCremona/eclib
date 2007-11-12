// FILE SPLITTER.H: Redundant 

int startswith(const longlist& a, const longlist& b, long l, long from=0);

class splitter {
private:
  homspace* h1;
  long maxdepth, depth; int use_ws;
  SUBSP** nest;     // array of pointers to subspaces
  longlist aplist, plist;
  vector basisvector, bplus, bminus;
  int *havemat;
  matrix conjmat;
  matrix *opmats;
  int plusflag, dual, verbose;
public:
  splitter(homspace* h, int plus, long d, int dualflag, int v=0, int uw=1);
  ~splitter(void); 
  void splitoff(const longlist& apl);
  vector getbasis() const {return basisvector;}
  vector getbasisplus() const {return bplus;}
  vector getbasisminus() const {return bminus;}
};

class dsplitter {
private:
  homspace* h1;
  long maxdepth, depth;
  SUBSP **lnest, **rnest;
  longlist aplist, plist;
  vector lbasisvector, lbplus, lbminus, rbasisvector, rbplus, rbminus;
  int *havemat;
  matrix conjmat, dualconjmat;
  matrix *opmats, *dualmats;
  int plusflag, dualflag, verbose;
public:
  dsplitter(homspace* h, int plus, long d, int dflag, int v=0);
  ~dsplitter(void);
// {if(dualflag) {delete rnest; delete opmats;} 
//      	    delete lnest; delete havemat; delete dualmats;
//		  }
  void splitoff(const longlist& apl);
  vector getbasis() const {return rbasisvector;}
  vector getbasisplus() const {return rbplus;}
  vector getbasisminus() const {return rbminus;}
  vector getdbasis() const {return lbasisvector;}
  vector getdbasisplus() const {return lbplus;}
  vector getdbasisminus() const {return lbminus;}
  long getindex(int detail=0);
};

