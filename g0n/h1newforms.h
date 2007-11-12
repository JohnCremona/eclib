// FILE H1NEWFORMS.H: redundant

#include "xsplit.h"

class h1newforms;
class jumps;

class h1newform {
public:
  vec bplus, bminus; // DUAL eigenvectors
  vector<long> aplist, aqlist; 
  scalar type;            // 2 for rectangular, 1 for triangular 
  long index;          // splitting index
//The following for non-square N only:
  long lplus, lminus;  // primes
  long mplus, mminus;  // magic factors
//The following for all N:
  long a,b,c,d,dotplus,dotminus;
  long pdot,dp0;       // dp0=1+p0-ap0, pdot = maninvector(p0).bplus
  rational loverp;    // =L(f,1)/per* = |pdot|/2*dp0
  long ap0;            // Eigenvalue of first "good" p
  long sfe;            // sign of F.E.
  long degphi;         // degree of Weil parametrization
  h1newform(void) {;}
  h1newform(const h1newforms* h, const homspace* h1, const jumps* jumpinfo,
            const vec vplus, const vec vminus, 
            const vector<long>& ap, long ind =0);
  h1newform(const vector<long>& ap,  ifstream& data);
  /*
  h1newform(const h1newform& f) 
    :bplus(f.bplus),bminus(f.bminus),aplist(f.aplist),aqlist(f.aqlist),
     type(f.type),  index(f.index),
  lplus(f.lplus),lminus(f.lminus),mplus(f.mplus),mminus(f.mminus),
  a(f.a),b(f.b),c(f.c),d(f.d),dotplus(f.dotplus),dotminus(f.dotminus),
  pdot(f.pdot),dp0(f.dp0),loverp(f.loverp),sfe(f.sfe),ap0(f.ap0),degphi(f.degphi)
 {;}
  void operator=(const h1newform& f)
    {
      bplus=f.bplus;bminus=f.bminus;aplist=f.aplist;aqlist=f.aqlist;
      type=f.type; index=f.index;
      lplus=f.lplus;lminus=f.lminus;mplus=f.mplus;mminus=f.mminus;
      pdot=f.pdot;dp0=f.dp0;loverp=f.loverp;sfe=f.sfe;ap0=f.ap0;
      degphi=f.degphi;
      a=f.a;b=f.b;c=f.c;d=f.d;dotplus=f.dotplus;dotminus=f.dotminus;
    }
  */
  void display(int detail=0) const;
  void writedata(ofstream& data);
};

class h1newforms :public level, splitter_base  {
  mat opmat(int i, int d, int v=0) 
  {return h1->opmat(i,d, v);}
  mat opmat_restricted(int i, const subspace& s, int d, int v=0) 
  {return h1->opmat_restricted(i,s,d,v);}
  smat s_opmat(int i, int d, int v=0) 
  {return h1->s_opmat(i,d, v);}
  smat s_opmat_restricted(int i, const ssubspace& s, int d, int v=0) 
  {return h1->s_opmat_restricted(i,s,d,v);}
  long matdim(void) {return h1->h1dim();} 
  long matden(void) {return h1->h1denom();}
  vector<long> eigrange(int i) {return h1->eigrange(i);}
  void use(const vec& bplus, const vec& bminus, const vector<long> apl);
  long dimoldpart(const class vector<long> ) {return 0;}
protected:
  long nap,ntp,nwq; int cuspidal;
  eigdata *filedata;
  homspace* h1; jumps* jumpinfo; ofstream data;
  long index;  // gets passed to use()
  void createfromeigs(long n, int verbose=0, long firstform=0); 
  void createfromdata(long n, int verbose=0);
  void writedata(long firstform=0);
public:
  long n1ds, iform, jform;
  vector<h1newform> nflist;
  long p0;  vec mvp;
  h1newforms(long n, int cuspidalflag=0, 
	     int disp=0, int olddata=0, long firstform=0);
  ~h1newforms(void) {;}
  void display(int detail=0) const;
  Curve getcurve(long i, int method, bigfloat& rperiod, int verbose=0);
  Cperiods getperiods(long i, int method=-1, int verbose=0);
// latter two implemented in periods.cc
};

class jumps {
private:
  long ncusps;
  mat* jumpmats;  //One per cusp
public:
  jumps(homspace *h1);              //Constructor
  ~jumps();       //Destructor
  long degphi(const vec&, const vec&, long) const;  //Computes degphi for f
};

char* datafile(long d);    //returns filename for data at level d
