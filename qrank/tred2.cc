// TRED2.CC:   program for reduction of quartics to minimal I,J
//    
#include "marith.h"
#include "unimod.h"
#include "points.h"
#include "mquartic.h"
#include "transform.h"
#include "msoluble.h"
#include "minim.h"
#include "reduce.h"

int getquartic(quartic& g);
void out(const quartic& g)
{
  cout<<g.geta()<<" "<<g.getb()<<" "<<g.getcc()<<" "<<g.getd()<<" "<<g.gete();
}

int main()
{
#ifdef LiDIA
  long lidia_precision=50;
//  cout<<"Enter number of decimal places: "; cin>>lidia_precision;
  bigfloat::precision(lidia_precision);
#else
  cout.precision(15);
#endif

  cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)
  
  int verb=0; 
  //   cout << "Verbose? "; cin >> verb;
  initprimes("PRIMES",verb);
  double hlim=8;
  //  cout << "Limit on height? "; cin >> hlim;

  quartic g;

  getquartic(g);
    {
//      out(g);
      bigint I = g.getI(), J=g.getJ();
//      cout<<"\t"<<I<<"\t"<<J<<endl;
      bigint ga=g.geta(), gb=g.getb(), gc=g.getcc(), gd=g.getd(), ge=g.gete();
      bigint p, badp;
      bigintArray plist = pdivs(g.getdisc());
      scaled_unimod m;

      bigint newa(ga), newb(gb), newc(gc), newd(gd), newe(ge);
      minim_all(newa,newb,newc,newd,newe,I,J,plist,m,0,0);
      quartic newg(newa,newb,newc,newd,newe);
      plist = pdivs(newg.getdisc());

      int locsol = locallysoluble(newg, plist, badp);
      if(!locsol)
	{
	  ;
	}
      minim_all(newa,newb,newc,newd,newe,I,J,plist,m,1,0);
      newg.assign(newa,newb,newc,newd,newe);
      plist = pdivs(newg.getdisc());
      if(check_transform(ga,gb,gc,gd,ge,m,newa,newb,newc,newd,newe))
	{
	  unimod m1;
	  reduce(newa,newb,newc,newd,newe,m1);
	  newg.assign(newa,newb,newc,newd,newe,newg.getroots(),0,I, J, 4*pow(I,3)-J*J);
	  out(newg);
	  cout<<"\t"<<I<<"\t"<<J<<endl;
	  m *= m1;
	  if(check_transform(ga,gb,gc,gd,ge,m,newa,newb,newc,newd,newe))
	    {
//	      cout << "OK\n";
	    }
	  else
	    {
	      cout << "check_transform fails after reduction!\n";
	    }
	}
      else
	{
	  cout << "check_transform fails after minimalization!\n";
	}
    }
}

int getquartic(quartic& g)
{
  bigint a, b, c, d, e;
  
//  cout << "Enter quartic coefficients a,b,c,d,e ?" << endl;
  char ch; cin>>ch;
  if((ch=='(')||(ch=='[')) cin>>a>>ch>>b>>ch>>c>>ch>>d>>ch>>e>>ch;
     else 
     {
       cin.putback(ch);
       cin >> a >> b >> c >> d >> e;
     }
     
  if (sign(a)==0&&sign(b)==0&&sign(c)==0&&sign(d)==0&&sign(e)==0)
     return 0;
    
  g=quartic(a,b,c,d,e);  // will set its own invariants, roots and type

  return 1;
}
