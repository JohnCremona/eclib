/* SANSONE.CC:  MWRANK for Rene Schoof's Cassels-Sansone curves  */
/*                                */
#include "points.h"  // from qcurves library
#include "mwprocs.h" //  ""     "       "
#include "mquartic.h"
#include "mrank1.h"
#include "mrank2.h"

int getcurve(Curvedata& CD, Curvedata& CD_orig, 
	     bigint& u, bigint& r, bigint& s, bigint& t, 
	     int& change, int verb);

void post_process(const Point& p); // on [a,0,1,0,0]

int main()
{
#ifdef LiDIA
  long lidia_precision=40;
  cout<<"Enter number of decimal places: "; cin>>lidia_precision;
  bigfloat::precision(lidia_precision);
#else
  cout.precision(15);
#endif
  initprimes("PRIMES",0);
  int verbose=1, traceabc=0, traceequiv=0, ptl=0, success;
  long i, sha2, sha2dash, rank,  naux=5;
// Increase the following if you need to search further for rational
// points on quartics (logarithmic!)
  long hlimq=8, hlimc=8; 
  cout << "Verbose mode? (0/1)\n"; cin >> verbose;
  cout << "Bound on height in quartic point search (e.g. 5)? "; cin >> hlimq;
  cout <<"\n";

  cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)
  
  Curvedata CD, CD_orig;
  bigint u,r,s,t;             // transform CD_orig -> CD
  int change;                 // flags whether CD != CD_orig
  PointArray plist;

  while ( getcurve(CD,CD_orig,u,r,s,t,change,verbose)	 )
    {
      rank2 r2(&CD, verbose, 20, hlimq);
 
      if (r2.ok())
        {   
          cout << "Rank = " << r2.getrank() << "\n";
	  plist = r2.getpoints();
        }
      else
        {
          rank1 r1(&CD, verbose, traceabc, traceequiv, 20, hlimq);
 
          if (r1.ok())
            {    
              cout << "Rank = " << r1.getrank() << "\n";
	      plist = r1.getpoints();
            }
          else cout << "Failed to compute rank\n";
        }

      bigfloat ht, maxht = 0;
      for (i=0; i<plist.getlength(); i++)
	{ 
	  ht  = height(plist[i]);
	  if(ht>maxht) maxht=ht;
	}

      mw mwbasis(&CD, verbose);
      mwbasis.process(plist);
      rank = mwbasis.getrank();
      if(verbose) cout<<"Rank of points found is "<<rank<<"\n\n";
      PointArray b = mwbasis.getbasis();
      if(verbose)
	for (i=0; i<rank; i++)
	  { Point P = b[i];
	    cout << "Generator "<<(i+1)<<" is "<<P<<"; ";
	    cout << "height "<<height(P)<<endl;
	  }
      if(verbose) cout<<"\nregulator is "<<mwbasis.regulator()<<endl;

      if(rank>0) // else don't bother to do extra search
	{
	  if(change) 
	    cout<<"Points on curve " <<(Curve)(CD_orig)<<":";
	  for (i=0; i<rank; i++)
	    { Point P = b[i];
	      if(change) P = shift(P,&CD_orig,u,r,s,t,1);
	      cout << "\nGenerator "<<(i+1)<<" is "<<P<<"; ";
	      cout << "height "<<height(P);
	      if(!P.isvalid()) cout<<" --warning: NOT on curve!";
	      cout << "\n";
	      post_process(P);
	    }
	}
      cout<<endl;
    }
}

int getcurve(Curvedata& CD, Curvedata& CD_orig, 
	     bigint& u, bigint& r, bigint& s, bigint& t, 
	     int& change, int verb)
{
  Curve C0;
  int sing=1; bigint a;
  bigint zero, one; zero=0; one=1;
  while(sing)
    {
      if(verb) cout  << "Enter parameter a (a=3 to exit): ";
      cin >> a;
      if(a==3) return 0;
      cout << "a = " << a << endl;
      C0 = Curve(a,zero,one,zero,zero);
      if(C0.isnull()) return 0;  // quitting condition
      CD_orig = Curvedata(C0,0); // DON'T change coords
      CD = CD_orig.minimalize(u,r,s,t);
      if (verb) cout << endl;
      sing=CD.isnull();
      if(sing) // input curve was singular, non-null
	{
	  cout<<"Curve "<<C0<<" is singular\n";
	}
    }
  if((Curve)CD!=C0)
    {
      if(verb)
	{
	  cout<<"Input curve "<<C0<<"\n";
	  cout<<"Working with minimal curve "<<(Curve)CD<<"\n";
	  cout<<"\t[u,r,s,t] = ["<<u<<","<<r<<","<<s<<","<<t<<"]\n";
	}
      change=1;
    }
  else 
    {
      if(verb) cout << "Curve "<< (Curve)CD<<" :  ";
      change=0;
    }
  if (verb) cout << endl;
  return 1;
}

void post_process(const Point& p) // on [a,0,1,0,0]
{
  bigint x,y,z;
  p.getcoordinates(x,y,z);  // unweighted proj coords
  z=gcd(x,z);
  x/=z;                     // weighted: p=(x/z^2,y/z^3) now
  bigint u = x*x, v=y*z, w=-x*z*z;

  bigint g = gcd(u,gcd(v,w));
  u/=g; v/=g; w/=g;
// Should now have u/v+v/w+w/u = a, i.e. (uuw+vvu+wwv)/uvw = a
  cout << "[u,v,w] = ["<<u<<","<<v<<","<<w<<"]";
  bigint an = u*u*w+v*v*u+w*w*v;
  bigint ad = u*v*w;
  if(div(ad,an))
    cout << " satisfy u/v+v/w+w/u = " << an/ad << endl;
  else
    cout << " --error: u/v+v/w+w/u is not integral!" << endl;
    
}
