// holzer.cc:   solution of legendre equations with Holzer reduction


//              Solves aX^2+bY^2+cZ^2=0 (if possible), 
//              then uses Mordell's reduction trick to ensure that 
//              the inequalities Z^2 < ab etc. all hold.

#include "marith.h"
#include "conic.h"

#ifndef VERBOSITY
#define VERBOSITY 1
#endif

int main()
{
  cout<<"Solving ax^2 + by^2 + cz^2 = 0\n\n";

  bigint a_in, b_in, c_in;
  bigint a,b,c,d,g,x0,y0,z0,x1,y1,z1;
  bigint xfac,yfac,zfac;

  while(1) {
  cout << "Enter coefficients a b c: ";
  cin >> a >> b >> c;
  cout<<a<<" "<<b<<" "<<c<<endl; 
  if(a*b*c==0) {break;}
  a_in=a; b_in=b; c_in=c;

  //Step 1: check signs and permute if necessary

  int sa=sign(a), sb=sign(b), sc=sign(c);
  int sx=1, sy=1, sz=1;
  int perm=0;  // =0 if none, 1 if (a,c), 2 if (b,c);
  if((sa==sb)&&(sb==sc))
    {
      cout<<"No solution: signs equal!\n";
      break;
    }
  if(sa==sb)
    {
      if(sc==1) // swap all signs, no permutation needed
	{
	  a=-a; b=-b; c=-c;
	}
    }
  else if(sa==sc)
    {
      if(sb==1) // swap signs
	{
	  a=-a; b=-b; c=-c;
	}
               // swap b<->c
	{
	  d=b; b=c; c=d; 
	  perm=2;
	}
    }
  else if(sb==sc)
    {
      if(sa==1) // swap signs
	{
	  a=-a; b=-b; c=-c;
	}
               // swap a<->c
	{
	  d=a; a=c; c=d; 
	  perm=1;
	}
    }

  //  Now a>0, b>0, c<0
  if(VERBOSITY)
    cout<<"After sorting signs, (a,b,c)=("<<a<<","<<b<<","<<c<<")\n";

  // Step 2: make coprime:
  g=gcd(a,b); if(g>1) g=gcd(g,c);
  if(g>1)
    {
      a/=g; b/=g; c/=g;
    }

  //  Now gcd(a,b,c)=1
  if(VERBOSITY)
    cout<<"After making coprime, (a,b,c)=("<<a<<","<<b<<","<<c<<")\n";

  // Step 3: remove square factors:
  bigint a1,a2,b1,b2,c1,c2;
  bigintArray aplist, bplist, cplist;
  sqfdecomp(a,a1,a2,aplist); a=a1;
  sqfdecomp(b,b1,b2,bplist); b=b1;
  sqfdecomp(c,c1,c2,cplist); c=c1;
  xfac=b2*c2;
  yfac=a2*c2;
  zfac=a2*b2;

  // Now a, b, c are square-free
  if(VERBOSITY)
    cout<<"After removing squares, (a,b,c)=("<<a<<","<<b<<","<<c<<")\n";

  //Step 4: Make pairwise coprime

  g=gcd(a,b);
  if(g>1) {a/=g; b/=g; c*=g; zfac*=g;}
  g=gcd(b,c);
  if(g>1) {a*=g; b/=g; c/=g; xfac*=g;}
  g=gcd(c,a);
  if(g>1) {a/=g; b*=g; c/=g; yfac*=g;}

  // Now abc is square-free
  if(VERBOSITY)
    {
      cout<<"After making pairwise coprime, (a,b,c)=("<<a<<","<<b<<","<<c<<")\n";
      cout<<"gcd(a,b)="<<gcd(a,b)<<endl;
      cout<<"gcd(b,c)="<<gcd(b,c)<<endl;
      cout<<"gcd(c,a)="<<gcd(c,a)<<endl;
    }

  //Step 5: Solve

  int res = solve_conic(a,0,b,-c,x0,z0,y0); // in that order!
  if(VERBOSITY)
    {
      cout<<"Solution 0 to reduced equation: = ("<<x0<<":"<<y0<<":"<<z0<<")\n";
    }
  testsol(a,0,b,-c,x0,z0,y0,VERBOSITY);

  // Check back in original equation:
  x1=x0*xfac;
  y1=y0*yfac;
  z1=z0*zfac;
  if(perm==1) swap(x1,z1); else if(perm==2) swap(y1,z1);
  if(VERBOSITY)
    {
      cout<<"Solution 0 to original equation = ("<<x1<<":"<<y1<<":"<<z1<<")\n";
    }
  testsol(a_in,0,b_in,-c_in,x1,z1,y1,VERBOSITY);
  
  //Step 6: Check if Holzer's conditions are satisfied 
  //        and if not use Mordell reduction

  conic_mordell_reduce(a, b, c, x0, y0, z0,  VERBOSITY);

  
  // Check back in original equation:
  x1=x0*xfac;
  y1=y0*yfac;
  z1=z0*zfac;
  if(perm==1) swap(x1,z1); else if(perm==2) swap(y1,z1);
  if(VERBOSITY)
    {
      cout<<"Reduced solution to original equation = ("
	  <<x1<<":"<<y1<<":"<<z1<<")\n";
    }
  testsol(a_in,0,b_in,-c_in,x1,z1,y1,VERBOSITY);
  }
}

