// kbessel.cc: implementation of K-Bessel function for arbitrary real
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2005 John Cremona
// 
// This file is part of the mwrank package.
// 
// mwrank is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// mwrank is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with mwrank; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// 
//////////////////////////////////////////////////////////////////////////
 //             parameter nu. Adapted from the PARI function kbessel.

#include <iostream>  // for debugging only
#include <cmath>
using namespace std;
#include "kbessel.h"
#ifndef M_LN2
#define M_LN2   0.69314718055994530942
#endif
#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif
#define LOG2 M_LN2
#define PI M_PI
#define expo(x) ((x)==0?-1000:(long)floor(log(fabs((x)))/LOG2))
#define lbin -64 // = bit-precision of doubles (or should be)

double kbessel(double nu, double gx, int debug)
{
  double  x,y,p1,p2,zf,zz,s,t,q,r,u,v,e,f,c,d,ak;
  double  nu2,w;
  long    k,k2,n2,n,ex;

  if(debug) {
    cout << "\nCalled kbessel("<<nu<<","<<gx<<").\n";
  }

  k=0;x=gx; 
  nu2=-4*nu*nu;
  n=(long)(32*LOG2+PI*nu/2);
  n2=(n<<1);
  if (x<n)
    {
      if(debug) {
        cout << "In the x<n case.\n";
      }
      zf=sqrt(PI/n2);
      zz=1.0/(n2<<2);
      s=1.0;t=0.0;k2=2*n2+1;
      for (k=n2;k>0;--k)
        {
          k2-=2;
          p1=k2*k2+nu2;
          ak=-p1*zz/k;  
          s=1+ak*s;
          t=k2+ak*t;
        }
      u=s*zf;  t/=2.0;
      v=-(t*zf+u*nu)/n2;
      q=n2;     r=x+x;

      if(debug) {
        cout << "Finished k loop.  lbin = "<<lbin<<endl;
        cout << "t,u,v,q,r = "<<t<<", "<<u<<", "<<v<<", "<<q<<", "<<r<<endl;
      }

      do
        {
          p1=5.0/q;
          if (expo(p1)>= -1) p1=0.5;
          p2=1.0-r/q;   if (p1>p2) p1=p2;
          c=-p1; k=1; d=1; e=u; f=v; 
          if(debug) {
            cout << "...outer loop: p1 = "<<p1<<", expo(p1) = "<<expo(p1)<<endl;
          }
          do
            {
              w=(((k-0.5)*u)+((q-k)*v));
              w+=nu*(u-(v+v));
              u=q*v/k;
              v=w/k;
              d*=c;
              e+=d*u; 
              f+=(p1=(d*v));
              k++;
              ex=expo(p1)-expo(f);
              if(debug) {
                cout << "......inner loop: f = "<<f<<", expo(f) = "<<expo(f)<<endl;
                cout << "                  p1= "<<p1<<", expo(p1) = "<<expo(p1)<<", ex = "<<ex<<endl;
              }
            }
          while(ex>lbin);
          if(debug) {
            cout << "...finished inner loop -- ex = "<<ex<<endl;
          }
          p1=u;u=e;e=p1;p1=v;v=f;f=p1;
          q*=(1+c);
          p1=q-r; ex=expo(p1);
          if(debug) {
            cout << "...outer loop -- ex = "<<ex<<endl;
            cout << "u,q,r,p1 = " << u<<", "<<q<<", "<<r<<", "<<r<< endl;
          }
        }
      while(ex>lbin);
      y=u*pow((x/n),nu);
    }
  else
    {
      if(debug) {
        cout << "In the x>=n case.\n";
      }
      p2=2*x;
      zf=sqrt(PI/p2);
      zz=1.0/(4*p2);
      s=1.0; k2=2*n2+1;
      for (k=n2;k>0;--k)
        {
          k2-=2;
          p1=k2*k2+nu2;
          ak=(p1*zz)/k;
          s=1.0-(ak*s);
        }
      y=s*zf;
    }
  y/=exp(x);
  if(debug) {
    cout << "kbessel returns " << y << endl;
  }
  return y;
}

