// FILE HEEGNER.CC: (Needs upgrading) partial Heegner point computation
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2007 John Cremona
// 
// This file is part of the mwrank/g0n package.
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
//
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "oldforms.h"
#include "newforms.h"
#include "cperiods.h"
#include "h1newforms.h"
#include "periods.h"
#include "points.h"
#include "mwprocs.h"
#include "lambda.h"


#define RAT_APPROX_METHOD
#define LAMBDA_METHOD
#define SINGLE
#define DEBUG

long dsearch(const level& N);       //  returns 0 if fails
bigcomplex getz(long n, long d);

void ratapprox2(const bigfloat& x, long max_denom_digits, bigint& a, bigint& b);

// new version doesn't work,theory wrong!
bigcomplex newgetz(long n, long d);  

long geta(long d, long n);

int process(Curvedata* CD, vector<bigfloat> lambdas, long nl, 
             const bigfloat& ra1, const bigfloat& ra2, 
             const bigfloat& ra3, Cperiods& cp, const bigcomplex& z,
             long maxdd, const bigfloat& h);
int newprocess(mw* mwb, const bigfloat& ra1, const bigfloat& ra2, 
             const bigfloat& ra3, Cperiods& cp, const bigcomplex& z,
             const bigfloat& hlim, const bigfloat& h);

int main(void)
{
  set_precision("Enter number of decimal places");
 long maxdiv=1;  bigfloat hlim;  long maxdd=20;
#ifdef RAT_APPROX_METHOD
#ifndef LAMBDA_METHOD  
 cout << "Enter max digits of denominator: "; cin>>maxdd;
#endif
#else
 cout << "Enter max height of point to look for: "; cin>>hlim;
#endif

 int n=1; 
  while (n>0) { cout<<"\n\nEnter level: "; cin>>n;
 if (n>0)
{
 cout << ">>>\nLevel " << n << "\n";
 int usedata=1, verb=1;
 long d = dsearch(level(n));
 if(d==0)
   {
     cout<<"No suitable discriminant!"<<endl; continue;
   }
 else
   cout << "Discriminant = " << d << endl;

 bigcomplex z0 = getz(n,d);
#ifdef DEBUG
 cout << "z0 = " << z0 << endl;
#endif
 bigfloat x0=TWOPI*real(z0);
 bigfloat y0=TWOPI*imag(z0);
 cout << "x0 = " << x0 << "\ny0 = " << y0 << endl;
#ifdef DEBUG
 bigcomplex q0(exp(-y0)*cos(x0),exp(-y0)*sin(x0));
 cout << "exp(2*pi*i*z0) = " << q0 << endl;
#endif
 h1newforms nf(n,1,usedata);

 int i;
#ifdef SINGLE
 cout<<"Which form number? "; cin>>i; i--;
#else
 for(i=0; i<nf.n1ds; i++)
#endif
   {
     cout << "\nForm number " << i+1 << ": " << endl;
     ldash1 ld1(&nf,i);
     ld1.compute();
     int r = ld1.rank();  // The anlytic rank
     if(r!=1) {cout << "Rank = " << r << ", not 1, skipping.\n"; continue;}
     cout << "Rank = 1\n";

     part_period pp(&nf,&(nf.nflist[i]));
     cout << "Computing Heegner point ... " << endl;
     pp.compute(z0);
     bigcomplex z1 = pp.getperiod();
     cout << "z mod L =  " << z1 << endl;

     Cperiods cp = nf.getperiods(i, -1, 0); 
     Curve C = trans_to_curve(cp);
     cp = Cperiods(C);   // recompute to greater precision!
     cout<<"periods are:\n" << cp<<endl;
     bigcomplex w1,w2;
     cp.getwi(w1,w2);
     int ncc = getnorm_code(cp);  // = Number of Connected Components

     Curvedata CD(C,1);  // The 1 causes minimalization
     if(getdiscr(C)!=getdiscr(CD))
       {
	 cout << "Non-minimal curve = \t" << C << ", minimal curve = \t";
       }
     else if(verb) cout << "Curve = \t";
     cout << (Curve)CD << "\t";
     CurveRed CR(CD);
     cout << "N = " << getconductor(CR) << endl;

     long nl;
     vector<bigfloat> Lambda_bad = lambda_bad(CR,nl,0);

     bigint a1,a2,a3,a4,a6;
     CD.getai(a1,a2,a3,a4,a6);
     bigfloat ra1=I2bigfloat(a1), ra2=I2bigfloat(a2), ra3=I2bigfloat(a3);

     bigfloat x,y,rp;
     rp=x=real(w1);y=imag(w2);
     if(ncc==1)	 x/=2;   // L=[2x,x+yi], else L=[x,yi]
     long k = longify(round(ncc*imag(z1)/y));
     bigfloat t = real(z1-(double)k*(w2/(double)ncc));
     t -= rp*floor(t/rp);
     int z_on_other_comp = (ncc==2)&&(odd(k));
     cout << "Real period = " << rp << "\n" << "t = " << t;
     if(z_on_other_comp) cout<<"; on other component";
     cout<<endl;

     bigfloat lf1 = ld1.value();
     long nt = ntorsion(CD);
     long pcp = prodcp(CR);
     bigfloat RS = (lf1*nt*nt) / (2*x*pcp);
     cout << "Expected height of generator = h = " << RS << endl;

     mw mwbasis(&CD,1);  // to process any points found

     cp.norm_region();
     //cout<<"After norm_region(), periods are:\n" << cp<<endl;
     cout << "Enter max value to divide Heegner point by: "; cin>>maxdiv;
     cout << "Now dividing z mod L by 1, ..., "<<maxdiv<<endl;

     int found=0;
     for(k=1; (k<=maxdiv) /* &&(!found) */ ; k++)
       {
	 if(z_on_other_comp&&even(k)) continue;
	 for(int m=0; (m<k) /* &&(!found) */ ; m++)
	   {
	     bigcomplex z2 = (t+(double)m*w1)/(double)k;
	     if((ncc==1)||!z_on_other_comp)
	       {
		 cout<<"(k,m)=("<<k<<","<<m<<"); z2 = "<<z2<<endl;
#ifdef RAT_APPROX_METHOD
                 found=process(&CD,Lambda_bad, nl, ra1,ra2,ra3,cp,z2,maxdd,RS);
#else
                 found=newprocess(&mwbasis,ra1,ra2,ra3,cp,z2,hlim,RS);
#endif
	       }
	     z2+=w2/(double)2;
	     if( /* (!found)&& */
		(ncc==2)&&
		((z_on_other_comp&&odd(k))||(!z_on_other_comp&&even(k))))
	       {
		 cout<<"(k,m)=("<<k<<","<<m<<")b; z2 = "<<z2<<endl;
#ifdef RAT_APPROX_METHOD
                 found=process(&CD,Lambda_bad, nl, ra1,ra2,ra3,cp,z2,maxdd,RS);
#else
                 found=newprocess(&mwbasis,ra1,ra2,ra3,cp,z2,hlim,RS);
#endif
	       }
	   }
       }
   }    // end of newform loop
}       // end of if(n>0)
}       // end of while()
}       // end of main()

long dsearch(const level& N)  // uses global level::modulus; returns 0 if fails
{
  const vector<long>& plist = N.plist;
  int found=0;
  long i,j,d = 0;

// The following should be in REVERSE order of size
  long dlist[] = {-163,-67,-43,-19,-11,-8,-7,-4,-3};
   for(i=0; (i<9)&&!found; i++)
   { 
     d = dlist[i];
     int okay = 1;
     for(j=0; (j<plist.size())&&okay; j++)
       okay = (kronecker(d,plist[j]) == +1);
     found = okay;
   }
  if(!found)d=0;
  return d;
}
 
// to solve x^2-dy^2=m with y maximal, d<0:

int nsolve(long d, long m, long& x, long& y)
{
  int ok=0;
  y=longify(floor(sqrt(-bigfloat(m)/d))); x=-1;
  while((!ok)&&(y>=0))
    {
      long x2=m+d*y*y;
      x=(long)(sqrt((double)x2)+0.1);
      if(x*x==x2) ok=1;
      else y--;
    }
  return ok;
}

bigcomplex getz(long n, long d)
{
  long a = geta(d,n);
  cout << "a = " << a << endl;
  if(d%4) return bigcomplex(1-2*a,sqrt(-bigfloat(d)))/(double)(2*n);
  else    return bigcomplex( -a, sqrt(-bigfloat(d)/4))/(double)(n);
}

bigcomplex newgetz(long n, long d)
{
  if(d%4)
    {
      long x,y;
      nsolve(d,4*n,x,y);
      return bigcomplex(x,y*sqrt(bigfloat(-d)))/(double)(2*n);
    }
  else
    {
      long x,y,d0 = d/4;
      nsolve(d0,n,x,y);
      return bigcomplex(x,y*sqrt(bigfloat(-d0)))/(double)n;
    }
}

long geta(long d, long n) // Solve Norm(a-w) = k*n 
{ 
  long eps = posmod(d,4);
  long d0 = (d-eps)/4;
  long ans = 0;
  for(long a=1; (a<=n)&&(ans==0); a++)
    if ( (((a*(a-eps)-d0) % n) == 0) ) ans = a;
  return ans;
}
 
int newprocess(mw* mwb, const bigfloat& ra1, const bigfloat& ra2, 
             const bigfloat& ra3, Cperiods& cp, const bigcomplex& z,
             const bigfloat& hlim, const bigfloat& h)
{
 bigcomplex X,Y;
 bigfloat RX,RY;

 cp.XY_coords(X, Y, z);
 RX = real(X) - (ra1*ra1 + 4*ra2)/12;
 RY = (real(Y) - ra1*RX - ra3)/2;
 cout<<"Real point [x,y] = ["<<RX<<","<<RY<<"]"<<endl;
 
 bigfloat x_eps = 1.0e-10;

 mwb->search_range(RX-x_eps,RX+x_eps,hlim);
 return (mwb->getrank())>0;
}

int process(Curvedata* CD, vector<bigfloat> lambdas, long nl, 
             const bigfloat& ra1, const bigfloat& ra2, 
             const bigfloat& ra3, Cperiods& cp, const bigcomplex& z,
             long maxdd, const bigfloat& h)
{
  int res=0;
 bigcomplex X,Y;
 bigfloat RX,RY;

 cp.XY_coords(X, Y, z);
 RX = real(X) - (ra1*ra1 + 4*ra2)/12;
 RY = (real(Y) - ra1*RX - ra3)/2;

 bigint nx,dx,ny,dy;
 Point P(CD);
#ifdef LAMBDA_METHOD
 res=make_point_from_x_and_ht(CD,lambdas,RX,h,&P);
#else
 res=make_point_from_x(CD,RX,maxdd,&P);
#endif
 if(res)
   {
    bigfloat hp=height(P);
    if(is_approx_zero(hp))
      {
	return 0;
      }
    bigfloat r = hp/h;
    bigint k = Iround(sqrt(r));
    cout << "\n***!!!*** Rational P = " << P 
         << " with height " << hp
         << "\n = h * "<< r
         << ", P = "<< k <<" * generator\n\n";
    return 1;
   }
 cout<<"Real point [x,y] = ["<<RX<<","<<RY<<"]"<<endl;
 return 0;
}

int is_small_enough(const bigfloat& x)
{
  return abs(x)<1.0e-30;
}

void ratapprox2(const bigfloat& x, long max_denom_digits, bigint& a, bigint& b)
{
cout<<"In ratapprox2 with x = " << x << endl;
  bigint c, x0, x1, x2, y0, y1, y2;
  bigfloat xx, diff, one; one=1;
  xx = x; x0 = 0; x1 = 1; y0 = 1; y1 = 0;
  diff = 1;
  bigint maxdenom = pow(BIGINT(10),max_denom_digits);
  while ( !is_small_enough(diff) && (y2<maxdenom))
    { c = Iround( xx );
      x2 = x0 + c*x1; x0 = x1; x1 = x2;
      y2 = y0 + c*y1; y0 = y1; y1 = y2;
      diff = abs( x - (I2bigfloat(x2)/I2bigfloat(y2)) );
cout<<"x2,y2,diff = " << x2 << ", " << y2 << ", " << diff << endl;
      if ( is_approx_zero(abs(xx - I2bigfloat(c))) ) diff = 0;
      else xx = one/(xx - I2bigfloat(c));
    }
  a = x2; b = y2;
  if ( b < 0 ) {a=-a; b=-b; }
}


