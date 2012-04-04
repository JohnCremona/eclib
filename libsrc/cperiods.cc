// cperiods.cc: implementations of class Cperiods and period lattice functions
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
 
#include "cperiods.h"

//#define DEBUG 1

#ifndef MPFP
void swap(bigcomplex& a, bigcomplex& b)
{
  bigcomplex c(a); a=b; b=c;
}
#define SMALL(x) is_zero((x))
#else
#define SMALL(x) is_approx_zero((x))
//#define SMALL(x) (abs(x)<1.0e-14)
//#define SMALL(x) is_zero((x))
#endif

// Reorders 3 complex nos so real parts are decreasing
void reorder1(bigcomplex& a, bigcomplex& b, bigcomplex& c)
{
  if (real(a) < real(c)) swap(a,c);
  if (real(a) < real(b)) swap(a,b);
  else if (real(b) < real(c)) swap(b,c);
}
 
//reorders 3 complex nos so e1 is real if any (
void reorder2(bigcomplex& e1, bigcomplex& e2, bigcomplex& e3)
{
#if(0)
  if (is_real(e1)) return;
  else if (is_real(e2)) {swap(e1,e2); return;}
  else if (is_real(e3)) {swap(e1,e3); return;}
#endif
#if(0)
  cout<<"Entering reorder2() with \n";
  cout<<"e1="<<e1<<"\n";
  cout<<"e2="<<e2<<"\n";
  cout<<"e3="<<e3<<"\n";
#endif
  if(abs(imag(e1))>abs(imag(e3))) {swap(e1,e3);}
  if(abs(imag(e1))>abs(imag(e2))) {swap(e1,e2);}
  else if(abs(imag(e2))>abs(imag(e3))) {swap(e2,e3);}
#if(0)
  cout<<"Leaving reorder2() with \n";
  cout<<"e1="<<e1<<"\n";
  cout<<"e2="<<e2<<"\n";
  cout<<"e3="<<e3<<"\n";
#endif
}
 
bigcomplex cagm1(const bigcomplex& a, const bigcomplex& b);

//Computes periods of a curve given the 3 2-division points (i.e. the three
//roots of the cubic)

// For real curves, here either the ei are real with e1<e2<e3 so a,b,c
// are real, w1 is real and w2 pure imaginary; or e3 is real and
// e1=conj(e2), in which case agm1 is real and w1 is real

void eiperiods(bigcomplex e1, bigcomplex e2, bigcomplex e3, 
               bigcomplex& w1, bigcomplex& w2)
{
  bigcomplex a(sqrt(e3-e1)); 
  bigcomplex b(sqrt(e3-e2)); 
  bigcomplex c(sqrt(e2-e1)); 
#ifdef DEBUG
  cout<<"In eiperiods with a = " << a << ", b = "<<b<< ", c = " << c << endl;
#endif
  bigcomplex agm1 = cagm1(a,b);
  bigcomplex agm2 = cagm1(a,c);
#ifdef DEBUG
  cout<<"agm1=agm(a,b)="<<agm1<<endl;
  cout<<"agm2=agm(a,c)="<<agm2<<endl;
#endif
  bigfloat pi = Pi();
  w1= bigcomplex(pi,to_bigfloat(0))/agm1;
  w2= bigcomplex(to_bigfloat(0),pi)/agm2;
#ifdef DEBUG
  cout<<"Leaving eiperiods with w1 = " << w1 << ", w2 = "<<w2 << endl;    
#endif
}

//#define DEBUG_CUBIC

bigcomplex* solve_nonsingular_cubic(const bigint& c1, const bigint& c2, const bigint& c3) 
//Returns an array of 3 complex roots.
{
#ifdef DEBUG_CUBIC
cout << "In solve_nonsingular_cubic with c1 = "<<c1<<", c2 = "<<c2<<", c3 = "<<c3<<"\n";
#endif
  bigfloat rc1=I2bigfloat(c1);
  bigfloat rc2=I2bigfloat(c2);
  bigfloat rc3=I2bigfloat(c3);
  static bigfloat three(to_bigfloat(3)), 
    two(to_bigfloat(2)), one(to_bigfloat(1));
  bigfloat third = one/three;
  bigcomplex w =  bigcomplex(to_bigfloat(-1), sqrt(three))/two;
  bigint p3=  3*c2 - c1*c1;
  bigint q = c1*(2*sqr(c1)-9*c2)+27*c3;
  bigfloat rq=I2bigfloat(q), rp3=I2bigfloat(p3);
  bigcomplex *roots = new bigcomplex[3];
  long i;
 
#ifdef DEBUG_CUBIC
  cout << "c1 = " << c1 << ", rc1 = " << rc1 << endl;
  cout << "p3 = " << p3 << ", rp3 = " << rp3 << endl;
  cout << "q  = "<<q<<", rq  = "<<rq<<"\n";
#endif

  if (is_zero(p3))       // pure cubic
    { 
#ifdef DEBUG_CUBIC
      //      cout << "In pure cubic case\n";
      //      cout << "About to take cube root of q = " << (q) 
      //      	   << " by pow(-,third) where third = " << third << endl;
#endif
      //      bigcomplex rootq = pow(bigcomplex(rq),bigcomplex(third));
#ifdef DEBUG_CUBIC
      cout << "In pure cubic case\n";
      cout << "About to take cube root of q = " << (q) 
	   << " by exp(log(-)/three) where three = " << three << endl;
#endif
      bigcomplex rootq = exp(log(bigcomplex(rq))/three);
#ifdef DEBUG_CUBIC
      cout << "returns result " << roots[0] << endl;
#endif
      roots[0]=-(rootq+rc1)/three;
      rootq*=w;
      roots[1]=-(rootq+rc1)/three;
      rootq*=w;
      roots[2]=-(rootq+rc1)/three;
    }
  else
    {
      //NB It is important to compute d EXACTLY and then convert to
      //floating point, rather than work with f.p. values for q and
      //p3, since otherwise bad cancellation can occur!
      bigint d =  q*q+ 4*p3*sqr(p3);
      bigfloat rd=I2bigfloat(d);
      bigcomplex t1cubed = (sqrt(bigcomplex(rd)) - rq)/two;
      bigcomplex t2cubed = (sqrt(bigcomplex(rd)) + rq)/two;
#ifdef DEBUG_CUBIC
      //      cout << "d = " << d << "\n";
      //      cout << "About to take cube root of t1cubed = " << t1cubed 
      //      	   << " by power(-,third) where third = " << third << endl;
#endif
      //      bigcomplex t1 = pow(t1cubed,third);
#ifdef DEBUG_CUBIC
      cout << "d = " << d << "\n";
      cout << "sqrt(d) = " << sqrt(bigcomplex(rd)) << endl;
      cout << "sqrt(d)-rq = " << sqrt(bigcomplex(rd))-rq << endl;
      cout << "sqrt(d)+rq = " << sqrt(bigcomplex(rd))+rq << endl;
      cout << "About to take cube root of t1cubed = " << t1cubed 
	   << " by exp(log(-)/three) where three = " << three << endl;
#endif
      bigcomplex t1 = exp(log(t1cubed)/three);
      bigcomplex t2 = exp(log(t2cubed)/three);
#ifdef DEBUG_CUBIC
      cout << "returns result " << t1 << endl;
      cout << "t1^3-t1cubed =  " << t1*t1*t1-t1cubed << endl;
#endif
      if(abs(t1)<abs(t2)) 
	{
	  t1=rp3/t2;
#ifdef DEBUG_CUBIC
	  cout << "resetting t1=p3/t2= " << t2 << endl;
#endif
	}
      roots[0] = (-rc1+t1-rp3/t1)* third;
      t1*=w;
      roots[1] = (-rc1+t1-rp3/t1)* third;
      t1*=w;
      roots[2] = (-rc1+t1-rp3/t1)* third;

      if(d<0) // then all roots should be real so we set this manually
	      // in case rounding disguises this:
	{
	  for(i=0; i<3; i++) roots[i]=real(roots[i]);
	}
    }
  int niter=3;
#ifdef DEBUG_CUBIC
  cout << "refining roots using Newton with " << niter << " iterations\n";
  cout << "unrefined roots: ";
  for(i=0; i<3;i++) cout << roots[i] << "\n";
#endif
  for(i=0; i<3; i++)
    {
      bigcomplex z = roots[i], fz, fdashz;
      for(int iter=0; iter<niter; iter++)
	{
	  fz = ((z+rc1)*z+rc2)*z+rc3;
	  fdashz = (three*z+two*rc1)*z+rc2;
	  if(!is_zero(fdashz)) z -= fz/fdashz;
	}
      roots[i] = z;
    }
#ifdef DEBUG_CUBIC
  cout << "refined roots: ";
  for(i=0; i<3;i++) cout << roots[i] << "\n";
#endif
  return roots;
}
 
// Gets the 3 2-division points given the coefficients 
void getei(const Curvedata& E, bigcomplex& e1, bigcomplex& e2, bigcomplex& e3)
{
  bigint b2,b4,b6,b8;    
  E.getbi(b2,b4,b6,b8);

#ifdef DEBUG
  cout<<"Solving monic cubic with coeffs "<<b2<<","<<(8*b4)<<","<<16*b6<<endl;
#endif
  bigcomplex* ei = solve_nonsingular_cubic(b2,8*b4,16*b6);
#ifdef DEBUG
  cout<<"ei = "<<ei[0]<<","<<ei[1]<<","<<ei[2]<<endl;
#endif
  bigfloat four(to_bigfloat(4));
  e1 = ei[0]/four;  e2 = ei[1]/four;  e3 = ei[2]/four;
#ifdef DEBUG
  cout<<"After rescaling,\n";
  cout<<"ei = "<<e1<<","<<e2<<","<<e3<<endl;
#endif
  delete [] ei;
}

Cperiods::Cperiods(const Curvedata& E)
{
  lattice_type = getconncomp(E);
#ifdef DEBUG
  cout<<"Lattice type = "<<lattice_type<<endl;
#endif
  getei(E,e1,e2,e3);
#ifdef DEBUG
  cout<<"Before reordering,\n";
  cout<<"ei = "<<e1<<","<<e2<<","<<e3<<endl;
#endif
  
  if (lattice_type==2) reorder1(e3,e2,e1); // if all real, make e1<e2<e3
  else reorder2(e3,e2,e1);                 // to make e3 real 
// this ordering ensures that eiperiods will give wR,wRI:
// wR real, and either (type 2) wRI pure imag or (type 1) Re(wRI)=wR/2
#ifdef DEBUG
  if(lattice_type==2) 
    {
      cout << "e1 = " << real(e1) << "\ne2 = " 
	   << real(e2) << "\ne3 = " 
	   << real(e3) << "\n";
      cout << "(all real, e1<e2<e3)"<<endl;
    }
  else
    {
      cout << "e1 = " << (e1) << "\ne2 = " 
	   << (e2) << "\ne3 = " 
	   << (e3) << "\n";
      cout << "(e3 real, e1=conj(e2))"<<endl;
    }
#endif
  
  eiperiods(e1,e2,e3,wR,wRI);
#ifdef DEBUG
  cout << "After eiperiods, \n";
  cout << "wR = " << wR << " (should be real)\n";
  cout << "wRI = " << wRI << " \n";
#endif  
  if(lattice_type==1)
    {
      while(real(wRI)/real(wR)<0) wRI+=wR;   
      while(real(wRI)/real(wR)>1) wRI+=wR;   
      wI = bigcomplex(to_bigfloat(0),2*imag(wRI));
    }
  else
    {
      wI=wRI;
    }
  w1=wR; w2=wRI;
#ifdef DEBUG
  cout << "Before lattice normalization, \n";
  cout << "wR = " << wR << " (should be real)\n";
  cout << "wI = " << wI << " \n";
  cout << "wRI = " << wRI << " \n";
  cout << "real(wRI)/real(wR) = "<<real(wRI)/real(wR)<<endl;
#endif  
  tau = normalize(w2,w1);  // NB reverse params;  from compproc.h
#ifdef DEBUG
  cout << "wR = " << wR << " (should be real)\n";
  cout << "wI = " << wI << " (should be imag)\n";
  if(lattice_type==1)
    cout << "wRI = " << wRI << " (real part should be half wR)\n";
  else
    cout << "wRI = " << wRI << " (real part should be 0)\n";
  cout << "w1 = " << w1 << "\n";
  cout << "w2 = " << w2 << "\n";
  cout << "tau       = "<<tau<<" (abs(tau)="<<abs(tau)<<")\n";
#endif  
  store_sums();
}

void Cperiods::store_sums()
{
  static bigfloat one(to_bigfloat(1));
  qtau = q(tau);
  if(abs(qtau)>0.99)
    {
      cout << "Warning from Cperiods::store_sums: qtau = " 
	   << qtau << " is not small!\n";
    }
  w1squared = w1*w1;
  w1cubed   = w1*w1squared;
  bigcomplex term = one, qtm = qtau;
  sum3=to_bigfloat(0);
  for (bigfloat m=to_bigfloat(1); ! SMALL(term); m+=1)
    { 
      term = qtm*m / (one - qtm);
      qtm *= qtau;
      sum3 += term;
#ifdef DEBUG
      cout<<"term = "<<term<<", sum3 = "<<sum3<<endl;      
#endif
    }
#ifdef DEBUG
  cout<<"final sum3 = "<<sum3<<endl;      
#endif
  sum3 = one/to_bigfloat(12) - to_bigfloat(2)*sum3;
#ifdef DEBUG
  cout<<"stored sum3 = "<<sum3<<endl;      
#endif
 }

//#define DEBUG_XY

bigcomplex Cperiods::X_coord(const bigcomplex& qz)
{
  static bigfloat one(to_bigfloat(1));
  bigcomplex sum(sum3), term(one),  qtm(one), w;
  while ( ! SMALL(term/sum) )
    { w = qtm*qz;
      term = w / pow((one - w), 2);
      qtm *= qtau;
      sum += term;
#ifdef DEBUG_XY
      cout<<"qtm = "<<qtm<<", X-term = "<<term<<", sum = "<<sum<<endl;      
#endif
    }
  term = one; qtm = qtau;
  while ( ! SMALL(term/sum) )
    { w = qtm / qz;
      qtm *= qtau;
      term = w / pow((one - w), 2);
      sum += term;
#ifdef DEBUG_XY
      cout<<"qtm = "<<qtm<<", X-term = "<<term<<", sum = "<<sum<<endl;      
#endif
    }
  bigcomplex ans = sum*TWOPIEYE*TWOPIEYE;
#ifdef DEBUG_XY
  cout<<"X_coord returning ans = "<<ans<<endl;
#endif
  return ans;
}

bigcomplex Cperiods::Y_coord(const bigcomplex& qz)
{
  static bigfloat one(to_bigfloat(1));
  bigcomplex sum(to_bigfloat(0)), term(one), qtm(one), w;
// we are summing w*qt^m(1+w*qt^m)/(1-w*qt^m)^3 for m in Z
// m=0  term
  w=qz;
  sum = w*(one + w) / pow((one - w), 3);
  qtm *= qtau;
#ifdef DEBUG_XY
  cout<<"m=0 gives sum = "<<sum<<endl;      
#endif
// positive m terms;  qtm=qt^m
  while ( ! SMALL(term/sum) )
    { w = qtm*qz;
      term = w*(one + w) / pow((one - w), 3);
      qtm *= qtau;
      sum += term;
#ifdef DEBUG_XY
      cout<<"Y-term = "<<term<<", sum = "<<sum<<endl;      
#endif
    }
// negative m terms;  qtm=qt^n where n=-m
  qtm = qtau; term = one;
  while ( ! SMALL(term/sum) )
    { w = qtm / qz;
      term = w*(one + w) / pow((w - one), 3);
      qtm *= qtau;
      sum += term;
#ifdef DEBUG_XY
      cout<<"Y-term = "<<term<<", sum = "<<sum<<endl;      
#endif
    }
  bigcomplex ans =  sum * TWOPIEYE * TWOPIEYE * TWOPIEYE;
#ifdef DEBUG_XY
  cout<<"Y_coord returning ans = "<<ans<<endl;
#endif
  return ans;
}


void Cperiods::XY_coords(bigcomplex& X, bigcomplex& Y, const bigcomplex& z)
{
#ifdef DEBUG_XY
  cout<<"In XY-coords with z = " << z << endl;
#endif
  // first adjust z w.r.t. lattice [wR,wI]:
  bigcomplex z1 = z;
  z1-=wR*floor(real(z1)/real(wR));
  z1-=wI*floor(imag(z1)/imag(wI));
  z1/=w1;
  bigcomplex qz = q(z1);
  //  while(abs(qz)>0.9) qz*=qtau;
  //  while(abs(qz)>1.1) qz*=qtau;
#ifdef DEBUG_XY
  cout<<"In XY-coords with z = " << z << ", z1 = " << z1 << endl;
  cout<<"qtau = " << qtau << ", qz = " << qz << endl;
  cout<<"abs(qtau) = " << abs(qtau) << ", abs(qz) = " << abs(qz) << endl;
  cout<<"w1 = " << w1 << ", w2 = " << w2 << endl;
#endif
  X = X_coord(qz) / w1squared;
  Y = Y_coord(qz) / w1cubed;
#ifdef DEBUG_XY
  cout<<"XY-coords returns X = " << X << ", Y = " << Y << endl;
#endif
  return;
}

vector<bigcomplex> Cperiods::ellztopoint(const bigcomplex& z, const bigcomplex& a1, const bigcomplex& a2, const bigcomplex& a3)
{
  vector<bigcomplex> xy(2);
  XY_coords(xy[0],xy[1],z);
  xy[0] -= (a1*a1+to_bigfloat(4)*a2)/to_bigfloat(12);
  xy[1] -= (a1*xy[0] + a3); xy[1]/=to_bigfloat(2);
  return xy;
}


//#define DEBUG_CAGM
bigcomplex cagm1(const bigcomplex& a, const bigcomplex& b)
{
  bigcomplex x=a, y=b, oldx;
#ifdef DEBUG_CAGM
  cout<<"cagm("<<x<<","<<y<<"):"<<endl;
#endif
  static bigfloat two=to_bigfloat(2);
  bigfloat theta, piby2=Pi()/two;
  while (1)
    { 
      oldx=x;
      x=(x+y)/two;
      y= sqrt(oldx*y);
      theta = arg(y/x);
      if ((theta>piby2) || (theta<=-piby2)) y=-y;   
#ifdef DEBUG_CAGM
      cout<<"x = "<<x<<"\ty = "<<y<<endl;
      cout<<"Relative error = "<<abs((x-y)/x)<<endl;
#endif
#ifdef MPFP
      if(is_approx_zero(abs((x-y)/x))) return x;
#else
      if(is_zero(abs((x-y)/x))) return x;
#endif
    }
  return x;
}
// end of file cperiods.cc
