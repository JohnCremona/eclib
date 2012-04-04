// periods.cc: test program for computing periods of an elliptic curve
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
 
#include "compproc.h"
 
void reorder1(bigcomplex& e1, bigcomplex& e2, bigcomplex& e3)
{  bigcomplex a=e1, b=e2, c=e3;
   if (real(a) < real(c)) swap(a,c);
   if (real(a) < real(b)) swap(a,b);
     else if (real(b) < real(c)) swap(b,c);
    e1=a; e2=b; e3=c;
}
 
void reorder2(bigcomplex& e1, bigcomplex& e2, bigcomplex& e3)
{ if (is_real(e2)) swap(e1,e2);
  else if (is_real(e3)) swap(e1,e3);
}
 
void eiperiods(bigcomplex e1, bigcomplex e2, bigcomplex e3, 
               bigcomplex& w1, bigcomplex& w2)
{
    bigcomplex a = root(e3-e1,2);
    bigcomplex b = root(e3-e2,2);
    bigcomplex c = root(e2-e1,2);
    w1= bigcomplex(M_PI,0)/cagm(a,b);
    w2= bigcomplex(0,M_PI)/cagm(a,c);
}
 
void getei(bigcomplex a1,bigcomplex a2,bigcomplex a3,bigcomplex a4,bigcomplex a6,
	   bigcomplex& e1, bigcomplex& e2, bigcomplex& e3)
{
  bigcomplex* ei = solvecubic(a2 + a1*a1/4,  a4 + a1*a3/2,  a6 + a3*a3/4);
  e1 = ei[0];
  e2 = ei[1];
  e3 = ei[2];
}
 
void fix(bigcomplex& z, int field)
{ bigfloat x = real(z), y = imag(z);
  if (field > 2) { y/=2; x+=y; }
  if (field > 1) { y*=sqrt(field); }
  z = bigcomplex(x,y);
}
 
void unfix(bigcomplex& z, int field)
{ bigfloat x = real(z), y = imag(z);
  if (field > 1) { y/=sqrt(field); }
  if (field > 2) { x-=y; y*=2; }
  z = bigcomplex(x,y);
}

int field;
bigfloat ra1,ia1,ra2,ia2,ra3,ia3,ra4,ia4,ra6,ia6;

int getfield(void)
{
 cout << "\nField ?  "; cin>>field;
 return   ((field==1) || (field==2) || (field==3) || (field==7) || (field==11)
    || (field==19) || (field==43) || (field==67) || (field==163));
}

int getai(void)
{
 cout << "\na1, a2, a3, a4, a6:";
 cin >> ra1 >> ia1 >> ra2 >> ia2 >> ra3 >> ia3 >> ra4 >> ia4 >> ra6 >> ia6;
 return   !((ra1== 0) && (ia1== 0) && (ra2== 0) && (ia2== 0) &&
	    (ra3== 0) && (ia3== 0) && (ra4== 0) && (ia4== 0) &&
	    (ra6== 0) && (ia6== 0));
}

int main(void)
{
 
 bigcomplex a1,a2,a3,a4,a6,w1,w2;
 bigcomplex tau, c4, c6, b2, b4, b6;
 bigcomplex e1, e2, e3;
 bigfloat area;
 int realcurve, allrealroots;
 
cout << "Computes periods of an elliptic curve given standard\n";
cout << "  coefficients a1,a2,a3,a4,a6 \n\n";
 
while (getfield())
{
while (getai())
  {
    a1=bigcomplex(ra1,ia1); fix(a1,field);
    a2=bigcomplex(ra2,ia2); fix(a2,field);
    a3=bigcomplex(ra3,ia3); fix(a3,field);
    a4=bigcomplex(ra4,ia4); fix(a4,field);
    a6=bigcomplex(ra6,ia6); fix(a6,field);
    
    getei(a1,a2,a3,a4,a6,e1,e2,e3);
    realcurve = is_real(a1) && is_real(a2) && is_real(a3) && is_real(a4) && is_real(a6);
    allrealroots = is_real(e1) && is_real(e2) && is_real(e3);
    
    if (realcurve) { if (allrealroots) reorder1(e3,e2,e1);
    else reorder2(e3,e2,e1);
		   }
    cout << "\ne1:" << e1 << "\ne2:" << e2 << "\ne3:" << e3 << "\n";
    
    eiperiods(e1,e2,e3,w1,w2);
    cout << "w1:" << w1 << "\nw2:" << w2 << "\n"; 
    tau = normalize(w1,w2);
    cout << "w1:" << w1 << "\nw2:" << w2 << "\n";
    cout << "tau:" << tau << "\n";
    area = imag((w1*conj(w2)));
    cout << "Area = " << area << "\n";
    
    getc4c6(w1,w2,c4,c6);  unfix(c4,field); unfix(c6,field);
    cout << "\nc4 and c6 from periods:\n";
    cout << "c4: " << c4 << "\n";
    cout << "c6: " << c6 << "\n";
    
    b2 = a1*a1 +  4*a2;
    b4 = a1*a3 +  2*a4;
    b6 = a3*a3 +  4*a6;
    c4=b2*b2 -  24*b4;  unfix(c4,field);
    c6=-b2*b2*b2 +  36*b2*b4 -  216*b6;  unfix(c6,field);
    
    cout << "\nc4 and c6 from coefficients:\n";
 cout << "c4: " << c4 << "\n";
 cout << "c6: " << c6 << "\n";
}
}
 cout << "\n\n\n";
}
