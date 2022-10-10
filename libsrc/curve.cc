// curve.cc: implementations of elliptic curve class Curve
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
// 
// This file is part of the eclib package.
// 
// eclib is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// eclib is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with eclib; if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
// 
//////////////////////////////////////////////////////////////////////////
 
// originally adapted from Elliptic.cc by Oisin McGuiness

#include <eclib/curve.h>

//Kraus' conditions:

int valid_invariants(const bigint& c4, const bigint& c6)  
{
  bigint disc = c4*c4*c4; disc -= c6*c6;
  if (sign(disc)==0) return 0;
  if (ndiv(1728,disc)) return 0;  // need c4^3-c6^2=1728D, with D|=0
  long x6= mod(c6,27);
  if((x6==9)||(x6==-9)) return 0; // need c6 != +-9 (mod 27)
  x6 = mod(c6,4);
  if(x6==-1) return 1;            // OK if c6 = -1 (mod 4)
  if(ndiv(16,c4)) return 0;       // else need c4=0 (mod 16)
  x6=mod(c6,32);                  //      and
  return ((x6==0) || (x6==8));    //      c6 = 0,8 (mod 32).
}

void c4c6_to_ai(const bigint& c4, const bigint& c6, 
                bigint& a1, bigint& a2, bigint& a3, bigint& a4, 
                bigint& a6, 
                bigint& b2, bigint& b4, bigint& b6, bigint& b8)
{
//  cout<<"In c4c6_to_ai() with c4="<<c4<<" and c6="<<c6<<endl;
  bigint I12; I12=12;
  b2 = mod(-c6,I12);               //  cout<<"...b2="<<b2<<endl;
  const bigint& b22 = b2*b2;
  b4 = (b22-c4)/24;                //  cout<<"...b4="<<b4<<endl;
  b6 = (-b2*b22+36*b2*b4-c6)/216;  //  cout<<"...b6="<<b6<<endl;
  b8 = (b2*b6 - b4*b4) / 4;        //  cout<<"...b8="<<b8<<endl;

  a1 = (odd(b2) ? 1 : 0);
  a3 = (odd(b6) ? 1 : 0);
  a2 = (b2-a1)/4;         // N.B. a1 == a1*a1 (= 0,1)
  a4 = (b4-a1*a3)/2;
  a6 = (b6-a3)/4;         // N.B. a3 == a3*a3 (= 0,1)
  //  cout<<"...returning a1="<<a1<<", a2="<<a2<<", a3="<<a3
  //      <<", a4="<<a4<<", a6="<<a6<<endl;
}

void c4c6_to_ai(const bigint& c4, const bigint& c6, 
                bigint& a1, bigint& a2, bigint& a3, bigint& a4, 
                bigint& a6)
{
  bigint b2, b4, b6, b8;
  c4c6_to_ai(c4,c6,a1,a2,a3,a4,a6,b2,b4,b6,b8);
}

void minimise_c4c6(const bigint& c4, const bigint& c6, const bigint& discr, 
                   bigint& newc4, bigint& newc6, bigint& newdiscr, bigint& u)
{
  bigint p,g; long a,b,d;
  u = 1; int u_is_1 = 1;
  newc4=c4; newc6=c6;
  const bigint& c62 = sqr(c6);
  newdiscr = (sqr(c4)*c4-c62)/1728; // this must be set before returning
  g=gcd(c4,c6); if(is_one(g)) return;
  g = gcd( c62, newdiscr );  if(is_one(g)) return;
  const vector<bigint>& p_list = pdivs(g);
//  cout<<"g = "<<g<<endl;
  vector<bigint>::const_iterator pr = p_list.begin();
  while ( pr!=p_list.end() )
  {
    p = *pr++;
    d = (long)floor(val(p,g)/12.0);
//    cout<<"With p="<<p<<", initial d="<<d<<endl;
    if (p==2)
      {
	a = mod(c4 >> (4*d) , 16);
	b = mod(c6 >> (6*d) , 32); if(b<0) b+=32;
//	cout<<"a="<<a<<", b="<<b<<endl;
	if (( (b%4)!=3) && !( (a==0) && (( b==0) || (b==8) )))
	  {
	    d--;
	  }
      }
    else if (p==3) if (val(3,c6)==(6*d + 2)) d--;
    if(d>0) {u *= pow(p,d); u_is_1 = 0;}
//    cout<<"With p="<<p<<", final d = "<<d<<", u="<<u<<endl;
  }
  if(u_is_1) return;
  bigint u2, u4, u6, u12;
  mulx(u,u,u2); mulx(u2,u2,u4); mulx(u2,u4,u6); mulx(u6,u6,u12);
  newc4 = c4 / u4;
  newc6 = c6 / u6;
  newdiscr /=  u12;
}

//constructor for curve with invariants as argument
Curve::Curve(const bigint& c4, const bigint& c6)
{
  if (valid_invariants(c4, c6))
    {
      c4c6_to_ai(c4,c6,a1,a2,a3,a4,a6);
    }
  else 
    {
      // cout << " ## attempt to call Curve constructor"
      //      << " with invalid invariants c4 = "<<c4<<", c6 = "<<c6
      //      << ": reading as null curve\n";
        a1=0; a2=0; a3=0; a4=0; a6=0;
    }
}

Curve::Curve(const bigrational& j) // one curve with this j-invariant
{
  if (is_zero(num(j)))
    {
      a1=0; a2=0; a3=1; a4=0; a6=0; // 27a3
    }
  else
    if (num(j)==1728*den(j))
      {
        a1=0; a2=0; a3=0; a4=-1; a6=0; // 32a2
      }
  else
    {
      a1=0; a2=0; a3=0;
      bigint n = num(j);
      bigint m = n-1728*den(j);
      a4 = -3*n*m;
      a6 = -2*n*m*m;
    }
}

void Curve::input(istream& is)
{
  char c;  // to eat commas and detect [ from {;
           // `{' flags curve input by invariants, eg {1,3}
           // `[' by coeffs a1--a6, eg [1,2,3,4,6]
           // seperators and terminators must then be present
           // (any nonnumeric will do after the first { or [ )
           // else assumes a1 a2 a3 a4 a6 separated by whitespace
  is>>skipws;
  is>>c;
  //  cout<<"First char read = "<<c<<"\n";
  switch (c) {
  case '{':
        {
	  //  cout<<"Reading {c4,c6}...\n";
         bigint c4, c6;
         is >> c4 >> c; 
	 if(c!=',')
	   {
	     cout << "syntax error on curve input" << endl;
             return;
	   }
	 is >> c6 >> c;
	 if(c!='}')
	   {
	     cout << "syntax error on curve input" << endl;
             return;
	   }
         if (valid_invariants(c4, c6))
           {
             const bigint& b2 = BIGINT(mod(-c6,12));
             const bigint& b22 = b2*b2;
             const bigint& b4 = (b22-c4)/24;
             const bigint& b6 = (-b2*b22+36*b2*b4-c6)/216;
             a1 = (odd(b2) ? 1 : 0);
             a3 = (odd(b6) ? 1 : 0);
             a2 = (b2-a1*a1)/4;
             a4 = (b4-a1*a3)/2;
             a6 = (b6-a3*a3)/4;
           }
         else 
           {
           cout << " ## invalid invariants, reading as null curve\n";
             a1=0; a2=0; a3=0;a4=0; a6=0;
           }
          }
         break;
       case '[':
	 //	 	 cout<<"Reading [a1,a2,a3,a4,a6]...\n";
         is >> a1 >> c;
	 if(c!=',')
	   {
	     cout << "syntax error on curve input" << endl;
             return;
	   }
	 is >> a2 >> c; 
	 if(c!=',')
	   {
	     cout << "syntax error on curve input" << endl;
             return;
	   }
	 is >> a3 >> c; 
	 if(c!=',')
	   {
	     cout << "syntax error on curve input" << endl;
             return;
	   }
	 is >> a4 >> c; 
	 if(c!=',')
	   {
	     cout << "syntax error on curve input" << endl;
             return;
	   }
	 is >> a6 >> c; 
	 if(c!=']')
	   {
	     cout << "syntax error on curve input" << endl;
             return;
	   }
	 //	 cout<<"["<<a1<<","<<a2<<","<<a3<<","<<a4<<","<<a6<<"]"<<endl;
         break;
       default:
	 //	 	 cout<<"Reading a1 a2 a3 a4 a6 ...\n";
	 is.unget();
         is >> a1 >> a2 >> a3 >> a4 >> a6;
	 //	 cout<<"["<<a1<<","<<a2<<","<<a3<<","<<a4<<","<<a6<<"]"<<endl;
       }
}


// puts out TeX-ed equation of curve
void Curve::tex_print(ostream &os) const
{
        os << "$y^2" ;
        if(a1==0){
                ;
        } else {
                if(a1==1) os << " + xy" ;
                else if(a1==-1) os << " - xy" ;
                else if(a1 > 0) os << " +" << a1 << "xy" ;
                else os << " " << a1 << " xy" ;
        }
        if(a3==0){
                ;
        } else {
                if(a3==1) os << " + y" ;
                else if(a3==-1) os << " - y" ;
                else if(a3 > 0) os << " +" << a3 << "y" ;
                else os << " " << a3 << " y" ;
        }
        os << " = x^3" ;
        if(a2==0){
                ;
        } else {
                if(a2==1) os << " + x^2" ;
                else if(a2==-1) os << " - x^2" ;
                else if(a2 > 0) os << " +" << a2 << "x^2" ;
                else os << " " << a2 << " x^2" ;
        }
        if(a4==0){
                ;
        } else {
                if(a4==1) os << " + x" ;
                else if(a4==-1) os << " - x" ;
                else if(a4 > 0) os << " +" << a4 << "x" ;
                else os << " " << a4 << " x" ;
        }
        if(a6==0){
                ;
        } else {
                if(a6==1) os << " + 1" ;
                else if(a6==-1) os << " - 1" ;
                else if(a6 > 0) os << " +" << a6  ;
                else os << " " << a6 ;
        }
        os << "$" ;
        return ;
}


// end of file: curve.cc

//
// The following functions from OM's code are never used
//
// /* number of solutions to y^2+ay=b mod p */
// long quadroots(const bigint& a, const bigint& b, long p ) 
// {
//         if (p == 2) { 
//            if (!odd(a))  return(1) ; 
//            if (!odd(b))  return(2) ;
//            return(0) ;
//                 } 
//         else {
//          long d = I2long((a*a+4*b)%p);
//          if (d==0) return 1;
//          if (legendre(d,p)==1) return 2; 
//          return 0;
//               }
// }

// // for finding number of points mod 2 and 3
// long pointsmod(long p, const Curve& E)
// {
//   bigint a1,a2,a3,a4,a6;
//   E.getai(a1,a2,a3,a4,a6);
//   if (p == 2)
//     return( 1 + quadroots(a3, a6, 2) +
//            quadroots(a1 + a3, 1 + a2 + a4 + a6, 2)) ;
//   if (p ==3)
//     return( 1 + quadroots(a3, a6, 3) +
//            quadroots(a1 + a3, 1 + a2 + a4 + a6, 3)
//            + quadroots(-a1 + a3, -1 + a2 - a4 + a6, 3)) ;
//  // Now p>3
//   long count=0;
//   for (long x=0;x<p;x++) count+=quadroots(x*a1+a3,x*(x*(x+a2)+a4)+a6, p);
//   return  count;
// }
