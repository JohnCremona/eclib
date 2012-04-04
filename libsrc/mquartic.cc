// mquartic.cc:   Implementation of class quartic and related functions
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
 
#include "mquartic.h"

// constructors
quartic::quartic() 
{
  have_zpol=0; equiv_code=0;
  roots=new bigcomplex[4];
  //cout<<"Quartic constructor #1: " << this << endl;
}

quartic::quartic(const bigint& qa, const bigint& qb, const bigint& qc, 
		 const bigint& qd, const bigint& qe, 
		 bigcomplex* qr,	 int qt,
		 const bigint& qi,const bigint& qj,const bigint& qdisc)
:a(qa),b(qb),c(qc),d(qd),e(qe),type(qt),ii(qi),jj(qj),disc(qdisc)
{ 
  have_zpol=0; equiv_code=0;
  roots=new bigcomplex[4]; 
  for(int i=0; i<4; i++) roots[i] = qr[i];
  //cout<<"Quartic constructor #2: " << this << ", roots="<<roots<< endl;
}

quartic::quartic(const bigint& qa, const bigint& qb, const bigint& qc, 
		 const bigint& qd, const bigint& qe)
:a(qa),b(qb),c(qc),d(qd),e(qe)
{
  have_zpol=0; equiv_code=0;
  roots=new bigcomplex[4]; 
  set_roots_and_type();
}

//#define DEBUG_ROOTS

void quartic::set_roots_and_type()
{
  ii = II(a,b,c,d,e);
  jj = JJ(a,b,c,d,e);
  disc = 4*pow(ii,3)-sqr(jj);
  bigint  H = H_invariant(a,b,c), R = R_invariant(a,b,c,d);
  bigint Q = H*H-16*a*a*ii;  // = 3*Q really
  bigfloat xH = I2bigfloat(H);
#ifdef DEBUG_ROOTS
  bigint diff = -H*H*H + 48*a*a*ii*H - 64*a*a*a*jj - 27*R*R;
  cout<<"H = "<<H<<", R = "<<R<<", diff = "<<diff<<endl;
  if(is_zero(diff)) cout<<"Syzygy satisfied.\n";
  else cout<<"Syzygy NOT satisfied.\n";
#endif
  int nrr;
  if(disc<0) 
    {type=3; nrr=2;}       // 2 real roots
  else 
    {
      if((sign(H)<0)&&(sign(Q)>0)) 
	{type=2; nrr=4;}   // 4 real roots
      else 
	{type=1; nrr=0;}   // 0 real roots
    }
#ifdef DEBUG_ROOTS
  cout<<"Type = " << type << " ("<<nrr<<" real roots)\n";
#endif
  bigcomplex c1(to_bigfloat(0)), c2(-3*I2bigfloat(ii)), c3(I2bigfloat(jj));
  vector<bigcomplex> cphi = solvecubic( c1, c2, c3);
#ifdef DEBUG_ROOTS
  cout<<"Roots of cubic are "<<cphi<<endl;
#endif
  bigfloat a4=4*I2bigfloat(a), xb=I2bigfloat(b);
  bigfloat oneover4a = to_bigfloat(1)/a4;
#ifdef DEBUG_ROOTS
  cout<<"a4 = "<<a4<<", xb = "<<xb<<endl;
#endif
  
  if(type<3)
    {
#ifdef DEBUG_ROOTS
      cout<<"Positive discriminant\n";
#endif
      // all the phi are real;  order them so that a*phi[i] decreases
      bigfloat phi1 = real(cphi[0]);
      bigfloat phi2 = real(cphi[1]);
      bigfloat phi3 = real(cphi[2]);
      if(a>0)      orderreal(phi1,phi2,phi3); 
      else         orderreal(phi3,phi2,phi1); 

#ifdef DEBUG_ROOTS
      cout<<"phi = "<<phi1<<", "<<phi2<<", "<<phi3<<"\n";
      cout<<"xH = " << xH << ", a4*phi3 = " << a4*phi3 <<"\n";
#endif
      if(type==2) // all roots are real
	{
#ifdef DEBUG_ROOTS
	  cout<<"Type 2\n";
#endif
	  bigfloat r1 = safe_sqrt((a4*phi1-xH)/3);
	  bigfloat r2 = safe_sqrt((a4*phi2-xH)/3);
	  bigfloat r3 = safe_sqrt((a4*phi3-xH)/3);
	  if(R<0)  r3 = -r3;
#ifdef DEBUG_ROOTS
	  cout<<"r_i = "<<r1<<", "<<r2<<", "<<r3<<"\n";
	  cout<<"product = "<<r1*r2*r3<<", R = "<<R<<endl;
#endif
	  roots[0] = bigcomplex(( r1 + r2 - r3 - xb) * oneover4a);
	  roots[1] = bigcomplex(( r1 - r2 + r3 - xb) * oneover4a);
	  roots[2] = bigcomplex((-r1 + r2 + r3 - xb) * oneover4a);
	  roots[3] = bigcomplex((-r1 - r2 - r3 - xb) * oneover4a);
	  // Those are all real and in descending order of size
	}
      else // no roots are real
	{
#ifdef DEBUG_ROOTS
	  cout<<"Type 1\n";
#endif
	  bigfloat r1  = safe_sqrt((a4*phi1-xH)/3);
	  bigfloat ir2 = safe_sqrt((xH-a4*phi2)/3);
	  bigfloat ir3 = safe_sqrt((xH-a4*phi3)/3);
	  if(R>0)  r1 = -r1;
#ifdef DEBUG_ROOTS
	  cout<<"r_i = "<<r1<<", "<<ir2<<"i, "<<ir3<<"i\n";
	  cout<<"product = "<<-r1*ir2*ir3<<", R = "<<R<<endl;
#endif
	  roots[0] = bigcomplex( r1-xb,  ir2 - ir3) * oneover4a;
	  roots[1] = conj(roots[0]); // bigcomplex( r1-xb,  ir2 - ir2) * oneover4a;
	  roots[2] = bigcomplex(-r1-xb,  ir2 + ir3 ) * oneover4a;
	  roots[3] = conj(roots[2]); // bigcomplex(-r1-xb, -ir2 - ir3 ) * oneover4a;
	}
    }
  else // disc < 0
    {
#ifdef DEBUG_ROOTS
      cout<<"Negative discriminant\nType 3\n";
#endif
      bigfloat realphi;       // will hold the real root, which will be cphi[2]
      if (is_real(cphi[1])) 
	{
	  realphi=real(cphi[1]);
	  cphi[1]=cphi[2];
	  cphi[2]=realphi;
	}
      else 
	if (is_real(cphi[2])) 
	  {
	    realphi=real(cphi[2]);
	  }
	else 
	  {
	    realphi=real(cphi[0]);
	    cphi[0]=cphi[2];
	    cphi[2]=realphi;
	  }
#ifdef DEBUG_ROOTS
      cout<<"Sorted roots of cubic (real one last) are \n";
      cout<<cphi[0]<<"\n"<<cphi[1]<<"\n"<<cphi[2]<<endl;
#endif
      bigfloat three(to_bigfloat(3));
      bigcomplex r1 = sqrt((a4*cphi[0]-xH)/three);
      bigfloat r3   = safe_sqrt((a4*realphi-xH)/three);
      if(R<0)  r3 = -r3;

#ifdef DEBUG_ROOTS
      cout<<"r_i = "<<r1<<", "<<conj(r1)<<", "<<r3<<"\n";
      cout<<"product = "<<r1*conj(r1)*r3<<", R = "<<R<<endl;
#endif
      roots[0] = bigcomplex( r3 - xb, 2*imag(r1) ) * oneover4a;
      roots[1] = conj(roots[0]);
      roots[2] = bigcomplex(( 2*real(r1) - r3 - xb)) * oneover4a;
      roots[3] = bigcomplex((-2*real(r1) - r3 - xb)) * oneover4a;
      // roots[2] and roots[3] are real
    }
#ifdef DEBUG_ROOTS
  cout << "finished setting roots of quartic.\n";
  dump(cout);
#endif  
}

quartic::~quartic() 
{ 
  //cout<<"Quartic destructor: " << this << ", roots="<<roots<<endl;
  delete[] roots;
}

quartic::quartic(const quartic& q)
:a(q.a),b(q.b),c(q.c),d(q.d),e(q.e),type(q.type),ii(q.ii),jj(q.jj),disc(q.disc)
{ 
  have_zpol=0; equiv_code=q.equiv_code;
  roots=new bigcomplex[4];
  for(int i=0; i<4; i++) roots[i] = q.roots[i];
  //cout<<"Quartic constructor #3: " << this << endl;
}

// member functions & operators

void quartic::assign(const bigint& qa, const bigint& qb, const bigint& qc, 
		     const bigint& qd, const bigint& qe)
{
  have_zpol=0; equiv_code=0;
  a=qa; b=qb; c=qc; d=qd; e=qe; 
  set_roots_and_type();
}

void quartic::assign(const bigint& qa, const bigint& qb, const bigint& qc, 
		     const bigint& qd, const bigint& qe, 
		     bigcomplex* qr,	     int qt,
		     const bigint& qi,const bigint& qj,const bigint& qdisc)
{ 
  have_zpol=0; equiv_code=0;
  a=qa; b=qb; c=qc; d=qd; e=qe; 
  for(int i=0; i<4; i++) roots[i] = qr[i];
  type=qt; ii=qi; jj=qj; disc=qdisc;
//  cout<<"Quartic assign, now: "; dump(cout);
}

void quartic::operator=(const quartic& q)
{ 
  have_zpol=0; equiv_code = q.equiv_code;
  //cout<<" Quartic op=, LHS was: "; dump(cout);
  //cout<<" RHS = "; q.dump(cout);
  a=q.a; b=q.b; c=q.c; d=q.d; e=q.e; 
  for(int i=0; i<4; i++) roots[i] = q.roots[i];
  type=q.type; ii=q.ii; jj=q.jj; disc=q.disc;
  //cout<<" Quartic op=, LHS now: "; dump(cout);
}

int quartic::trivial() const // Checks for a rational root 
{
//    cout << "Checking triviality...\n";
      int i,found;
      bigint num; bigfloat realroot;
      int start = (type==1)? 5 :
                  (type==2)? 1 : 3;
//    cout << "Start = " << start << "\n"; 
      bigint ac = a*c, a2d = a*a*d, a3e = a*a*a*e;
      bigfloat ra = I2bigfloat(a);
      for (i = start, found=0; i<=4 &&  ! found; i++)
         { //  cout << "i= "<<i<<", root = "<<(roots[i-1])<<endl;
            realroot = real((roots)[i-1]);
            num = Iround(ra*realroot);
            found = (((((num+b)*num+ac)*num+a2d)*num+a3e)== 0);
         }
//    cout << "...returning "<<found<<endl;
      return(found);
}

void quartic::make_zpol()
{
  if(have_zpol) return;
  bigint b2 = sqr(b);
  asq=sqr(a);
  p = -H_invariant(a,b,c);
  psq = sqr(p);
  r = R_invariant(a,b,c,d);
  have_zpol=1;
}

// find the number of roots of aX^4 + bX^3 + cX^2 + dX + e = 0 (mod p)
// except 4 is returned as 3 so result is 0,1,2 or 3.
long quartic::nrootsmod(long p) const
{
#ifdef TEST_EQCODE
  cout << "Counting roots mod " << p << " of " << (*this) << "\n";
#endif
  long ap = mod(a,p);
  long bp = mod(b,p);
  long cp = mod(c,p);
  long dp = mod(d,p);
  long ep = mod(e,p);
#ifdef TEST_EQCODE
  cout << "reduced coefficients: " << ap <<","<< bp <<","<< cp <<","<< dp <<","<< ep << "\n";
#endif

  long nroots = (ap==0);  // must count infinity as a root!
  for (long i = 0; (i < p)&&(nroots<3) ; i++)
    {
      long temp = ((((ap*i+bp)*i + cp)*i + dp)*i + ep);
      if ((temp%p)==0) {nroots++;} 
    }
  if(nroots==4) return 3; 
#ifdef TEST_EQCODE
  cout << "returning code " << nroots << "\n";
#endif
  return nroots;
}

unsigned long quartic::set_equiv_code(const vector<long>& plist)
{
#ifdef TEST_EQCODE
  cout << "Setting equiv_code for " << (*this) << "\n";
#endif
  equiv_code=0;
#ifdef NEW_EQUIV // else leave all codes 0, i.e. disable this test
  for(unsigned long i=0; i<plist.size(); i++)
    {
      int code = nrootsmod(plist[i]);
      equiv_code |= (code<<(2*i));
    }
#endif
#ifdef TEST_EQCODE
  cout << "Final code = " << equiv_code << "\n";
#endif
  return equiv_code;
}
