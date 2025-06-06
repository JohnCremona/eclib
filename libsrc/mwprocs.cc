// mwprocs.cc: implementation of class mw for Mordell-Weil basis
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2023 John Cremona
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


//#define DEBUG_QSIEVE

#include <eclib/interface.h>
#include <eclib/compproc.h>

#include <eclib/method.h>

#include <eclib/points.h>
#include <eclib/polys.h>
#include <eclib/curvemod.h>
#include <eclib/pointsmod.h>
#include <eclib/ffmod.h>
#include <eclib/divpol.h>
#include <eclib/tlss.h>
#include <eclib/elog.h>
#include <eclib/saturate.h>

#include <eclib/sieve_search.h>

#include <eclib/mwprocs.h>

//
// some locally called general functions, belong in library maybe:
//

// unlikely to be called by anything but find_inf:
vector<bigcomplex> roots_of_cubic(const Curve& E)
{
  bigint a1,a2,a3,a4,a6;
  E.getai(a1,a2,a3,a4,a6);

  bigfloat ra1=I2bigfloat(a1),
           ra2=I2bigfloat(a2),
           ra3=I2bigfloat(a3),
           ra4=I2bigfloat(a4),
           ra6=I2bigfloat(a6);

  bigfloat c1 = ra2 + ra1*(ra1/4);
  bigfloat c2 = ra4 + ra1*(ra3/2);
  bigfloat c3 = ra6 + ra3*(ra3/4);
  return solvecubic(c1,c2,c3);
}

bigfloat min_real(vector<bigcomplex> array)
{
//cout<<"In min_real() with array:\t"<<array<<endl;
  bigfloat minr, r; int first=1; minr=0;
  for (unsigned int i=0; i<array.size(); i++)
    { if(abs(imag(array[i]))<0.001)  // then the root is regarded as real
	{
	  r = real(array[i]);
	  if (first||(minr > r)) {minr = r; first=0;}
	}
    }
//cout<<"minr finally " << minr << "\n";
  return minr;
}

int order_real_roots(vector<double>& bnd, vector<bigcomplex> roots);
//checks (and returns) how many roots are actually real, and puts those in
//bnd, in increasing order, by calling set_the_bound
int set_the_bounds(vector<double>& bnd, bigfloat x0, bigfloat x1, bigfloat x2);
//This transforms (if possible) x0, x1 and x1 into double;  the search
//should be made on [x0,x1]U[x2,infty] so if x1 or x2 overflows, the search
//is on [x0,infty].  The function returns 3 in the first case, 1 in the second.
//If x0 overflows, it returns 0.  A warning is printed out.

bigfloat det(const vector<vector<bigfloat>>& m, long m_size);
   // fwd declaration: det and detminor jointly recursive

vector<vector<bigfloat>> get_minor(const vector<vector<bigfloat>>& m, long m_size, long i0, long j0)
{
  vector<vector<bigfloat>> the_minor(MAXRANK, vector<bigfloat>(MAXRANK));
  for (long i=0; i<m_size-1; i++)
    {
      long ii = (i>=i0? i+1 : i);
      for (long j=0; j<m_size-1; j++)
      {
	long jj = (j>=j0? j+1 : j);
	the_minor[i][j] = m[ii][jj];
      }
    }
  return the_minor;
}

bigfloat det_minor(const vector<vector<bigfloat>>& m, long m_size, long i0, long j0)
{
  return det(get_minor(m,m_size,i0,j0), m_size-1);
}

bigfloat det(const vector<vector<bigfloat>>& m, long m_size)
{
  bigfloat one = to_bigfloat(1);
  switch (m_size) {
  case 0:
    return one; break;
  case 1:
    return m[0][0]; break;
  case 2:
    return m[0][0]*m[1][1] - m[0][1]*m[1][0];
    break;
  default:
    // Make a copy since we will change m before recursing
    vector<vector<bigfloat>> m1 = m;
    int neg = 0;
    long i,j,i0;
    bigfloat pivot=m1[0][0], piv, temp, eps=to_bigfloat(1.0e-6);
    for(i0=0; i0<m_size && abs(pivot)<eps; i0++)
      pivot=m1[i0][0];
    if(i0==m_size)
      return to_bigfloat(0); // first column all 0
    if(i0>0)  // swap rows 0, i0:
      {
	neg = 1;
	for(j=0; j<m_size; j++)
	  {
	    temp=m1[i0][j]; m1[i0][j]=m1[0][j]; m1[0][j]=temp;
	  }
      }
    // eliminate first column
    pivot=m1[0][0];
    for(i=1; i<m_size; i++)
      {
	piv=m1[i][0]/pivot;
	for(j=0; j<m_size; j++)
	  m1[i][j] -= m1[0][j]*piv;
      }
    if (neg) pivot = -pivot;
    return pivot*det_minor(m1,m_size,0,0);
  }
  return one;  // shouldn't get here
}

//#define DEBUG 1

vector<long> cleardenoms(vector<bigfloat> alpha)
{
  long len = alpha.size();
  vector<long> nlist(len);  // returned
  vector<long> dlist(len);
  long i, lcmd = 1;
  bigfloat x, last=alpha[len-1];
  for (i=0; i < len-1; i++)  // i doesn't include rank (new value)
    { x = alpha[i] / last;
      ratapprox(x, nlist[i], dlist[i]);
      lcmd = (lcmd*dlist[i]) / ::gcd(lcmd, dlist[i]);    // ie lcm(d, dlist[i])
                                        // ie we find the lcm of whole of dlist
#ifdef DEBUG
      cout<<"ratapprox: of "<< x <<" is "<<nlist[i]<<" / "<<dlist[i]<<endl;
      cout<<"  lcm of denoms so far: "<<lcmd<<endl;
#endif
    }
  for (i=0; i < len-1; i++)
    nlist[i] *= (lcmd / dlist[i]);  // clear the denominators
  nlist[len-1] = lcmd;
  return nlist;
}

mw::mw(Curvedata *EE, int verb, int pp, int maxr)
  :E(EE), rank(0), maxrank(maxr), reg(to_bigfloat(1)), verbose(verb), process_points(pp), satsieve(EE,1,verb)
{
#ifdef DEBUG
  verbose=1;
#endif
  height_pairs.resize(MAXRANK, vector<bigfloat>(MAXRANK));
}

// NB We cannot use the default parameter mechanism as this must fit
// the template for the virtual function declared in class
// point_processor!
int mw::process(const bigint& x, const bigint& y, const bigint& z)
{
  return process(x,y,z,MAXSATPRIME);
}

int mw::process(const bigint& x, const bigint& y, const bigint& z, int sat)
{
#ifdef DEBUG
  cout<<"mw::process with x = "<< x <<", y = "<<y<<", z = "<<z<<endl;
#endif
  bigint rz; isqrt(z,rz);
  bigint x1=x*rz, y1=y, z1=z*rz;
  if(iso)
    {
      y1 -= (a1*x1+4*a3*z1);
      x1 *= 2;
      z1 *= 8;
    }
  Point P(E, x1,y1,z1);
  if(P.isvalid()) return process(P,sat);

  // error:
  cerr<<"Raw point       x,y,z = "<<x<<", "<<y<<", "<<z<<endl;
  cerr<<"converted point x,y,z = "<<x1<<", "<<y1<<", "<<z1<<"\t";
  cerr<<"--not on curve!"<<endl;
 return 0;
}

int mw::process(const vector<Point>& Plist, int sat)
{
  // process the points without saturation, do that at the end

  if(verbose)
    cout<<"Processing "<<Plist.size()<<" points ..."<<endl;
  // cout<<"mw::E is "<<E<<endl;
  int flag=0; // flags whether rank has increased
  for( const auto& P : Plist)
    {
      // ensure all P are processed whatever the flag values
      int flag1 = process(P,0);
      flag = flag || flag1;
    }

  if(verbose)
    cout<<"Finished processing the points (which had rank "<<rank<<")"<<endl;

  if((sat>0)&&(rank>0))
    {
      if (verbose) cout<<"saturating up to "<<sat<<"..."<<flush;
      satsieve.set_points(basis);
      long index;
      vector<long> unsat;
      int sat_ok = satsieve.saturate(unsat, index, sat);
      if(verbose)
        {
          cout<<"done";
          if (!sat_ok)
            {
              cout<<" (saturation failed for "<<unsat<<")";
            }
          cout<<endl;
        }
      if(index>1)
	{
	  basis = satsieve.getgens();
	  if(verbose)
	    cout<<"Gained index "<<index<<", new generators = "<<basis<<endl;
	}
// compute the height pairing matrix and regulator
      int i, j;
      for (i=0; i < rank; i++)
	{
	  mat_entry(i,i) = height(basis[i]);
	  for (j=0; j < i; j++)
	    {
	      mat_entry(i,j)
		= mat_entry(j,i)
		= height_pairing(basis[i], basis[j]);
	    }
	}
      reg = det(height_pairs,rank);
      if(verbose)
	cout<<"Regulator =  "<<reg<<endl;
    }
  return flag;
}

int mw::process(const Point& PP, int sat)
{
#ifdef DEBUG
  cout<<"mw::process with P = "<< PP <<endl;
#endif
  Point P(E, PP);  // (1) copy so we can process const points
                   // (2) make sure P's curve is on this mw's curve
  if(!P.isvalid()) cerr << "###Runtime error### P="<<P<<" not on curve "<<(Curve)(*E)<< " in mw::process()\n";
  long ord = order(P);
  long i, j, rank1=rank+1;
#ifdef DEBUG
  cout<<"P = "<< P <<" has order "<<ord<<endl;
  bigfloat hP=height(P);
  cout<<"P = "<< P <<" has height "<<hP<<endl;
#endif

  if (verbose)
    {cout<<"P"<<rank1<<" = "<<P;
#ifdef DEBUG
    cout<<" (height "<<height(P)<<")";
#endif
    cout << "\t" << flush;
#ifdef DEBUG
    cout << "\n";
#endif
    }

  if (ord > 0)
    {
      if (verbose) cout<<" is torsion point, order "<<ord<<endl;
      return 0;
    } // we're not interested in torsion points

  if(!process_points)
    {
      basis.push_back(P); rank++;
      if(verbose) cout<<"P = "<<P<<", ht(P) = "<<height(P)<<endl;
      return 0;
    }

  if (rank==0)  // first non-torsion point
    {
      reg = height(P);
      height_pairs[0][0] = reg;
      basis.push_back(P); rank=1;
      if (verbose) cout<<"  is generator number 1\n";
#ifdef DEBUG
      cout << "first non-torsion point P = "<<P<<", reg so far = h(P) = "<<reg<<"\n";
      //      cout << "returning "<<(maxrank<2)<<endl;
#endif
      if(sat>0)
	{
          satsieve.set_points(basis);
          if (verbose) cout<<"saturating up to "<<sat<<"..."<<flush;
          long index;
          vector<long> unsat;
          int sat_ok = satsieve.saturate(unsat, index, sat);
          if(verbose)
            {
              cout<<"done";
              if (!sat_ok)
                {
                  cout<<" (saturation failed for "<<unsat<<")";
                }
              cout<<endl;
            }
          if(index>1)
            {
              basis = satsieve.getgens();
              if(verbose) cout<<"Gained index "<<index<<", new generator = "<<basis[0]<<endl;
              reg = height(basis[0]);
              mat_entry(0,0) = reg;
            }
	}
      return (maxrank<2); // 1 if max reached
    }

  // otherwise general procedure:

#ifdef DEBUG
  cout<<"additional non-torsion point..."<<endl;
#endif
  //  update the height pairing matrix (at least for now)
  // but don't add point yet
  mat_entry(rank,rank) = height(P); // also sets height in P
  for (i=0; i < rank; i++)
    {
      Point Q = basis[i];
      bigfloat hp = height_pairing(Q, P);
      mat_entry(i,rank) = hp;
      mat_entry(rank,i) = hp;
    }

  // compute cofactors of new last column
  vector<bigfloat> alpha(rank1);  // to store cofactors
  long detsign = ( odd(rank) ? +1 : -1 ); //set for flip before first use
  for (i=0; i < rank; i++)
    {
      detsign = -detsign;
      alpha[i] = det_minor(height_pairs, rank1, i, rank) * detsign;
#ifdef DEBUG
      cout<<"alpha["<<i<<"] = "<<alpha[i]<<"\n";
#endif
    }
  alpha[rank] = reg;  // ie the previous value, before this point
#ifdef DEBUG
  cout<<"alpha["<<rank<<"] = "<<alpha[rank]<<"\n";
#endif

  // find the new determinant
  bigfloat newreg = to_bigfloat(0);
  for (i=0; i <= rank; i++) newreg += mat_entry(i,rank) * alpha[i];
#ifdef DEBUG
  cout<<"After adding P, new height pairing matrix:\n";
  for(i=0; i<=rank; i++)
    {
      for(j=0; j<=rank; j++) cout << mat_entry(i,j) << "\t";
      cout << "\n";
    }
  cout<<"\nreg is now " << newreg << "\t";
#endif

  // test for simple case, new point is indep previous
  if ( abs(newreg/reg) > 1.0e-4 )
    {
      reg = newreg;
#ifdef DEBUG
      cout << "treating as NON-zero" << "\n";
#endif
      basis.push_back(P); rank=rank1;
      if (verbose) cout<<"  is generator number "<<rank<<endl;
      if(sat>0)
	{
          satsieve.reset_points(basis);
          if (verbose) cout<<"saturating up to "<<sat<<"..."<<flush;
          long index;
          vector<long> unsat;
          int sat_ok = satsieve.saturate(unsat, index, sat);
          if(verbose)
            {
              cout<<"done, index = "<<index<<".";
              if (!sat_ok)
                {
                  cout<<" (saturation failed for "<<unsat<<")";
                }
              cout<<endl;
            }
          if(index>1)
            {
              basis = satsieve.getgens();
              if(verbose) cout<<"Gained index "<<index<<", new generators = "<<basis<<endl;
              // completely recompute the height pairing matrix
              for (i=0; i < rank; i++)
                {
                  mat_entry(i,i) = height(basis[i]);
                  for (j=0; j < i; j++)
                    {
                      mat_entry(i,j)
                        = mat_entry(j,i)
                        = height_pairing(basis[i], basis[j]);
                    }
                }
              reg /= (index*index);
            }
	}
#ifdef DEBUG
      cout << "about to return, rank = "<<rank<<", maxrank = "<<maxrank<<": returning "<<(maxrank==rank)<<endl;
#endif
      return (maxrank==rank); // 1 if max reached
    }

  // otherwise, express new point as lin-comb of the previous Now that
  // we are saturating as we go, this should not happen (unless the
  // index is divisible by a prime > sat)
#ifdef DEBUG
  cout << "treating as ZERO" << "\n";
  cout << "Finding a linear relation between P and current basis\n";
#endif

  vector<long> nlist = cleardenoms(alpha);
  long index = nlist[rank];

#ifdef DEBUG
  cout<<"index = "<<index<<endl;
#endif

  // test simple case when new is just Z-linear comb of old
  if ( index == 1 )
    { if (verbose)
	{
	  cout<<" = "<<-nlist[0]<<"*P1";
	  for(i=1; i<rank; i++)
	     cout<<" + "<<-nlist[i]<<"*P"<<(i+1);
	  cout<<" (mod torsion)\n";
	}
// Check:
    Point Q = P;
    for(i=0; i<rank; i++) Q += nlist[i]*basis[i];
    int oQ = order(Q);
    if(oQ==-1) // infinite order, problem!
      {
	cout<<"Problem in mw::process(), bad linear combination!\n";
	cout<<"Difference = "<<Q<<" with height "<<height(Q)<<endl;
      }
    else if(verbose>1)
      {
	cout<<"Difference = "<<Q<<" with height "<<height(Q)<<endl;
      }

    return 0;  // with regulator and basis (and h_p matrix) unchanged
    } // end of if(index==1)

  // otherwise add P to the list now, compute to gain index
  basis.push_back(P);

  if (verbose)
    {
      for (i=0; i < rank; i++)
	cout<<nlist[i]<<"*P"<<(i+1)<<" + ";
      cout<<index<<"*"<<"P"<<(rank1)<<" = 0 (mod torsion)\n";
    }
// Check:
  Point Q = index*P;
  for(i=0; i<rank; i++) Q += nlist[i]*basis[i];
  int oQ = order(Q);
  if(oQ==-1) // infinite order, problem!
    {
      cout<<"Problem in mw::process(), bad linear combination!\n";
      cout<<"Difference = "<<Q<<" with height "<<height(Q)<<endl;
    }
  else if(verbose>1)
    {
      cout<<"Difference = "<<Q<<" with height "<<height(Q)<<endl;
    }


  // find minimum coeff. |ai|
  long min = 0, ni;
  long imin = -1;
  for (i=0; i <= rank; i++)
    { ni = abs(nlist[i]);
      if (  (ni>0) && ( (imin==-1) || (ni<min) )  )
	{ min = ni; imin = i; }
    }

  // find aj with ai ndiv aj, write aj = ai*q + r with 0<r<|ai|
  // since then ai*Pi + aj*Pj = ai(Pi + q*Pj) + r*Pj,
  // replace generator Pi with Pi + q*Pj
  // and replace aj by r, and i by j
  // after finite no. steps obtain some minimal ai = 1
  // then we can just discard Pi.

  long r, q;
  while ( min > 1 )
    {
      for (i=0; i <= rank; i++)
      {
	r = mod(nlist[i], min);
	q = (nlist[i] - r) / nlist[imin];
	if ( r!=0 )
	{
	  basis[imin] += q*basis[i];
	  imin=i; min=abs(r); nlist[imin]=r;
#ifdef DEBUG
	  if (verbose)
	    { for (j=0; j < rank; j++)
		cout<<nlist[j]<<"*"<<basis[j]<<" + ";
	      cout<<nlist[rank]<<"*"<<basis[rank]<<" = 0 (mod torsion)\n";
	    }
#endif
  	  break;  // out of the for loop, back into the while
	}
      } // ends for
    } // ends while

  if (verbose)
    cout<<"Gaining index "<<index<<"; ";

  basis.erase(basis.begin()+imin);

  // completely recompute the height pairing matrix
  for (i=0; i < rank; i++)  // 19/8/02: this was <=rank
    { mat_entry(i,i) = height(basis[i]);
      for (j=0; j < i; j++)
	{
	  mat_entry(i,j)
	    = mat_entry(j,i)
	    = height_pairing(basis[i], basis[j]);
	}
    }

  reg /= (index*index);
  if (verbose)
    {
      cout<<"\nNew set of generators: \n";
      for(i=0; i<rank; i++)
	{
	  if(i)cout<<", ";
	  cout<<"P"<<(i+1)<<" = "<< basis[i];
	}
      cout<<endl;
    }
  return 0; // rank did not increase
} // end of function mw::process(Point)

int mw::saturate(long& index, vector<long>& unsat, long sat_bd, long sat_low_bd)
{
  if (verbose) cout<<"saturating basis (reg="<<reg<<")..."<<flush;

  // This code does a dummy call to index_bound() in order to get
  // points from the search done there into the mwbasis.  But it was
  // decided to let the user do some searching if relevant instead
#if(0)
  vector<Point> pts;
  bigint ib = index_bound(E,basis,pts,0,(verbose>1));
  // Must make sure that the new points have the correct Curvedata pointer!
  for(unsigned int i=0; i<pts.size(); i++)
    pts[i].init(E,getX(pts[i]),getY(pts[i]),getZ(pts[i]));
  bigfloat oldreg = reg;
  int oldrank=rank;
  process(pts,0); //no saturation here!  This may update basis
  bigint ind = Iround(sqrt(oldreg/reg));
  if(verbose&&(ind>1))
    {
      cout<<"after search, gained index "<<ind
	  <<", regulator = "<<reg<<endl;
    }
  if((rank>oldrank))
    {
      cout<<"after search, rank increases to "<<rank
	  <<", regulator = "<<reg<<endl;
    }
#endif
  satsieve.set_points(basis);
  int ok = 1;
  if(rank>0) ok=satsieve.saturate(unsat,index,
                                  sat_bd, sat_low_bd,
                                  10);  // maxntries before dividing
  if(verbose) cout<<"done"<<endl;
  if(!ok)
    cout<<"Failed to saturate MW basis at primes "<<unsat<<endl;
  if(index>1)
    {
      basis = satsieve.getgens();
  // completely recompute the height pairing matrix
      for (int i=0; i < rank; i++)
	{
	  mat_entry(i,i) = height(basis[i]);
	  for (int j=0; j < i; j++)
	    {
	      mat_entry(i,j)
		= mat_entry(j,i)
		= height_pairing(basis[i], basis[j]);
	    }
	}
      reg /= (index*index);
      if(verbose)
	{
	  cout<<"Gained index "<<index<<endl;
	  cout<<"New regulator =  "<<reg<<endl;
	}
    }
  return ok;
}

void mw::search(bigfloat h_lim, int moduli_option, int verb)
{
#ifdef DEBUG
  cout<<"In mw::search;  maxrank = "<<maxrank<<endl;
#endif
  if(moduli_option)
    {
      sieve s(E, this, moduli_option, verb);
      s.search(h_lim);
    }
  else // use Stoll's sieve, Sophie Labour's conversion:
    {
      vector<bigint> c(4);
      E -> getai(a1,a2,a3,a4,a6);
      iso = !((a1==0)&&(a3==0));
      if(iso)
	{
	  c[0]=16*getb6(*E);
	  c[1]= 8*getb4(*E);
	  c[2]=   getb2(*E);
	  c[3]=1;
	}
      else
	{
	  c[0]=a6;
	  c[1]=a4;
	  c[2]=a2;
	  c[3]=1;
	}
      if(iso) h_lim+=2.08;
//    if(iso) cout<<"Adding log(8) to h_lim, increasing it to "<<h_lim<<endl;
      qsieve s(this, 3, c, h_lim, verb);
      bigcomplex c1(I2bigfloat(c[2])),
                 c2(I2bigfloat(c[1])),
                 c3(I2bigfloat(c[0]));
      vector<bigcomplex> rts=solvecubic(c1,c2,c3);
      vector<double> bnd(3);
      int nrr=order_real_roots(bnd,rts);
#ifdef DEBUG_QSIEVE
      cout<<endl;
      cout<<"cubic "<<c[0]<<" "<<c[1]<<" "<<c[2]<<" "<<c[3]<<endl;
      cout<<"coeff "<<c1<<" "<<c2<<" "<<c3<<endl;
      cout<<"roots "<<rts<<endl;
      cout<<"bnd "<<bnd<<endl;
      cout<<"smallest "<<bnd[0]<<endl;
#endif
      s.set_intervals(bnd,nrr,1);
      s.search(); //searches and processes
    }
}

void mw::search_range(bigfloat xmin, bigfloat xmax, bigfloat h_lim,
		    int moduli_option, int verb)
{
  sieve s(E, this, moduli_option, verb);
  s.search_range(xmin,xmax,h_lim);
}

//#define DEBUG_SIEVE

sieve::sieve(Curvedata * EE, mw* mwb, int moduli_option, int verb)
: E(EE), mwbasis(mwb), verbose(verb)
{
  E->getai(a1,a2,a3,a4,a6);
  int ncomp = getconncomp(*E);
  posdisc = ncomp==2;
  long i, j;

// find pt of order two in E(R) with minimal x-coord

  vector<bigcomplex> rts = roots_of_cubic(*E);
  if(posdisc)
    {
      x1=real(rts[0]);
      x2=real(rts[1]);
      x3=real(rts[2]);
      orderreal(x3,x2,x1);  // so x1<x2<x3
      xmin=x1;
    }
  else
    x3=xmin = min_real(rts);

  if (verbose)
    {
      cout << "sieve: real points have ";
      if(posdisc) cout<<x1<<" <= x <= " << x2 << " or ";
      cout << x3 << " <= x; xmin =  " << xmin << endl;
    }
// set up list of auxiliary moduli
// and a list of which residues modulo each of these are squares

  switch(moduli_option) {
  case 1:
    num_aux = 10;
    auxs = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31};
    break;

  case 2:// the following taken from Gebel's scheme
    num_aux = 3;
    auxs = {5184, 5929, 4225}; // 5184 = (2^6)*(3^4), 5929 = (7^2)*(11^2), 4225 = (5^2)*(13^2)
    break;

  case 3:
  default:
    num_aux = 9;
    auxs = {32, 9, 25, 49, 11, 13, 17, 19, 23};
    break;
  }

#ifdef DEBUG_SIEVE
  if(verbose)
    cout<<"Using "<<num_aux<<" sieving moduli:\n" << auxs <<endl;
#endif

  xgood_mod_aux.resize(num_aux);
  squares.resize(num_aux);
  amod.resize(num_aux);

  for (i = 0; i < num_aux; i++)
    {
      long aux = auxs[i];
      squares[i].resize(aux, 0);
      for (j = 0; 2*j < aux; j++) squares[i][(j*j)%aux]=1;
      xgood_mod_aux[i].resize(aux);
    }  // end of aux loop

#ifdef DEBUG_SIEVE
  if(verbose)
    {
      cout<<"squares lists:\n";
      for(i=0; i<num_aux; i++)
	{
	  cout << auxs[i]<<":\t";
	  for(j=0; j<auxs[i]; j++) if(squares[i][j]) cout<<j<<"\t";
	  cout<<endl;
	}
    }
#endif

// variables for collecting efficiency data:

  modhits.resize(num_aux, 0);
  ascore=0; npoints=0;
}

void sieve::search(bigfloat h_lim)
{
// N.B. On 32-bit machines, h_lim MUST be < 21.48 else exp(h_lim)>2^31
//      and overflows
//      On 64-bit machines, h_lim must be < 43.668.

  long i,j;

  // set initial bounds for point coefficients
  alim = I2long(Ifloor(exp(h_lim)));
  clim = clim1 = clim2 = clim0 = I2long(Ifloor(exp(h_lim / 2)));

  if(posdisc)
    {
      if(x2<-1)
	{
	  long temp = I2long(Ifloor(sqrt(alim/(-x2))));
	  if(clim1>temp) clim1=temp;
	}
      if(x1>1)
	{
	  long temp = I2long(Ifloor(sqrt(alim/x1)));
	  if(clim1>temp) clim1=temp;
	}
    }
  if (x3>1)
    {
      long temp = I2long(Ifloor(sqrt(alim/x3)));
      if(clim2>temp) clim2=temp;
    }
  clim=clim2;
  if(posdisc) if(clim1>clim2) clim=clim1;

  if (verbose)
    cout<< "sieve::search: trying a up to "<<alim<<" and c up to "<<clim<<endl;

// declare and initialize other loop variables
  long pd1,pd2,pd3,pd4,pd6, csq, aux;
  cflag.resize(10000);  // max needed; only used up to c each time,
                        // and only when c<=10000 (use_cflag==1)

//
// MAIN LOOP
//

#define FIRSTC 1   // for debugging etc

  for (c = FIRSTC; c <= clim; c++)
    {
      // some preliminary calculations of multiples of c etc.
      csq = c*c /* long */;
      c2 = csq  /* bigint */;
      c3 = c*c2; c4 = c2*c2; c6 = c2*c4;
      d1 = a1*c; d2 = a2*c2; d3 = a3*c3; d4 = a4*c4; d6 = a6*c6;

#ifdef DEBUG_SIEVE
  if(verbose)
    {
      cout<<"c = "<<c<<"\n";
      cout<<"d1,...,d6 = "<<d1<<", "<<d2<<", "<<d3<<", "<<d4<<", "<<d6<<"\n";
    }
#endif
      use_cflag = (c<=10000);
      if(use_cflag)
	{
// set up flag array of residues coprime to c
	  cflag[0]=(c==1);
	  for(i=1; 2*i<=c; i++)
            cflag[i] = cflag[c-i] = (::gcd(i,c)==1);
	}

      // set the main flag matrix
      for (long index = 0; index < num_aux; index++)
	{
	  aux = auxs[index];
          {
            pd1 = posmod(d1 , aux);
            pd2 = posmod(d2 , aux);
            pd3 = posmod(d3 , aux);
            pd4 = posmod(d4 , aux);
            pd6 = posmod(d6 , aux);

            long dddf= 24%aux;
            long ddf = posmod(2*(pd1*pd1)%aux + 8*pd2+24 , aux);
            long df  = posmod(pd1*(pd1+2*pd3)%aux + 4*(pd4+pd2+1) , aux);
            long f   = posmod(((pd3*pd3)%aux+4*pd6) , aux);

            auto flag = xgood_mod_aux[index].begin();
            auto sqs = squares[index].begin();
            long x=aux;
#ifdef DEBUG_SIEVE
            if(verbose)
              {
                cout<<"aux = "<< aux <<"\n ";
                cout<<"pd1,...,pd6 = "<<pd1<<", "<<pd2<<", "<<pd3<<", "<<pd4<<", "<<pd6<<"\n";
              }
#endif
            while(x--)
              {
                *flag++ = sqs[f];
#ifdef DEBUG_SIEVE
                if(verbose)
                  {
                    cout<<"x = "<< aux-x-1 <<", f(x) = "<<f<<": flag = "<<sqs[f]<<"\n";
                  }
#endif
                f   +=  df; if(f  >=aux) f  -=aux;
                df  += ddf; if(df >=aux) df -=aux;
                ddf +=dddf; if(ddf>=aux) ddf-=aux;
              }
          }  // end of default case
 	}  // end of aux loop
#ifdef DEBUG_SIEVE
  if(verbose)
    {
      for(i=0; i<num_aux; i++)
	{
	  cout<<"possible x mod "<<auxs[i]<<": ";
	  for(j=0; j<auxs[i]; j++) if(xgood_mod_aux[i][j]) cout<<j<<"\t";
	  cout<<endl;
	}
    }
#endif
  if(verbose>1)
    {
      for(i=0; i<num_aux; i++)
	{
	  int n=0;
	  cout<<"Number of possible x mod "<<auxs[i]<<": ";
	  for(j=0; j<auxs[i]; j++) if(xgood_mod_aux[i][j]) n++;
	  cout<<n<<" ("<<100.0*n/auxs[i]<<" percent)";
	  cout<<endl;
	}
    }

      // set up for a-loop(s)

      if(posdisc&&(c<=clim1))
	{
	  long amin = -alim, amax = alim;
	  long temp = I2long(Ifloor(csq*x1));
	  if(temp>amin) amin=temp;
	  temp = I2long(Ifloor(csq*x2));
	  if(temp<amax) amax=temp;
	  a_search(amin,amax);
	}

      if(c<=clim2)
	{
	  long temp = I2long(Ifloor(csq*x3));
	  long amin = -alim, amax = alim;
	  if(temp>amin) amin=temp;
	  a_search(amin,amax);
	}
    } // ends c- loop
} // end of sieve::search()

void sieve::search_range(bigfloat xmin, bigfloat xmax, bigfloat h_lim)
{
// N.B. h_lim MUST be < 21.48 else exp(h_lim)>2^31 and overflows
  long i;

  // set initial bounds for point coefficients
  alim = I2long(Ifloor(exp(h_lim)));
  clim = clim1 = clim2 = clim0 = I2long(Ifloor(exp(h_lim / 2)));
  long temp;

  if(xmax<-1)
    {
      temp = I2long(Ifloor(sqrt(alim/(-xmax))));
      if(clim1>temp) clim1=temp;
    }
  if(xmin>1)
    {
      temp = I2long(Ifloor(sqrt(alim/xmin)));
      if(clim1>temp) clim1=temp;
    }
  clim=clim2;
  if(clim1>clim2) clim=clim1;

  if (verbose)
    cout<< "sieve::search: trying a up to "<<alim<<" and c up to "<<clim<<endl;

// declare and initialize other loop variables
  long pd1,pd2,pd3,pd4,pd6, csq, aux;
  cflag.resize(10000);  // max needed; only used up to c each time,
                        // and only when c<=10000 (use_cflag==1)

//
// MAIN LOOP
//

  for (c = 1; c <= clim; c++)
    {
      if(c>clim1) continue;

      // some preliminary calculations of multiples of c etc.
      csq = c*c /* long */;
      c2 = csq  /* bigint */;
      c3 = c*c2; c4 = c2*c2; c6 = c2*c4;
      d1 = a1*c; d2 = a2*c2; d3 = a3*c3; d4 = a4*c4; d6 = a6*c6;

      long amin = -alim, amax = alim;
      temp = I2long(Iceil(csq*xmin));
      if(temp>amin) amin=temp;
      temp = I2long(Ifloor(csq*xmax));
      if(temp<amax) amax=temp;
cout<<"amin = " << amin << ", amax = " << amax <<  endl;

      if(amin>amax) continue;   // skip this c

      if((amax-amin)<10)        // don't both sieving for a
	{
	  a_simple_search(amin,amax);
	}
      else
	{
// set up flag array of residues coprime to c
	  use_cflag = (c<=10000);
	  if(use_cflag)
	    {
	      cflag[0]=(c==1);
	      for(i=1; 2*i<c; i++)
                cflag[i] = cflag[c-i] = (::gcd(i,c)==1);
	    }

// set the main flag matrix
	  for (long index = 0; index < num_aux; index++)
	    {
	      aux = auxs[index];
	      pd1 = posmod(d1 , aux);
	      pd2 = posmod(d2 , aux);
	      pd3 = posmod(d3 , aux);
	      pd4 = posmod(d4 , aux);
	      pd6 = posmod(d6 , aux);

	      long dddf= 24%aux;
	      long ddf = posmod(2*(pd1*pd1)%aux + 8*pd2+24 , aux);
	      long df  = posmod(pd1*(pd1+2*pd3)%aux + 4*(pd4+pd2+1) , aux);
	      long f   = posmod(((pd3*pd3)%aux+4*pd6) , aux);

              auto flag = xgood_mod_aux[index].begin();
	      auto sqs = squares[index].begin();
	      long x=aux;
	      while(x--)
		{
		  *flag++ = sqs[f];
		  f   +=  df; if(f  >=aux) f  -=aux;
		  df  += ddf; if(df >=aux) df -=aux;
		  ddf +=dddf; if(ddf>=aux) ddf-=aux;
		}
	    }  // end of aux loop
	  a_search(amin,amax);
	}
    } // ends c- loop
} // end of sieve::search() version 2 (explicit xmin, xmax)

void sieve::a_search(const long& amin, const long& amax)
{
  bigint pb,qb,db,rdb,rdb2,b,ac;
  long i, a=amin;
  a--;
  if (verbose) cout<<"sieve::search: trying c = "<<c<<"\t"
                   <<"("<<amin<<" <= a <= "<<amax<<")"<<endl;

  for (i=0; i < num_aux; i++)  amod[i] = posmod(a, auxs[i]);
  amodc = posmod(a,c);
#ifdef DEBUG_SIEVE
  if(verbose)
    {
      cout<<"Initial a =  "<<a<<" modulo moduli = \t";
      for(i=0; i<num_aux; i++) cout << amod[i]<<"\t";
      cout<<endl;
    }
#endif

  while (a < amax)
    {
      a++;
      // check that a is good for all the auxiliaries
      amodc++; if(amodc==c) amodc=0;
      int try_x;
      if(use_cflag) try_x = cflag[amodc];
      else          try_x = (::gcd(a,c)==1);

      if(try_x)
	{
	  ascore++;
	}
      // DON'T add "else continue; (with next a)
      // as the amod[i] are not yet updated!
      for ( i=num_aux-1; (i>=0); i--)
	{ long& amodi = amod[i];
	  amodi++;
	  if (amodi == auxs[i]) amodi = 0;
	  if(try_x)
	    {
	      try_x = xgood_mod_aux[i][amodi];
	      if(!try_x) modhits[i]++;
	    }
	}
      if (!try_x) continue;

      pb=a; pb*=d1; pb+=d3;              //    pb = a*d1 + d3;
                                         //    qb = d6 + a*(d4 + a*(d2 + a));
      qb=a; qb+=d2; qb*=a; qb+=d4; qb*=a; qb+=d6;
      db = sqr(pb); db += (4*qb);
      if(isqrt(db,rdb))
	{
	  b = rdb-pb; b/=2; ac = a*c;
	  Point P(*E, ac, b, c3);
	  mwbasis->process(P);
	  npoints++;
	}
    } // ends a-loop
}

void sieve::a_simple_search(const long& amin, const long& amax)
{
  bigint pb,qb,db,rdb,rdb2,b,ac;
  long a;
  if (verbose) cout<<"sieve::search: trying c = "<<c<<"\t"
                   <<"("<<amin<<" <= a <= "<<amax<<")\n";

  for (a=amin; a<amax; a++)
    {
//    pb = a*d1 + d3;
      pb=a; pb*=d1; pb+=d3;
//    qb = d6 + a*(d4 + a*(d2 + a));
      qb=a; qb+=d2; qb*=a; qb+=d4; qb*=a; qb+=d6;
      db = sqr(pb); db += (4*qb);
      if(isqrt(db,rdb))
	{
	  b = rdb-pb; b/=2; ac = a*c;
	  Point P(*E, ac, b, c3);
	  mwbasis->process(P);
	  npoints++;
	}
    } // ends a-loop
}


void sieve::stats(void)
{
  cout << "\nNumber of points found: "<<npoints<<"\n";
  cout << "\nNumber of a tested: "<<ascore<<"\n";
  cout<<"Numbers eliminated by each modulus:\n";
  long nmodhits=0;
  for(long i=0; i<num_aux; i++)
    {
      cout<<auxs[i]<<": "<<modhits[i]<<"\n";
      nmodhits+=modhits[i];
    }
  cout<<"Number eliminated by all moduli: "<<nmodhits<<"\n";
  bigfloat eff = to_bigfloat(nmodhits*100.0)/(ascore-npoints);
  cout<<"Sieve efficiency: "<<eff<<"\n";
}

int order_real_roots(vector<double>& bnd, vector<bigcomplex> roots)
{//checks (and returns) how many roots are actually real, and puts those in
 //bnd, in increasing order, by calling set_the_bound
  long i,nrr=0;
  vector<bigfloat> real_roots;

  for (i=0;i<3;i++)
    {
      if (is_approx_zero(roots[i].imag()))
	{
	  real_roots.push_back(roots[i].real());
	  if (is_approx_zero(real_roots[nrr]))  real_roots[nrr]=0;
	  nrr++;
	}
    }
  //  cout<<"nrr = "<<nrr<<endl;
  //  cout<<"real_roots = "<<real_roots<<endl;
  //  cout<<"Now ordering them..."<<endl;
  switch (nrr)
    {
    case 1: // possible overflow in assignment from bigfloat to double
	return !doublify(real_roots[0],bnd[0]);
	break;
    case 3:
      orderreal(real_roots[2],real_roots[1],real_roots[0]);
      return set_the_bounds(bnd,real_roots[0],real_roots[1],real_roots[2]);
    default:
      cout<<"mw_info::set_the_bounds: two real roots for the elliptic curve...\n";
    }
  return 0; //we should not get here...
}

//This transforms (if possible) x0, x1 and x2 into double; the search
//should be made on [x0,x1]U[x2,infty] so if x1 or x2 overflows, the
//search is made on [x0,infty].  The function returns 3 in the first
//case, 1 in the second.  If x0 overflows, it returns 0.  A warning is
//printed out.
int set_the_bounds(vector<double>& bnd, bigfloat x0, bigfloat x1, bigfloat x2)
{
  if (doublify(x0,bnd[0]))
    {
      cout<<"##WARNING##: lowest bound "<<x0<<" is not a double.\n";
      cout<<"Search will be made over [-height,height]."<<endl;
      return 0;
    }
  else
    {
      if (doublify(x1,bnd[1]) || doublify(x2,bnd[2]))
	{
	  cout<<"##WARNING##: second or third root is not a double.\n";
	  cout<<"]x2,x3[ not excluded in search."<<endl;
	  return 1;
	}
      else
	{
	  return 3;
	}
    }
}


//end of file mwprocs.cc
