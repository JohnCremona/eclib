// saturate.cc: implementation of class saturator for sieving E(Q)/pE(Q)
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
 
#include <eclib/matrix.h>
#include <eclib/subspace.h>

#include <eclib/points.h>
#include <eclib/polys.h>
#include <eclib/curvemod.h>
#include <eclib/pointsmod.h>
#include <eclib/ffmod.h>
#include <eclib/divpol.h>
#include <eclib/tlss.h>
#include <eclib/elog.h>
#include <eclib/sieve_search.h>
#include <eclib/mwprocs.h>
#include <eclib/saturate.h>
#include <eclib/egr.h>
#include <eclib/htconst.h>

// If point search bound is greater than this, output a warning
// message and reduce to this value:
const int max_search_bound = 18;

// How many auxiliary primes to use without the image increasing in
// rank before attempting to enlarge the subgroup:
const int n_aux_stuck = 20;

void saturator::reset_points(const vector<Point>& PP)
{
  Plist=PP;
  Plistx=PP;
  unsigned int i;
  for(i=0; i<Plistp.size(); i++) Plistx.push_back(Plistp[i]);
  rank=Plistx.size();
  TLimage=mat(0,rank); //  holds TL image in echelon form
  TLrank=0;
  qvar.init(); qvar++; qvar++;   // skip past 2 and 3
  stuck_counter=0;
  Eqptr=Eqlist.begin();
  newq=0;
}

int saturator::test_saturation(int pp, int ms)
{      
  p=pp;
  // We add a basis for the torsion/p to the given points:
  Plistx=Plist;
  Plistp = pCoTorsion(AllTorsion,p);
  int npcot = Plistp.size();
  if(npcot>0) 
    {
      if(verbose>1) 
	cout<< "saturator: adding "<<npcot<<" extra points before sieving: "
	    <<Plistp<<endl;
      int i;
      for(i=0; i<npcot; i++) Plistx.push_back(Plistp[i]);
    }
  rank=Plistx.size();
  TLimage=mat(0,rank); //  holds TL image in echelon form
  TLrank=0;
  
  if(use_div_pols)
    {
      pdivpol = makepdivpol(E,p);
      //cout<<p<<"-division poly = "<<pdivpol<<endl;
    }

  qvar.init(); qvar++; qvar++;   // skip past 2 and 3
  stuck_counter=0;
  log_index=0;
  Eqptr=Eqlist.begin();
  newq=0;
  while((TLrank<rank)&&(stuck_counter<ms)) nextq();
  return rank==TLrank;
}

int saturator::test_saturation_extra(int pp, int ms)
{
  stuck_counter=0;
  while((TLrank<rank)&&(stuck_counter<ms)) nextq();
  return rank==TLrank;
}

void saturator::nextq()
{
  int ntp=0;
  TLSS sieve; bigint q;
  while (ntp==0) /* (ntp<2) */
    {
      qvar++; q=qvar; 
      if(!qvar.ok()) 
	{
	  if(verbose>1)
	    cout<<"*** not enough precomputed primes for saturation (max = "
		<<maxprime()<<"), computing more primes..."<<flush;
	  the_primes.init(2*maxprime());
	  if(verbose>1)
	    cout<<"done, now max prime = "<<maxprime()<<endl;
	  qvar.init();
	  while(qvar.value()<=q) qvar++;
	  q=qvar; 
	  if(verbose>1)
	    cout<<"Continuing with q="<<q<<endl;
	}
      while(div(q,disc)) {qvar++; q=qvar; }
      if(verbose>1) cout<<"Trying q="<<q<<endl;
      if(newq||(Eqptr==Eqlist.end()))
	{
	  newq=1;
	  if(verbose>2) cout<<"Initializing q =  "<<q<<": "<<endl;
	  curvemodqbasis Eq(*E,q);  //,(p>10));
	  Eqlist.push_back(Eq);
	  sieve.assign(Eq);
	}
      else
	{
	  sieve.assign(*Eqptr);
	  Eqptr++;
	  newq=(Eqptr==Eqlist.end());
	  //	  cout<<"Using stored reduced curve mod "<<q<<": "<<Eq<<endl;
	}
      if(q==p) continue;
      if(use_div_pols) sieve.init(p,pdivpol,verbose);
      else sieve.init(p,verbose);
      ntp=sieve.get_rank();
    }

  if(verbose>1) cout<<"Using q = "<<q<<endl;
  mat TLim = sieve.map_points(Plistx);
  if(verbose>2) 
    {
      cout<<"Adding "<<ntp<<" rows to TL matrix;\n";
      cout<<TLim<<endl;
      cout<<"Now reducing to echelon form..."<<endl;
    }
  vec pcols, npcols; // not used
  long newTLrank, ny;
  mat newTLmat = echmodp(rowcat(TLimage,TLim),pcols, npcols, newTLrank, ny, p);
  if(verbose>2) 
    {
      cout<<"New rank = "<<newTLrank<<endl;
      cout<<"New TL matrix = "<<newTLmat<<endl;
    }
  if(newTLrank==TLrank)
    {
      stuck_counter++;
      if(verbose>1) 
	cout<<"Stuck at rank "<<TLrank<<" for the last "<<stuck_counter<<" primes"<<endl;
    }
  else 
    {
      stuck_counter=0;
      if(verbose>1) 
	cout<<"rank increases by "<<(newTLrank-TLrank)<<" to "<<newTLrank<<endl;
      TLimage=newTLmat;
      TLrank=newTLrank;
    }
  if(verbose>1) cout<<endl;
}

vec saturator::kernel_vector()
{
  if(TLrank==rank) return vec(0); // should not be called in this case!
  // Now we assume that TLimage is in echelon form
  mat ker = basis(pkernel(TLimage, p));
  return ker.col(1);
}

// enlarge basis if dim(kernel)>0:
int saturator::enlarge()
{
  if(TLrank==rank) return 0; // no enlargement;  should not be called in this case
  vec ker = basis(pkernel(TLimage, p)).col(1);
  if(verbose>0) cout<<"possible kernel vector = "<<ker<<endl;
  Point Q(E); int i, ci, keepi=-1;
  for(i=0; i<rank; i++) 
    {
      if((ci = mod(ker[i+1],p))) 
	{
	  if((keepi<0)&&(abs(ci)==1)) keepi=i;
	  Q+=ci*Plistx[i];
	}
    }
  if(verbose>0) cout<<"This point may be in "<<p<<"E(Q): "<<Q<<endl;
  vector<Point> Pi;
  long prec, original_prec;
  if(order(Q)==-1) // non-torsion point
    {
      // cout << "[attempting to divide by "<<p<<" using bit precision "
      //      <<bit_precision()<<"]"<<endl;
      Pi=division_points(*E,Q,p);
      if(Pi.size()==0)
        {
          prec = original_prec = bit_precision();
          //cout << "Saving bit precision "<<prec<<endl;
          prec *= 2;
          set_bit_precision(prec);
          // cout << "[attempting to divide by "<<p<<" using bit precision "
          //      <<prec<<"]"<<endl;
          Pi=division_points(*E,Q,p);
          set_bit_precision(original_prec);
          //cout << "Restoring bit precision "<<bit_precision()<<endl;
        }
    }
  if(Pi.size()==0)
    {
      if(verbose>0) cout<<"...but it isn't! "
			<<"(this may be due to insufficient precision: decimal precision "
                        <<prec<<" was used)"<<endl;
      return 0;
    }
  if(verbose>0) cout<<"...and it is! "<<endl;
  //  cout<<Pi<<endl;
  Q=Pi[0];
  if(verbose>0) cout<<"Replacing old generator #"<<(keepi+1)
		  <<" with new generator "<<Q<<endl;
  Plist[keepi]=Q;
  Plistx[keepi]=Q;
  log_index++;
  // reset TL matrix and q iteration
  TLimage=mat(0,rank); //  holds TL image in echelon form
  TLrank=0;
  qvar.init(); qvar++; qvar++;   // skip past 2 and 3
  stuck_counter=0;
  Eqptr=Eqlist.begin();
  newq=0;
  return 1;
}

// repeat testing saturation and enlarging until done:
// returns log_p of index
int saturator::do_saturation(int pp, int maxntries)
{
  p=pp;
  if(verbose>1) 
    cout<<"Testing "<<p<<"-saturation..."<<endl;
  if(test_saturation(p,n_aux_stuck)) return 0;
  if(verbose>1) 
    cout<<"Points not (yet) proved to be "<<p
	<<"-saturated, attempting enlargement..."<<endl;
  int n=0;
  while(1)
    {
      if(enlarge()) {n=0;}
      else
	{
	  if(verbose>1) cout<<" enlargement failed!"<<endl;
	  n++;
	  if(n==maxntries) // give up
	    {
	      cout<<"After "<<n<<" attempts at enlargement, giving up!\n";
	      cout<<"--points not proved "<<p<<"-saturated,"<<endl;
	      return -1;
	    }
	}
      if(test_saturation_extra(p,n_aux_stuck)) return log_index;
      if(verbose>1) cout<<"Points not (yet) proved to be "<<p
	  <<"-saturated, attempting enlargement..."<<endl;
    }
}

int saturator::do_saturation_upto(int maxp, int maxntries)
{
  int pi, p, index=1;
  primevar pvar;  p=pvar;
  while(p<=maxp)
    {
      if(verbose) cout<<"Checking "<<p<<"-saturation "<<endl;
      pi = do_saturation(p,maxntries);
      if(verbose&&(pi>=0))
	{
	  cout<<"Points have successfully been "<<p
	      <<"-saturated (max q used = "<<get_q()<<")"<<endl;
	  if(pi>0) cout<<"Index gain = "<<p<<"^"<<pi<<endl;
	}
      if(pi>0) while(pi--) index *= p;
      pvar++;
      p=pvar;
    }
  return index; 
}

int l2i(long i) {return (int)i;}
vector<int> lv2iv(const vector<long>& v)
{
  vector<int> ans;
  transform(v.begin(),v.end(),inserter(ans,ans.end()),ptr_fun(l2i));
  return ans;
}
int i2l(int i) {return (long)i;}
vector<long> iv2lv(const vector<int>& v)
{
  vector<long> ans;
  transform(v.begin(),v.end(),inserter(ans,ans.end()),ptr_fun(i2l));
  return ans;
}

int saturator::do_saturation(vector<long> plist, 
			     bigint& index, vector<long>& unsat, 
			     int maxntries)
{
  vector<int>iplist = lv2iv(plist), iunsat;
  int ans =  do_saturation(iplist,index,iunsat,maxntries);
  unsat = iv2lv(iunsat);
  return ans;
}

int saturator::do_saturation(vector<int> plist, 
			     bigint& index, vector<int>& unsat, 
			     int maxntries)
{
  unsigned int i; int pi, p;
  int success=1;
  index=1;
  if(verbose) cout<<"Checking saturation at "<<plist<<endl;
  for(i=0; i<plist.size(); i++)
    {
      p = plist[i];
      if(verbose) cout<<"Checking "<<p<<"-saturation "<<endl;
      pi = do_saturation(p,maxntries); // = log_index if >=0, -1 if failed
      if(pi<0) 
	{
	  cout<<p<<"-saturation failed!"<<endl;
	  unsat.push_back(p);
	  success=0;
	}
      else
	{
	  if(verbose)
	    {
	      if(pi>0) 
		{
		  cout<<"Points have successfully been "<<p
		      <<"-saturated (max q used = "<<get_q()<<")"<<endl;
		  cout<<"Index gain = "<<p<<"^"<<pi<<endl;
		}
	      if(pi==0) 
		{
		  cout<<"Points were proved "<<p
		      <<"-saturated (max q used = "<<get_q()<<")"<<endl;
		}
	    }
	  while(pi--) index *= p;
	}
    }
  return success;
}

int saturator::saturate(vector<long>& unsat, bigint& index, long sat_bd, 
			int egr, int maxntries, int odd_primes_only)
{
  // Determine the primes at which saturation is necessary: all those
  // up to index bound (but truncated at sat_bd unless this is -1),
  // and also the "Tamagawa primes" if the egr option is set

  vector<long> satprimes;
  primevar pr; 
  if(odd_primes_only) pr++;  // useful after a 2-descent
  int p=pr.value();
  bigint ib = index_bound(E,Plist,egr,(verbose>1));
  bigint ib0=ib;
  if(sat_bd==-1) sat_bd=SAT_MAX_PRIME;
  int bound_too_big = (ib>sat_bd);
  if(verbose)
    cout<<"Saturation index bound = "<<ib<<endl;
  if(bound_too_big)
    {
      if(!verbose) cout<<"Saturation index bound = "<<ib<<endl;
      cout<<"WARNING: saturation at primes p > "<<sat_bd
	  <<" will not be done;  \n"
	  <<"points may be unsaturated at primes between "<<sat_bd
	  <<" and index bound"<<endl;
      ib = sat_bd;
    }
  while(p<=ib)
    {
      satprimes.push_back(p);
      pr++; p=pr.value(); 
    }
  // In principle we should add these primes to unsat list, but in
  // practice we will not as there are likely to be too many!
#if(0)
  if(bound_too_big)
    while(p<=ib0)
      {
	unsat.push_back(p);
	pr++; p=pr.value(); 
      }
#endif
  if(egr)
    satprimes=vector_union(satprimes,tamagawa_primes(*E));

  // do the saturation:

  int sat_ok = do_saturation(satprimes, index, unsat);
  return (!bound_too_big) && sat_ok;
}

// This function returns a list of 0,1 or 2 points which generate 
// torsion modulo p*torsion:  
//
// 0 if p ndiv #torsion; else
// 1 (a generator) if torsion is cyclic; else
// 2 (a point of max order and an independent 2-torsion point)

vector<Point> pCoTorsion(const vector<Point>& AllTorsion, int p)
{
  long i, maxorder=0, ntorsion = AllTorsion.size();
  vector<Point> ans;

  // Case 0:

  if(ndivides(p,ntorsion)) return ans; // empty

  //  find point Q of maximal order:

  Point P,Q;
  for(i=0; (i<ntorsion)&&(maxorder<ntorsion); i++) 
    {
      P = AllTorsion[i];
      if(order(P)>maxorder) 
	{
	  Q=P; maxorder=order(Q);
	}
    }
  ans.push_back(Q);
  
  // Case 1: 

  if((maxorder==ntorsion)||(p>2)) return ans; // p-torsion is cyclic, return Q (generator)

  // Now order is 4, 8, 12 and torsion is not cyclic: add either
  // 2-torsion point not a multiple of Q:

  Q = (maxorder/2)*Q; // the 2-torsion point to avoid
  for(i=0; i<ntorsion; i++) 
    {
      P = AllTorsion[i];
      if((order(P)==2) && (Q != P)) 
	{
	  ans.push_back(P);     
	  return ans;
	}
    }
  return ans; // not necessary except to keep -Wall happy
}

int saturate_points(Curvedata& C, vector<Point>& points, 
		    bigint& index, vector<long>& unsat,  
		    long sat_bd, int egr, int verbose)
{
  saturator sieve(&C,verbose);
  sieve.set_points(points);
  int ans = sieve.saturate(unsat, index, sat_bd, egr, (verbose));
  points = sieve.getgens();
  return ans;
}

// Bound for the index of saturation for the given set of points. If
// egr is set it determines the egr subgroup of the group the points
// generate and only searches for points with egr, This might be faster
// in some cases...
//
// Flaw: the point search carried out in order to find a lower bound
// for the height of non-torsion points might find points which gain
// some index on the input points, but we do not use this.  Of course
// the caller can do their own points search first, in which case
// there is no (or less) loss, except that the searching has been done
// twice.  Some redesign would be needed to optimize this -- for
// example, index_bound could be part of the mw class.
//
// New version 08/08/06: uses class CurveHeightConst and strategy in
// ANTS7 paper to find a lower bound for the height of egr non-torsion
// points
//

bigint index_bound(Curvedata* C, vector<Point>& points, 
		   int egr, int verbose)
{
  if(verbose) 
    cout<<"Entering index_bound("<<(Curve)(*C)<<")"<<endl;
  int npts = points.size();

  bigfloat reg = regulator(points);
  if(verbose) 
    cout<<"Regulator of input points = "<<reg<<endl;

  bigfloat gamma=lattice_const(npts);
  if(verbose) 
    cout<<"Lattice constant = "<<gamma<<endl;

  // If egr==1, find regulator of egr subgroup

  bigfloat index = to_bigfloat(1), egr_reg=reg;
  if(egr)
    {
      index   = I2bigfloat(egr_index(points));
      egr_reg = index*index*reg;
      if(verbose)
	{
	  cout<<"Index of egr points = "<<index<<endl;
	  cout<<"Regulator of egr points  = "<<egr_reg<<endl;
	}
    }
  // else we'll divide lambda later instead
#ifdef USE_SEARCHING_STRATEGY
  bigfloat lambda=index_bound(C,points,egr,verbose);
  if(verbose) cout<<"lambda (via search) = "<<lambda<<endl;
#else // use ANTS7 strategy instead to get lower bound for egr height
  CurveHeightConst CHC(*C);
  CHC.compute();
  bigfloat lambda=CHC.get_value();
  if(verbose) cout<<"lambda (via ANTS7) = "<<lambda<<endl;
#endif
  if(!egr) 
    {
      bigfloat tam = I2bigfloat(Tamagawa_exponent(*C));
      lambda/=(tam*tam);
    }
  bigfloat ib = index*sqrt(reg*pow(gamma/lambda,npts));
  if(verbose) 
    cout<<"raw index bound = "<<ib <<endl;
  bigint ans = Ifloor(ib+0.1);  // be careful about rounding errors!
  if(ans<2) ans=1;  // In case 0.9999 has rounded down to 0
  if(verbose) 
    cout<<"Saturation index bound = "<<ans <<endl;
  return ans;
}  // end of index_bound()

// Tamagawa primes: primes dividing any Tamagawa number
vector<long> tamagawa_primes(const Curvedata& C)
{
  CurveRed CR(C);
  vector<bigint> badp = getbad_primes(CR);
  vector<long> tp;
  for(unsigned int i=0; i<badp.size(); i++)
    {
      tp = vector_union(tp,  pdivs(getc_p(CR,badp[i])));
    }
  return tp;
}


//end of file saturate.cc

#if(0)
  // Find optimally x-shifted curve for better point searching...

  bigint x_shift;
  Curvedata C_opt = opt_x_shift(*C,x_shift);
  int shift_flag = !is_zero(x_shift);
  if(shift_flag&&verbose) 
    cout<<"Using shifted model "<<(Curve)C_opt<<" for searching"<<endl;

  double hc; 
  if(egr)
    hc = egr_height_constant(C_opt);
  else  
    hc = height_constant(C_opt);
  if(verbose) 
    cout<<"height bound constant for shifted curve  = "<<hc<<endl;
  double hc1; 
  doublify(reg,hc1); 
  hc1 = (hc+hc1/(3.9));
  if(hc1>12) hc1=12;    // so hc1 = min(12,R/4+ht.const.)
  double hcx = hc1-hc; // = min(12-ht.const., R/4)
  if(hcx<0) {hcx=0.1; hc1=hcx+hc;}
  if(verbose)
    {
      if(egr)
	cout<<"Searching for egr points to naive height "<<hc1<<endl;
      else
	cout<<"Searching for all points to naive height "<<hc1<<endl;
    }
  if(hc1>max_search_bound) 
    {
      cout<<"\n***Warning: search bound of "<<hc1
	  <<" reduced to "<<max_search_bound
	  <<" -- points may not be saturated***"<<endl;     
      hc1=max_search_bound;
    }
  point_min_height_finder pmh(&C_opt,egr,verbose);
  pmh.search(to_bigfloat(hc1));
  bigfloat lambda=pmh.get_min_ht();
  newpoints = pmh.points();
  //  cout<<"Before shifting, newpoints = "<<newpoints<<endl;
  if(shift_flag) 
    for(unsigned int i=0; i<newpoints.size(); i++)
      newpoints[i] = transform(newpoints[i],C,BIGINT(1),
			   x_shift,BIGINT(0),BIGINT(0),1);
  //  cout<<"After shifting, newpoints = "<<newpoints<<endl;
  Point Pmin = pmh.get_min_ht_point();
  if(lambda==0)
    {	  
      lambda=hcx;
      if(verbose) 
	cout<<"No points found, lambda = "<<lambda<<endl;
    }
  else
    {
      if(verbose) 
	cout<<"Min height of points found = "<<lambda<<" (point "<<Pmin<<")"<<endl;
      if(lambda>hcx) lambda=hcx;
      if(verbose) 
	cout<<"Using lambda = "<<lambda<<endl;
    }
#endif
