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
 
#include <eclib/subspace.h>
#include <eclib/saturate.h>
#include <eclib/tlss.h>
#include <eclib/htconst.h>
#include <eclib/divpol.h>

// If point search bound is greater than this, output a warning
// message and reduce to this value:
const int max_search_bound = 18;

void saturator::reset_points(const vector<Point>& PP)
{
  Plist=PP;
  Plistx=PP;
  unsigned int i;
  for(i=0; i<Plistp.size(); i++) Plistx.push_back(Plistp[i]);
  rank=Plistx.size();
  TLimage=mat_l(0,rank); //  holds TL image in echelon form
  TLrank=0;
  qvar.init(); qvar++; qvar++;   // skip past 2 and 3
  stuck_counter=0;
  the_index_bound = BIGINT(0);
}

// initialize index bound
void saturator::set_index_bound()
{
  the_index_bound = index_bound(Plist, egr_flag, (verbose>1));
}

// return current index bound (compute if necessary)
bigint saturator::get_index_bound()
{
  if (is_zero(the_index_bound))
    set_index_bound();
  return the_index_bound;
}

// test whether p is less than the saturation index, or a Tamagawa prime
int saturator::trivially_saturated(long p)
{
  return ((p>the_index_bound)
          &&
          (find(tam_primes.begin(), tam_primes.end(), p) == tam_primes.end()));
}

int saturator::test_saturation(int pp, int ms)
{
  p=pp;
  if (trivially_saturated(p))
    return 1; // success
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
  TLimage=mat_l(0,rank); //  holds TL image in echelon form
  TLrank=0;
  
  if(use_div_pols)
    {
      pdivpol = division_polynomial(E,p);
      //cout<<p<<"-division poly = "<<pdivpol<<endl;
    }

  qvar.init(); qvar++; qvar++;   // skip past 2 and 3
  stuck_counter=0;
  log_index=0;
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
      if(q==p) continue;

      if(verbose>2) cout<<"Trying q="<<q<<endl;

      // First just check the order of E mod q, skip this q if not a multiple of p

      map<bigint,bigint>::iterator Eqoi = Emodq_order.find(q);
      bigint order_mod_q;
      if(Eqoi==Emodq_order.end())
	{
	  if(verbose>2) cout<<"Computing order mod q =  "<<q<<": "<<endl;
          curvemodq Eq(*E,q);
          if(0) // use orders of some random points as a proxy
            {
              bigint upper, lower; // bounds on group order
              set_hasse_bounds(q,lower,upper);
              pointmodq P1 = Eq.random_point();
              bigint n1 = my_order_point(P1,lower,upper);
              if (verbose>2)
                cout<<"q="<<q<<"\tn1 = "<<n1<<endl;
              order_mod_q = n1;
            }
          else
            {
              order_mod_q = Eq.group_order();
            }

          Emodq_order[q] = order_mod_q;
          if (verbose>2)
            cout<<"Setting order mod "<<q<<" to "<<order_mod_q<<endl;
	}
      else
	{
          order_mod_q = Eqoi->second;
          if (verbose>2)
            cout<<"reusing order mod "<<q<<" as "<<order_mod_q<<endl;
	}
      if ((order_mod_q%p)!=0)
        {
          if (verbose>2)
            cout<<"Order mod "<<q<<" is "<<order_mod_q<<", not a multiple of p="<<p<<endl;
          continue;
        }
      else
        {
          if (verbose>1)
            cout<<"*** using q="<<q<<" with order "<<order_mod_q<<" a multiple of p="<<p<<endl;
        }
      // next compute the structure of q if not yet known

      map<bigint,curvemodqbasis>::iterator Eqi = Emodq.find(q);
      if(Eqi==Emodq.end())
	{
	  if(verbose>2) cout<<"Initializing q =  "<<q<<": "<<endl;
	  curvemodqbasis Eq(*E,q);  //,(p>10));
	  Emodq[q] = Eq;
	  sieve.assign(Eq);
          q_tally[q] = 0;
	}
      else
	{
	  sieve.assign(Eqi->second);
	  //	  cout<<"Using stored reduced curve mod "<<q<<endl;
	}
      if(use_div_pols) sieve.init(p,pdivpol,verbose);
      else sieve.init(p,verbose);
      ntp=sieve.get_rank();
    }

  if(verbose>1) cout<<"Using q = "<<q<<endl;
  q_tally[q] += 1;
  if (q>maxq)
    {
      maxq=q;
      maxp=p;
    }
  mat_l TLim = sieve.map_points(Plistx);
  if(verbose>2) 
    {
      cout<<"Adding "<<ntp<<" rows to TL matrix;\n";
      cout<<TLim<<endl;
      cout<<"Now reducing to echelon form..."<<endl;
    }
  vec_l pcols, npcols; // not used
  long newTLrank, ny;
  mat_l newTLmat = echmodp(rowcat(TLimage,TLim),pcols, npcols, newTLrank, ny, p);
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

vec_l saturator::kernel_vector()
{
  if(TLrank==rank) return vec_l(0); // should not be called in this case!
  // Now we assume that TLimage is in echelon form
  mat_l ker = basis(pkernel(TLimage, p));
  return ker.col(1);
}

// enlarge basis if dim(kernel)>0:
int saturator::enlarge()
{
  if(TLrank==rank) return 0; // no enlargement;  should not be called in this case
  vec_l ker = basis(pkernel(TLimage, p)).col(1);
  if(verbose>0) cout<<"possible kernel vector = "<<ker<<endl;
  Point Q(E), newQ(E); int flag, i, ci, keepi=-1;
  for(i=0; i<rank; i++) 
    {
      if((ci = mod(ker[i+1],p))) 
	{
	  if((keepi<0)&&(abs(ci)==1)) keepi=i;
	  Q+=ci*Plistx[i];
	}
    }
  if(verbose>0) cout<<"This point may be in "<<p<<"E(Q): "<<Q<<endl;

  flag = !Q.is_torsion();
  if (flag)
    {
      // this used elog method:
      //divide_point(*E, Q, p, newQ);
      // this uses division polynomials (exact):
      vector<Point> newQlist = Q.division_points(p);
      if (newQlist.size()>0)
        {
          newQ = newQlist[0];
        }
      else
        {
          flag = 0;
        }
    }

  if(!flag)
    {
      if(verbose>0)
        {
          cout<<"...but it isn't! (this may be due to insufficient precision)";
        }
      return 0;
    }
  if(verbose>0) cout<<"...and it is! "<<endl;
  if(verbose>0) cout<<"Replacing old generator #"<<(keepi+1)
		  <<" with new generator "<<newQ<<endl;
  Plist[keepi]=newQ;
  Plistx[keepi]=newQ;
  log_index++;

  // Now the points we have have a regulator which is reduced by a
  // factor of p^2, the saturation index bound can often be reduced by
  // a factor of p, but not necessarily so.  For simplicity we just
  // recompute the index bound.

  bigint old_index_bound = the_index_bound;
  set_index_bound();
  if(verbose)
    if (the_index_bound < old_index_bound)
      cout << "Reducing index bound from " << old_index_bound <<" to " << the_index_bound << endl;
    else
      cout << "The index bound " << the_index_bound << " has not changed"<<endl;

  // reset TL matrix and q iteration
  TLimage=mat_l(0,rank); //  holds TL image in echelon form
  TLrank=0;
  qvar.init(); qvar++; qvar++;   // skip past 2 and 3
  stuck_counter=0;
  return 1;
}

// repeat testing saturation and enlarging until done:
// returns log_p of index
int saturator::do_saturation(int pp, int maxntries)
{
  p=pp;
  if(verbose>1)
    cout<<"Testing "<<p<<"-saturation..."<<endl;
  if (trivially_saturated(p))
    return 0;  // index=1, log=0
  if(test_saturation(p, maxntries))
    return 0;
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
      if(test_saturation_extra(p, maxntries)) return log_index;
      if(verbose>1) cout<<"Points not (yet) proved to be "<<p
	  <<"-saturated, attempting enlargement..."<<endl;
    }
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
			     long& index, vector<long>& unsat, 
			     int maxntries)
{
  vector<int>iplist = lv2iv(plist), iunsat;
  int ans =  do_saturation(iplist,index,iunsat,maxntries);
  unsat = iv2lv(iunsat);
  return ans;
}

int saturator::do_saturation(vector<int> plist, 
			     long& index, vector<int>& unsat, 
			     int maxntries)
{
  unsigned int i; int pi, p;
  int success=1;
  index=1;
  if(verbose) cout<<"Checking saturation at "<<plist<<endl;
  for(i=0; i<plist.size(); i++)
    {
      p = plist[i];
      if (trivially_saturated(p))
        continue; // to the next prime in the list
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

int saturator::saturate(vector<long>& unsat, long& index,
                        long sat_bd, long sat_low_bd,
			int maxntries)
{
  // Determine the primes at which saturation is necessary: all those
  // up to index bound (but truncated at sat_bd unless this is -1),
  // and also the "Tamagawa primes" if the egr option is set

  vector<long> satprimes;
  primevar pr;
  while(pr.value()<sat_low_bd) pr++;
  int p=pr.value();

  bigint ib = get_index_bound();
  if(verbose)
    {
      cout<<"Saturation index bound ";
      if (egr_flag) cout << "(for points of good reduction) ";
      cout<< " = "<<ib<<endl;
    }

  if (sat_bd==-1) // no bound was specified, but we warn if the computed bound is large
    {
      if ((ib>SAT_MAX_PRIME) && verbose)
        {
          cout<<"Saturation index bound = "<<ib<<" is large, ";
          cout<<"and saturation will take a long time."<<endl;
        }
    }
  else // a bound was specified
    {
      if (ib<sat_bd) // we can reduce the specified bound as it was larger than necessary
        {
          if (verbose)
            {
              cout << "Reducing saturation bound from given value " << sat_bd;
              cout << " to computed index bound " << ib << endl;
            }
        }
      else // we'll use the specified bound, so cannot guarantee saturation
        {
          if (verbose)
            {
              cout << "Only p-saturating for p up to given value " << sat_bd << ".\n";
              cout << "The resulting points may not be p-saturated for p between this ";
              cout << "and the computed index bound " << ib << endl;
            }
          ib = sat_bd;
        }
    }
  while(p<=ib)
    {
      //cout<<"adding p="<<p<<" to saturation list"<<endl;
      satprimes.push_back(p);
      pr++; p=pr.value();
    }

  // under the egr option, the computed bound is on the saturation
  // index of the egr points (of everywhere good reduction) and we
  // must also be sure to saturate at any primes dividing any Tamagawa
  // number.  If an upper bound sat_bd has been specified manually we
  // do not add these primes (some may be under the given bound
  // anyway).  Before this step, satprimes contains all primes between
  // sat_low_bd and ib, so we just add any Tamagawa primes greater
  // than ib.

  if(egr_flag)
    {
      if (verbose)
        cout << "Tamagawa index primes are " << tam_primes << endl;
      for (vector<long>::iterator pi = tam_primes.begin(); pi!=tam_primes.end(); pi++)
        {
          p = *pi;
          if ((p > ib) && ((sat_bd==-1) || (p <= sat_bd)))
            {
              if (verbose)
                cout << "adding Tamagawa index prime " << p << " to saturation list" << endl;
              satprimes.push_back(p);
            }
        }
    }

  // do the saturation.  Will return ok iff we succeeded in saturating
  // at all p in satprimes, otherwise the failures will be in unsat.

  int sat_ok = do_saturation(satprimes, index, unsat, maxntries);
  return sat_ok;
}

void saturator::show_q_tally()
{
  cout << "Summary of auxiliary primes used" <<endl;
  int num_q_used = 0;
  map<bigint,int>::iterator qcount;
  cout << "Number of q used: " << q_tally.size() << endl;
  cout << "Maximum   q used: " << maxq << " (used for p="<<maxp<<")"<<endl;
  if (verbose<2)
    return;
  cout << "Counts of how many times each q was used:" << endl;
  bigint q;
  int c;
  for (qcount = q_tally.begin(); qcount!=q_tally.end(); qcount++)
    {
      q = qcount->first;
      c = qcount->second;
      if (c)
        cout << q << "\t" << c <<endl;
    }
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
		    long& index, vector<long>& unsat,
		    long sat_bd, long sat_low_bd,
                    int egr, int verbose)
{
  saturator sieve(&C, egr, verbose);
  sieve.set_points(points);
  int ans = sieve.saturate(unsat, index, sat_bd, sat_low_bd);
  points = sieve.getgens();
  if (verbose>0)
    sieve.show_q_tally();
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

bigint index_bound(vector<Point>& points,
		   int egr, int verbose)
{
  int npts = points.size();
  if (npts==0)
    return BIGINT(1);

  Curvedata C = Curvedata(points[0].getcurve(), 0);
  if(verbose)
    cout<<"Entering index_bound("<<(Curve)(C)<<", egr="<<egr<<")"<<endl;

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
  bigfloat lambda = lower_height_bound(C, egr);
  if(verbose) cout<<"lambda (via ANTS7) = "<<lambda<<endl;
#endif
  // no need to use egr_reg in next line since we multiply by index
  bigfloat ib = index*sqrt(reg*pow(gamma/lambda,npts));
  bigint ans = Ifloor(ib+0.1);  // be careful about rounding errors!
  if(ans<2) ans=BIGINT(1);  // In case 0.9999 has rounded down to 0
  if(verbose)
    cout<<"Saturation index bound " << ib << ", rounds down to "<<ans<<endl;
  return ans;
}  // end of index_bound()

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
