// sieve_search.cc: implementations of classes point_processor and qsieve
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
 
// File        : sieve_search.c
// Author      : John Cremona, Sophie Labour
//		 Adaption of M. Stoll's and J. Cremona's code
// Last change : SL, Apr 21 1999, first definitions
// Last change : JC, Dec 2 1999, adapted for direct use in mwrank et al.
// Last change : JC, Aug 19 2002, adapted for stdc++ library

#include <stdlib.h> // for qsort
#include "interface.h"
#include "marith.h"
#include "sieve_search.h"

//debug info printed out if DEBUG_QS is defined to 0; more info if defined to 1
//#define DEBUG_FORBIDDEN
//#define DEBUG_QS 0
//#define DEBUG_QS 1

int compare_entries(const void *a, const void *b)
{
  double diff = (((entry *)a)->r - ((entry *)b)->r);
  return (diff > 0) ? 1 : (diff < 0) ? -1 : 0;
}


// ***********************************************************
//   Class qsieve
// ***********************************************************

//char qsieve::init_made = 0;
long qsieve::prime[QS_NUM_PRIMES] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251};

double qsieve::ratio1_def = 1000;
double qsieve::ratios2[QS_MAX_DEGREE+1] = { 0, 2.5, 4, 5, 6.5, 7.5, 9, 10.5, 12, 14, 15.5 };
double qsieve::ratio1 = 0.0;
double qsieve::ratio2 = 0.0;

long qsieve::tech = QS_DEFAULT_TECH;

void qsieve::init_data()
{
  init_made=0;
  bits= new bit_array[QS_LONG_LENGTH];
  squares= new char*[QS_NUM_PRIMES];
  is_f_square= new char*[QS_NUM_PRIMES];
  pnn= new long[QS_NUM_PRIMES];
  kpa= new long[QS_NUM_PRIMES];
  check_denom=1;
  use_squares=0;
  odd_nums=0;
  num_surv1=0;
  num_surv2=0;
  long i;
  for (i=0;i<QS_NUM_PRIMES;i++)
    {
      squares[i]=new char[QS_MAX_PRIME];
      is_f_square[i]=new char[QS_MAX_PRIME];
    }
  init_all();//initialize bits[] and squares[][]
  // Set global variables to their default values 
  num_inter = 0;       // No interval up to now ([-hlim,hlim]
  sieve_primes1 = -1;  // automatic determination of this number 
  sieve_primes2 = -1;  // automatic determination of this number 
  no_check = 0;        // do the check by default 
  points_at_infty = 1; // also find points at infinity 
  b_low = 1;           // denominators go from 1 to h 
  b_high = height;
  use_opt = 1;         // use even denominator optimisation 
  array_size = QS_DEFAULT_SIZE;
}

qsieve::qsieve(point_processor* acurve, int deg, vector<bigint> coeff, int verb)
: curve(acurve), degree(deg), verbose(verb)
{
  long i;
  for (i=0;i<=degree;i++)
   c[i]=coeff[i];
#ifdef DEBUG_QS
  cout<<"In qsieve constructor, degree = "<<degree<<"\n";
  cout<<"coeff = [";
  for (i=0;i<=degree;i++) cout<<c[i]<<" ";
  cout<<"]"<<endl;
#endif  

  set_height(QS_DEFAULT_HEIGHT);
  init_data();
}

qsieve::qsieve(point_processor* acurve, int deg, vector<bigint> coeff, bigfloat h_limx, int verb)
  : curve(acurve), degree(deg), verbose(verb)
{
  double h_lim; doublify(h_limx,h_lim);
  long i;
  for (i=0;i<=degree;i++)
    c[i]=coeff[i];
#ifdef DEBUG_QS
  cout<<"In qsieve constructor, degree = "<<degree<<"\n";
  cout<<"coeff = [";
  for (i=0;i<=degree;i++) cout<<c[i]<<" ";
  cout<<"]"<<endl;
#endif  

  set_height(h_lim);
  init_data();
}
/*
qsieve::qsieve(point_processor* acurve, int deg, vector<bigint> coeff, double h_lim, double up, double low, int verb)
   : verbose(verb), curve(acurve), degree(deg)
{
  long i;
  for (i=0;i<=degree;i++)
    c[i]=coeff[i];
#ifdef DEBUG_QS
  cout<<"In qsieve constructor, degree = "<<degree<<"\n";
  cout<<"coeff = [";
  for (i=0;i<=degree;i++) cout<<c[i]<<" ";
  cout<<"]"<<endl;
#endif  

  init_all();
  height=(long)(floor(exp(h_lim)));
  w_height=(height-1)/QS_LONG_LENGTH + 1;
  upper_p=lower_p=1;
  upper=up;
  lower=low;
  w_upper=upper/QS_LONG_LENGTH;
  w_lower=lower/QS_LONG_LENGTH;
}*/
      
void qsieve::init_all()
{
  if (!init_made)
    {
      init_made=1;
      // intialize bits[] 
      bit_array bit = (bit_array)1;
      long n;
      for(n = 0; n < (long)QS_LONG_LENGTH; n++) { bits[n] = bit; bit <<= 1; }
      
      // initialize squares[][] 
      //      if (verbose)
      //	cout<<"Init squares"<<endl;
      long a, pn;
      long p;
      for(pn = 0; pn < QS_NUM_PRIMES; pn++)
	{
	  p = prime[pn];
	  for(a = 0 ; a < p ; a++) squares[pn][a] = 0;
	  for(a = 0 ; a < p; a += 2) squares[pn][(a*a)%p] = 1;
	  // if p==2 then this does not work! 
	}
    }
#ifdef DEBUG_QS
    cout<<"QS_LONG_LENGTH "<<QS_LONG_LENGTH<<"\nQS_LONG_SHIFT "<<QS_LONG_SHIFT<<endl;
#endif
}

qsieve::~qsieve()
{
  long i;
  delete[] bits;
  delete[] pnn;
  delete[] kpa;
  for (i=0; i<QS_NUM_PRIMES; i++)
    {
      delete[] squares[i];
      delete[] is_f_square[i];
    }
  delete[] squares;
  delete[] is_f_square;
}

void qsieve::set_intervals(vector<double> interv,int nb_bnd,int start_low,int pos_x_only)
{
#ifdef DEBUG_QS
  cout<<"Entering set_intervals with nb_bnd="<<nb_bnd<<", height = "<<height<<endl;
#endif

  long j=0;
  int uplow=0;
  num_inter=0; //reset
  double the_min;
  if (pos_x_only)
    the_min=0;
  else
    the_min=-height;
  
  while (j<nb_bnd && interv[j]<the_min) {j++;}
  
  if (j==nb_bnd)
    {
      domain[0].low=the_min;
      domain[0].up=height;
      num_inter=1;
      return;
    }

  if (start_low)
    {
      if (j&1)//upper bound
	{
	  domain[0].low=the_min;
	  domain[0].up=interv[j];
	  num_inter++;
	  uplow++; uplow++;
	}
      else //lower bound
	{
	  domain[0].low=interv[j];
	  uplow++;
	}
    }
  else
    {
      if (j&1)//lower bound
	{
	  domain[0].low=interv[j]; 
	  uplow++;
	}
      else //upper bound
	{
	  domain[0].low=the_min;
	  domain[0].up=interv[j];
	  num_inter++;
	  uplow++; uplow++;
	}
    }
  j++;

  while (j<nb_bnd)
    {
      if (uplow&1) //upper bound
	{
	  if (interv[j]<domain[num_inter].low)
	    cout<<"qsieve::set_intervals:interv[i]>interv[i+1]"<<endl;
	  if (interv[j]>=height)
	    j=nb_bnd; //done
	  else
	    {
	      domain[num_inter].up=interv[j];
	      num_inter++; uplow++; j++; //goto next lower bound
	    }
	}
      else  //lower bound
	{
	  if ((num_inter>0) && (interv[j]<domain[num_inter-1].up))
	    cout<<"qsieve::set_intervals:interv[i]>interv[i+1]"<<endl;
	  else
	    {
	      if (interv[j]>=height)
		j=nb_bnd; //done
	      else
		{
		  domain[num_inter].low=interv[j];
		  uplow++; j++; //goto next upper bound
		}
	    }
	}
    }
  
  if (uplow&1) //finish last interval
    {
      domain[num_inter].up=height;
      num_inter++;
    }

#ifdef DEBUG_QS
  cout<<"Number of intervals:"<<num_inter<<endl;
  long k;
  for(k=0;k<num_inter;k++)
    cout<<"["<<domain[k].low<<","<<domain[k].up<<"] ";
  cout<<endl;
#endif

}

void qsieve::set_sieve_primes1(long sp1)
{
  sieve_primes1=sp1;
  if (sieve_primes1 < 0) 
    sieve_primes1 = 0;
  else
    { 
      if (sieve_primes2 >= 0 && sieve_primes1 > sieve_primes2)
	sieve_primes1 = sieve_primes2;
      else 
	{ 
	  if(sieve_primes1 > QS_NUM_PRIMES) 
	    sieve_primes1 = QS_NUM_PRIMES; 
	}
    }
}

void qsieve::set_sieve_primes2(long sp2)
{
  sieve_primes2=sp2;
  if(sieve_primes2 < 0) 
    sieve_primes2 = 0;
  else
    { 
      if(sieve_primes1 >= 0 && sieve_primes1 > sieve_primes2)
	sieve_primes2 = sieve_primes1;
      else 
	{ 
	  if(sieve_primes2 > QS_NUM_PRIMES) 
	    sieve_primes2 = QS_NUM_PRIMES; 
	}
    }
}


void qsieve::init_f()
{
//    if (verbose)
//      cout<<"Entering init_f"<<endl;

  //init is_f_square
  {
    long n, fdc = 0;
    forbidden[fdc]=0;
    if(!(degree&1))
      {// check if we can exclude even numerators
	unsigned long t1, t2;
	t1 = posmod(c[0],8);
	t2 = posmod(c[1],4) + 2*posmod(c[2],2);
	switch(t1)
	  { 
	  case 2: case 3: case 6:
	  case 7: if((t2&1) == 0) odd_nums = 1; break;
	  case 5: if((t2&3) == 0) odd_nums = 1; break;
	  } 
      }
    
    if(check_denom)
      {//check if we can exclude denominators divisible by 
       //low powers of 2,3,5, or 7 
	unsigned long t, t1, t2, t3;
	//check powers of 2 first 
	t1 = posmod(c[degree],8);
	t2 = posmod(c[degree-1],4);
	t3 = posmod(c[degree-2],2);
	switch(t1)
	  { case 2: case 3: case 6: case 7:
	    if((t2&1) == 0)
	      { forbidden[fdc] = 2; fdc++; use_opt = 0; }
	    else
	      { forbidden[fdc] = 4; fdc++; }
	      break;
	  case 5:
	    if(((t2 + 2*t3)&3) == 0)
	      { forbidden[fdc] = 2; fdc++; use_opt = 0; }
	    else
	      { if((t2&1) == 0)
		{ forbidden[fdc] = 4; fdc++; }
	      else
		{ forbidden[fdc] = 8; fdc++; }
	      }
	    break;
	  }
#ifdef DEBUG_FORBIDDEN
	int i;
	cout<<"Forbidden divisors of the denominator (2): fdc="<<fdc<<endl;
	for (i=0;forbidden[i];i++)
	  cout<<forbidden[i]<<", ";
	cout<<endl;
#endif
	//check 3, 5, 7 
	t = posmod(c[degree],9);
	t1 = t%3;
	t2 = posmod(c[degree-1],3);
	if(t1 == 2)
	  { forbidden[fdc] = 3; fdc++; }
	else
	  { if(t == 3 || t == 6)
	    { forbidden[fdc] = (t2 == 0) ? 3 : 9; fdc++; }
	  } 
#ifdef DEBUG_FORBIDDEN
	cout<<"Forbidden divisors of the denominator (3): fdc="<<fdc<<endl;
	for (i=0;forbidden[i];i++)
	  cout<<forbidden[i]<<", ";
	cout<<endl;
#endif
        
	t = posmod(c[degree],25);
	t1 = t%5;
	t2 = posmod(c[degree-1],5);
	if((t1 == 2) || (t1 == 3))
	  { forbidden[fdc] = 5; fdc++; }
	else
	  { if(t != 0 && t1 == 0)
	    { forbidden[fdc] = (t2 == 0) ? 5 : 25; fdc++; }
	  }
#ifdef DEBUG_FORBIDDEN
	cout<<"Forbidden divisors of the denominator (5): fdc="<<fdc<<endl;
	for (i=0;forbidden[i];i++)
	  cout<<forbidden[i]<<", ";
	cout<<endl;
#endif
	
	t = posmod(c[degree],49);
	t1 = t%7;
	t2 = posmod(c[degree-1],7);
	if((t1 == 3) || (t1 == 5) || (t1 == 6))
	  { forbidden[fdc] = 7; fdc++; 
#ifdef DEBUG_FORBIDDEN
	    cout<<"Adding 7 to forbidden list, now fdc="<<fdc<<endl;
#endif
	  }
	else
	  { if((t != 0) && (t1 == 0))
	    { forbidden[fdc] = ((t2 == 0) ? 7 : 49); fdc++; 
#ifdef DEBUG_FORBIDDEN
	    cout<<"Adding 7 or 49 to forbidden list, now fdc="<<fdc<<endl;
#endif
	    }
	  }
      }
#ifdef DEBUG_FORBIDDEN
    int i;
	cout<<"Forbidden divisors of the denominator (7): fdc="<<fdc<<endl;
	for (i=0;forbidden[i];i++)
	  cout<<forbidden[i]<<", ";
	cout<<endl;
#endif
    
    long pn;  
    for (pn=0;pn<QS_NUM_PRIMES; pn++)
      {
	long p=prime[pn];
	
	//compute coefficients mod p
	long n, t;
	for(n=0; n<=degree; n++) 
	  {
	    t=posmod(c[n],p);
	    coeffs_mod_p[pn][n]=t;
	  } 
	long np=squares[pn][coeffs_mod_p[pn][0]];
	is_f_square[pn][0]=np;
	long a;
	for(a = 1 ; a < p; a++) 
	  {
	    //compute f(a) mod p with Horner method
	    long s = coeffs_mod_p[pn][degree];
	    for(n = degree - 1 ; n >= 0 ; n--)
	      {
		s *= a;
		s += coeffs_mod_p[pn][n];
		s %= p;
	      }
	    if((is_f_square[pn][a]=squares[pn][s])) 
	      np++; //np is the number of a / f(a) is a square mod p
	  }
	// Fill array with info for p 
	prec[pn].p = p;
	prec[pn].n = pn;
	if(degree&1 || squares[pn][coeffs_mod_p[pn][degree]])
	  { prec[pn].r = ((double)(np*(p-1) + p))/((double)(p*p)); }
	else
	  { prec[pn].r = (double)np/(double)p;
	  /* denominator divisible by p is excluded */
	  if(check_denom && (p > 7)) { forbidden[fdc] = p; fdc++; } 
	  }
      }
#ifdef DEBUG_FORBIDDEN
	cout<<"Forbidden divisors of the denominator (>7): fdc="<<fdc<<endl;
	for (i=0;forbidden[i];i++)
	  cout<<forbidden[i]<<", ";
	cout<<endl;
#endif

    // sort the array to get at the best primes 
    qsort(prec, QS_NUM_PRIMES, sizeof(entry), compare_entries);
    //Determine the number of sieving primes:
    //First stage: take smallest n such that prod(j<=n) prob(j) * ratio1 < 1.
    //Second stage: take smallest N such that (1-prob(N+1)) * ratio2 < 1. 
    if(sieve_primes1 < 0)
      { long n = 0, m = (sieve_primes2 >= 0) ? sieve_primes2-1 : QS_NUM_PRIMES-1;
      double prod;
      if(ratio1 == 0.0) ratio1 = ratio1_def;
      prod = ratio1;
      for(n = 0; n < m; n++)
	{ prod *= prec[n].r;
	if(prod < 1.0) break;
	}
      sieve_primes1 = n + 1;
      }
    if(sieve_primes2 < 0)
      { long n;
      if(ratio2 == 0.0) ratio2 = ratios2[degree];
      for(n = sieve_primes1; n < QS_NUM_PRIMES; n++)
	if(ratio2*(1.0 - prec[n].r) < 1.0) break;
      sieve_primes2 = n;
      }
    for(n = 0; n < sieve_primes2; n++)
      { pnn[n] = prec[n].n;
      sieves[n].p = prime[pnn[n]];
      }
    if (verbose)
      {
	cout<<"Using speed ratios "<<ratio1<<" and "<<ratio2<<endl;
	if (odd_nums)
	  cout<<"Even numerators are excluded."<<endl;
	cout<<sieve_primes1<<" primes used for first stage of sieving\n"
	  <<sieve_primes2<<" primes used for both stages of sieving together."<<endl;
	if (sieve_primes2>0)
	  {
	    cout<<"Sieving primes:\n First stage: ";
	    long i;
	    for(i = 0; i < sieve_primes1; i++)
	      cout<<prime[pnn[i]]<<", ";
	    cout<<endl<<" Second stage: ";
	    for(;i<sieve_primes2;i++)
	      cout<<prime[pnn[i]]<<", ";
	    cout<<endl;
	    cout<<"Probabilities: Min("<<prec[0].p<<") = "<<prec[0].r
	      <<", Cut1("<<prec[sieve_primes1-1].p<<") = "
		<<prec[sieve_primes1-1].r;
	    cout<<", Cut2("<<prec[sieve_primes2-1].p<<") = "
	      <<prec[sieve_primes2-1].r<<", Max("
		<<prec[QS_NUM_PRIMES-1].p<<") = "
		  <<prec[QS_NUM_PRIMES-1].r<<endl<<endl;
	  }
      }
    //terminate array of forbidden divisors 
    if(check_denom)
      { forbidden[fdc] = 0;
      if(fdc == 0) check_denom = 0;
      }
    if (verbose)
      {
	long i;
	cout<<"Forbidden divisors of the denominator:"<<endl;
	for (i=0;forbidden[i];i++)
	  cout<<forbidden[i]<<", ";
	cout<<endl<<endl;
      }
    
  }//end of init is_f_square


  //init sieve[][][]
  {
    //help is an array for temporarily storing the sieving information. Bit j of help[b][a0] says whether f(a/b) is a square mod p, where a = a0*QS_LONG_LENGTH + j. a runs from 0 to k*p - 1, where k*p is the least multiple of p exceeding QS_LONG_LENGTH.
    bit_array help[QS_MAX_PRIME][QS_MAX_PRIME / QS_LONG_LENGTH + 2];
    //help2 plays the same role when using the optimisation for even denoms 
    bit_array help2[QS_MAX_PRIME][QS_MAX_PRIME / QS_LONG_LENGTH + 2];

    sieve= new bit_array**[sieve_primes2];
    sieve2= new bit_array**[sieve_primes2];
    
    //  cout<<"Entering init_sieve"<<endl;
    long a, b, i, pn; 
    long p, kp, bb, aa;
    for(pn = 0; pn < sieve_primes2; pn++)
      {
	p = prime[pnn[pn]];
	sieve[pn]=new bit_array*[p+1];
	sieve2[pn]=new bit_array*[p+1];
	//kp = k*p such that k = 1 if p > tech and k*p maximal <= tech otherwise 
	if(p > tech)
	  kp = p;
	else
	  kp = (tech/p) * p;
	kpa[pn] = kp;
	//determine which quotients a/b are excluded mod p
	bb = (b_high >= p) ? p : b_high + 1;
	for (b=0;b<p+1;b++)
	  {
	    sieve[pn][b]=new bit_array[p+1];
	    sieve2[pn][b]= new bit_array[p+1];
	  }

	aa = p>>QS_LONG_SHIFT;
	if(aa == 0) aa = 1;
	if(!odd_nums)
	  {
	    //initialize help[][]
	    for(b = 0; b < bb; b++)
	      for(a = 0; a < aa; a++) 
		help[b][a] = bit_zero;
	    if(p < (long)QS_LONG_LENGTH) //small p
	      for(a = 0; a < p; a++)
		{
		  if(is_f_square[pnn[pn]][a])
		    for(b = 1; b < bb; b++)
		      { 
			long ab = (a * b) % p;
			long i;
			for(i = 0; i <= (long)QS_LONG_LENGTH/p; i++) //repeat the pattern
			  {
			    help[b][ab>>QS_LONG_SHIFT] |= bits[ab & QS_LONG_MASK];
			    ab += p;
			  } 
		      }
		}
	    else
	      for(a = 0; a < p; a++)
		{
		  if(is_f_square[pnn[pn]][a])
		    for(b = 1; b < bb; b++)
		      { 
			long ab = (a * b) % p;
			help[b][ab>>QS_LONG_SHIFT] |= bits[ab & QS_LONG_MASK];
		      }
		}
	  }
	
	if(odd_nums || use_opt)
	  { // prepare help2 
	    for(b = 0; b < bb; b++)
	      for(a = 0; a <= aa; a++) help2[b][a] = bit_zero;
	    if(p < (long)QS_LONG_LENGTH)
	      for(a = 0; a < p; a++)
		{ 
		  if(is_f_square[pnn[pn]][a])
		    for(b = 1; b < bb; b++)
		      { 
			long ab = (a * b) % p;
			ab = (ab&1) ? (ab-1)>>1 : (ab-1+p)>>1;
			for(i = 0; i <= (long)QS_LONG_LENGTH / p; i++)
			  {
			    help2[b][ab>>QS_LONG_SHIFT] |= bits[ab & QS_LONG_MASK];
			    ab += p;
			  } 
		      }
		}
	    else
	      for(a = 0; a < p; a++)
		{ 
		  if(is_f_square[pnn[pn]][a])
		    for(b = 1; b < bb; b++)
		      { 
			long ab = (a * b) % p;
			ab = (ab&1) ? (ab-1)>>1 : (ab-1+p)>>1;
			help2[b][ab>>QS_LONG_SHIFT] |= bits[ab & QS_LONG_MASK];
		      }
		}
	  }

	if(!odd_nums)
	  { 
	    //fill that bit pattern into sieve[pn][][] 
	    //sieve[pn][b][a0] has the same semantics as help[b][a0], but here, a0 runs from 0 to kp-1 and all bits are filled.
	    for(b = 1; b < bb; b++)
	      { 
		bit_array *si = sieve[pn][b];
		bit_array *he = help[b];
		long p1 = (QS_LONG_LENGTH/p + 1) * p;
		long diff_shift = p1 & QS_LONG_MASK;
		long diff = QS_LONG_LENGTH - diff_shift;
		bit_array diff_mask = ~(all_ones<<diff);
		long wp = p1>>QS_LONG_SHIFT;
		for(a = 0; a < wp; a++) 
		  si[a] = he[a];
		long a1;
		for(a1 = a ; a < p; a++)
		  {
		    he[a1] |= (he[(a1 == wp) ? 0 : a1 + 1] & diff_mask)<<diff_shift;
		    si[a] = he[a1];
		    if(a1 == wp) 
		      a1 = 0; 
		    else a1++;
		    he[a1] >>= diff;
		  } 
	      }
	  }

	if(odd_nums || use_opt)
	  { //fill in sieve2 for the case of only odd numerators */
	    for(b = 1; b < bb; b++)
	      { 
		bit_array *si = sieve2[pn][b];
		bit_array *he = help2[b];
		long p1 = (QS_LONG_LENGTH/p + 1) * p;
		long diff_shift = p1 & QS_LONG_MASK;
		long diff = QS_LONG_LENGTH - diff_shift;
		bit_array diff_mask = ~(all_ones<<diff);
		long wp = p1>>QS_LONG_SHIFT;
		//copy the first chunk from help2[b][] into sieve2[pn][b][] 
		for(a = 0; a < wp; a++) si[a] = he[a];
		//now keep repeating the bit pattern, rotating it in help2 
		long a1;
		for(a1 = a ; a < kp; a++)
		  { 
		    he[a1] |= (he[(a1 == wp) ? 0 : a1 + 1] & diff_mask)<<diff_shift;
		    si[a] = he[a1];
		    if(a1 == wp) a1 = 0; else a1++;
		    he[a1] >>= diff;
		  } 
	      }      
	  }

	// implement check for denominator divisble by p 
	{ 
	  bit_array pattern = bit_zero;
	  if((degree & 1L) || squares[pnn[pn]][coeffs_mod_p[pnn[pn]][degree]])
	    pattern = ~bit_zero;
	  if(!odd_nums)
	    for(a = 0; a < kp; a++) sieve[pn][0][a] = pattern;
	  if(odd_nums || use_opt)
	    for(a = 0; a < kp; a++) sieve2[pn][0][a] = pattern;
	}
      }    
  }
}

long qsieve::sift(long b)
     //In M. Stoll's code: "The following procedure is the heart of the matter. The overall speed of the program highly depends on the quality of the code the compiler produces for the innermost loops. It is advisable to force the compiler to put the variables surv, siv0, siv1 into registers. With gcc-2.7.2.1, I got a speedup of more than 30% ! -- MS"
     //For now, the restrictions on the registers used on a ix86 machine have been removed.  They consist in assigning surv to "%ecx", siv1 to "%ebx", i and siv0 to "%edi", j to "%esi". ("on an ix86 machine with sufficiently new gcc; this will produce much better code in the sift function, resulting in a considerable speedup (> 30% on my machine -- MS)")
{
#if DEBUG_QS>0
  cout<<"\nEntering sift b="<<b<<endl;
#endif
  long total=0;
  halt_flag=0;

  int use_odd_nums; //flag that tells whether to use only odd numerators
  use_odd_nums = odd_nums || (!(b&1) && use_opt);
  { 
    long n;
    if(use_odd_nums)
      for(n = 0; n < sieve_primes2; n++)
	{ 
	  long pn = pnn[n], p = prime[pn];
	  sieves[n].ptr = &sieve2[n][b%p][0];
	}
    else
      for(n = 0; n < sieve_primes2; n++)
      { 
	long pn = pnn[n], p = prime[pn];
        sieves[n].ptr = &sieve[n][b%p][0];
      }   
  }

  long low, high, range;
  long w_low, w_high;
  long w_low0, w_high0;
  compute_bc = 1;
  long k;

#if DEBUG_QS>0
  cout<<"num_inter="<<num_inter<<endl;
#endif

  for(k = 0; k < num_inter; k++)
    { //Determine relevant interval [low, high[ of numerators. 
      { 
	double hb = (double)height/(double)b;
	interval inter = domain[k];
#ifdef DEBUG_QS
	if (verbose)
	  cout<<"Search interval number "<<k<<": ["<<inter.low<<","<<inter.up<<"]"<<endl;
#endif
	if(inter.low <= -hb)
	  low = -height;
	else
	  { 
	    if(inter.low > hb) 
	      return(total); 
	    low = (long)floor(b*inter.low);
	  }
	if(inter.up >= hb)
	  high = height;
	else
	  { 
	    if(inter.up < -hb) 
	      high = -height-1;
	    else 
	      high = (long)ceil(b*inter.up);
	  }
	high++;
      }
      if(use_odd_nums)
	{ 
	  low >>= 1;
	  high--; high >>= 1;
	}
#ifdef DEBUG_QS
    cout<<"low = "<<low<<", high = "<<high<<endl;
#endif

      if(low < high)
	{ //Now the range of longwords (= bit_arrays) 
	  w_low = w_floor(low, QS_LONG_LENGTH);
	  w_high = w_ceil(high, QS_LONG_LENGTH);
	  for(w_low0 = w_low; w_low0 < w_high; w_low0 += array_size)
	    { //array_size is the size of the array of survivors
	      w_high0 = w_low0 + array_size;
	      if(w_high0 > w_high) 
		w_high0 = w_high;
	      range = w_high0 - w_low0;
	      {//initialise the bits 
		register bit_array *surv;
		register long i;
		surv = survivors; // &survivors[0];
		if(!use_odd_nums && !b&1)
		  for(i = range; i; i--) *surv++ = QS_HALF_MASK;
		else
		  for(i = range; i; i--) *surv++ = ~bit_zero;
	      }
	      if(w_low0 == w_low)
		survivors[0] &= (all_ones)<<(low - QS_LONG_LENGTH * w_low);
	      if(w_high0 == w_high)
		survivors[range-1] &= (all_ones)>>(QS_LONG_LENGTH * w_high - high);
	      { 
		long num = sift0(b, w_low0, w_high0, use_odd_nums);
		total += num;
#if DEBUG_QS>0
		cout<<"After calling sift0(), total = "<<total
		    <<" and halt_flag = "<<halt_flag<<endl;
#endif
		if(halt_flag) return(total);
	      } 
	    } 
	}
    }
  return(total);
}//end of sift

long qsieve::sift0(long b, long w_low, long w_high, int use_odd_nums)
{//now do the sieving (fast!)
#if DEBUG_QS>0
  cout<<"Entering sift0. w_low="<<w_low<<"; w_high="<<w_high<<endl;
#endif

  long total = 0, n, i, range;
  halt_flag=0;
  range = w_high - w_low;
  for(n = 0; n < sieve_primes1; n++)
    {
      bit_array *sieve_n = sieves[n].ptr;
      long kp = kpa[n];
      long p_low = w_ceil(w_low, kp), p_high = w_floor(w_high, kp);
      register bit_array *surv;
      surv = survivors;

      if(p_high < p_low)
	{ 
	  register bit_array *siv1;
	  register long i;
	  siv1 = &sieve_n[w_low - kp * p_high];
	  for(i = range; i ; i--) 
	    *surv++ &= *siv1++;
	}
      else
	{
	  register bit_array *siv1;
	  register long j;
	  j = kp * p_low - w_low;
	  siv1 = &sieve_n[kp-j];
	  if(j)
	    { 
	      register long i;
	      for(i = j; i; i--) 
		*surv++ &= *siv1++;
	    }
	  for(i = p_high - p_low; i; i--)
	    {
	      register bit_array *siv0;
	      siv0 = siv1;
	      siv1 -= kp;
	      while(siv1 != siv0) *surv++ &= *siv1++;
	    }
	  j = w_high - kp * p_high;
	  if(j)
	    {
	      register long i;
	      siv1 -= kp;
	      for(i = j; i; i--) 
		*surv++ &= *siv1++;
	    } 
	}
    }
  //Check the points that have survived the sieve if they really are points 
  { 
    bit_array *surv0 = survivors; // &survivors[0];
    bit_array nums=0;
    for(i = w_low; i < w_high; i++)
      {
	if((nums = *surv0++))
	  {
	    check_point(nums, b, i, &total, use_odd_nums);
#if DEBUG_QS>0
	    cout<<"After calling check_point(), total = "<<total
		<<" and halt_flag = "<<halt_flag<<endl;
#endif
	    if(halt_flag) return(total);
	  }
      }
  }
  return(total);
}

void qsieve::check_point(bit_array nums, long b, long i, long* total, int use_odd_nums)
     //MS: "Apart from the innermost loops of sift0, the speed of the following procedure influences the overall speed of the program to some extent. If it is faster, then it is possible to use fewer primes for sieving. Since it is possible to determine beforehand a good bound for the size of the multi-precision integers needed for the computation (the result and all intermediate results are at most of absolute value
     //(l^1-norm of the coeffs of f)*(height bound)^(degree[+1])  ), 
     //it might be possible to gain a little bit by using one's own fixed-length multi-precision integers here."
{//this checks a complete bit_array of possible survivors
#ifdef DEBUG_QS
  cout<<"Entering check_point, i="<<i<<";b="<<b<<";nums="<<(long)nums<<endl;
#endif
  sieve_spec *ssp = &sieves[sieve_primes1];
  num_surv1++;
  if(i < 0)
  { 
    long n;
    for(n = sieve_primes2 - sieve_primes1; n && nums; n--)
      { 
	long p = ssp->p;
	long ip = i%p;
	if(ip) ip += p;
	nums &= ssp->ptr[ip];
	ssp++;
      }
  }
  else
  { 
    long n;
    for(n = sieve_primes2 - sieve_primes1; n && nums; n--)
      { 
	long p = ssp->p;
	long ip = i%p;
	nums &= ssp->ptr[ip];
	ssp++;
      }
  }

  if(nums)
    {
      long a, d;
      bit_array test = 1;//used to isolate the corresponding bit in nums
      i <<= QS_LONG_SHIFT;
      if(!use_odd_nums)
	{ d = 1; a = i; }
      else
	{ d = 2; a = 2*i + 1; }
      
      for(test = 1 ; test; a += d, test <<= 1)
	{ //test one bit
	  if((nums & test) && (gcd(a,b)==1))
	    { //a/b is a surviving fraction so far. 
	      //Now check if it really gives a point. 
	      num_surv2++;
	      if(!no_check)
		{ //Compute F(a, b), where F is the homogenized version of f of smallest possible even degree 
		  long k;
		  if(compute_bc)
		    { //compute entries bc[k] = c[k] * b^(degree-k), k < degree
		      bigint tmp=BIGINT(1);
		      for(k = degree-1; k >= 0; k--)
			{
			  tmp *= b;
			  bc[k] = c[k]*tmp;
			}
		      compute_bc = 0;
		    }
		  bigint fff=c[degree];
		  for(k = degree-1; k >= 0; k--)
		    {
		      fff *= a;  
		      fff += bc[k];
		    }
		  if(degree&1 && !use_squares) 
		    fff *= b;
		  //If degree is odd and lcf is 1, b is a square anyway...
		  bigint y;
		  if(isqrt(fff,y))
		    { 
		      if (verbose)
			cout<<"x = "<<a<<"/"<<b<<" gives a rational point."<<endl;
		      halt_flag = (*curve).process(BIGINT(a),y,BIGINT(b));  
#if DEBUG_QS>0
		      cout<<"After calling process(), halt_flag = "<<halt_flag<<endl;
#endif
		      //process point where [a*b^n-1:y:b^n] are the true homogeneous coordinates,
		      //n=degree/2 if degree even (fff=b^n * y^2)
		      //n=degree/2 if degree odd & polyn monic: b is a square
		      //n=(degree+1)/2 otherwise
		      (*total)++;
		      if(halt_flag) return;
		    }
		} //endif(!no_check)
	      else
		{ 
		  if (verbose)
		    cout<<a<<"/"<<b<<" may be a point (no check)."<<endl;
		  halt_flag = (*curve).process(BIGINT(a),BIGINT(0),BIGINT(b));
#if DEBUG_QS>0
		  cout<<"After calling process(), halt_flag = "
		      <<halt_flag<<endl;
#endif
		  (*total)++;
		  if(halt_flag)	return;
		}
	    }//endif ((nums & test) && ...
	}//endfor (test = 1 ; test; a += d, test <<= 1)
    }//endif (nums)
  return;
}

void qsieve::dealloc_sieves()
{
  long i,j;
  long p;
  for (i=0;i<sieve_primes2;i++)
    {
      p=prime[pnn[i]];
      for (j=0;j<p+1;j++)
	{
	  delete[] sieve[i][j];
	  delete[] sieve2[i][j];
	}
      delete[] sieve[i];
      delete[] sieve2[i];
    }
  delete[] sieve;
  delete[] sieve2;
  delete[] survivors;
}

long qsieve::search()
{
  if (verbose)
    {
      cout<<"Entering qsieve::search: y^2 = ";
      long i;
      for (i=degree;i>0;i--)
	cout<<c[i]<<"x^"<<i<<" + ";
      cout<<c[0]<<endl;
    }
  
  if (num_inter==0)
    {
      num_inter=1;
      domain[0].low=-height;
      domain[0].up=height;
    }

  long s = 2*w_ceil(height, QS_LONG_LENGTH);
  array_size <<= 13 - QS_LONG_SHIFT; /* from kbytes to longs */
  if(s < array_size) array_size = s;
  
  long total = 0,b;
  halt_flag=0;
  int proc_infty = 0;
  bigint t;
  
  if(degree&1 || isqrt(c[degree],t))
    { 
      if(points_at_infty) proc_infty = 1;
      check_denom = 0;
    }
  //Can use only squares as denoms if degree is odd and poly is monic 
  if(degree&1 && (c[degree] == 1)) 
    use_squares = 1;
  
  if (proc_infty)
    {
      halt_flag = (*curve).process(BIGINT(0),BIGINT(1),BIGINT(0));
#if DEBUG_QS>0
      cout<<"After calling process(), halt_flag = "<<halt_flag<<endl;
#endif
      total++;
      if (halt_flag)	return total;
    }
  
  init_f();//sieves get allocated
  survivors=new bit_array[array_size];
  
  if(sieve_primes2 > 0 && prec[0].r == 0.0)
    { 
      cout<<"sieve_primes2 = "<<sieve_primes2<<endl;
      cout<<"prec[0].p = "<<prec[0].p<<endl;
      cout<<"prec[0].n = "<<prec[0].n<<endl;
      cout<<"prec[0].r = "<<prec[0].r<<endl;
      if(verbose) 
	cout<<"No solution (prob=0)"<<endl;
      dealloc_sieves();
      return 0; 
    }
  
  //Otherwise, try to find the points
  if (verbose)
    {
      cout<<"Try to find the points up to height "<<height<<endl;
    }
  
  if (use_squares)
    {
      long limit = (long)floor(sqrt((double)b_high));
#ifdef DEBUG_QS
      cout<<"search with use_squares with b_low "<<b_low<<"to limit "<<limit<<endl;
#endif
      for(b = b_low; b <= limit; b++)
	{ 
	  long num = sift(b*b);
	  total+=num;
#if DEBUG_QS>0
	  cout<<"After calling sift(), total = "<<total
	      <<" and halt_flag = "<<halt_flag<<endl;
#endif
	  if (halt_flag)
	    {
	      dealloc_sieves();
	      return total;
	    }
	}
    }
  else
    { 
#ifdef DEBUG_QS
      cout<<"search without use_squares with b_low "<<b_low<<"to b_high "<<b_high<<endl;
#endif
      if(check_denom)
	{ 
#ifdef DEBUG_QS
	  cout<<"check_denom=1="<<check_denom<<endl;
#endif
	  long *forb;
	  for(b = b_low; b <= b_high; b++)
	    { 
	      //check if denominator is excluded 
	      for(forb = &forbidden[0] ; *forb && (b % (*forb)); forb++)
		{};
	      if(*forb == 0)
		{ 
#ifdef DEBUG_QS
	  cout<<"going to sift with b="<<b<<endl;
#endif
		  long num = sift(b);
		  total += num;
#if DEBUG_QS>0
		  cout<<"After calling sift(), total = "<<total
		      <<" and halt_flag = "<<halt_flag<<endl;
#endif
		  if (halt_flag)
		    {
		      dealloc_sieves();
		      return total;
		    }
		}
	    } 
	}
      else
	for(b = b_low; b <= b_high; b++)
	  {
	    long num = sift(b);
#if DEBUG_QS>0
	    cout<<"After calling sift(), total = "<<total
		<<" and halt_flag = "<<halt_flag<<endl;
#endif
	    total+=num;
	    if (halt_flag)
	      {
		dealloc_sieves();
		return total;
	      }
	  }
    }
  dealloc_sieves();
  return(total); //number of points found 
}

long qsieve::search(double h_lim)
{
  set_height(h_lim);
  return search();
}

/*
long qsieve::search_range(double xmin, double xmax, double h_lim)
{
  set_lower(xmin);
  set_upper(xmax);
  set_height(h_lim);
  return search();
}
*/

