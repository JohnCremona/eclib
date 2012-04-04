// sieve_search.h: declarations of classes point_processor and qsieve
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
 
// File        : sieve_search.h
// Author      : Sophie Labour
//		 Adaption of M. Stoll's code
// Last change : SL, Apr 21 1999, first definitions
// Last change : JC, Dec 2 1999, adapted for direct use in mwrank et al.
// Last change : JC, Aug 19 2002, adapted for stdc++ library

#ifndef SIEVE_SEARCH_H
#define SIEVE_SEARCH_H

class point_processor { //An abstract class to be used as interface in qsieve
 public:
  virtual int process(const bigint& x,const bigint& y,const bigint& z)=0;
  // where x/z is the point found and y the square root of f(a,b);
  // #!# different if the degree is even or odd, poly monic or not;
  // true homogeneous coordinates are [a*b^n-1:y:b^n] where n=degree/2 if 
  // degree even or degree odd & polyn monic (b is a square so b^n is 
  // (sqrt(b))^degree), n=(degree+1)/2 otherwise 
  // if no_check option, process(x,0,z) is called
  // the int return value is to signal the search to stop (if set to 1)
  virtual ~point_processor() {;}
};

class point_printer : public point_processor {
 public:
  point_printer() {};
  ~point_printer() {};
  int process(const bigint& xx, const bigint& yy, const bigint& zz) 
    {cout<<"x= "<<xx<<" y= "<<yy<<" z= "<<zz<<" is a point."<<endl; return 0;}
};

class point_counter : public point_processor {
  int tally;
 public:
  point_counter() {tally=0;};
  ~point_counter() {};
  int process(const bigint& xx, const bigint& yy, const bigint& zz) 
    {tally++; return 0;}
  int get_tally() {return tally;}
};

// *****************************************************
//          class qsieve
// *****************************************************

#define QS_VERSION "This is an adaptation of ratpoints-1.4 by M. Stoll"

#define QS_MAX_DEGREE 10  // max. degree of f(x) 
#define QS_DEFAULT_SIZE 10     // Default value for the array_size option 
#define QS_DEFAULT_TECH 1      // Default value for the tech option 

#define QS_LONG_LENGTH (8 * sizeof(unsigned long))
#define QS_LONG_SHIFT ((QS_LONG_LENGTH == 16) ? 4 : \
                    (QS_LONG_LENGTH == 32) ? 5 : \
		    (QS_LONG_LENGTH == 64) ? 6 : 0)
#define QS_LONG_MASK (~(-1L<<QS_LONG_SHIFT))
#define QS_HALF_MASK ((bit_array)(~((unsigned long)(-1L) / 3)))

     //#define PERCENT 0.6
     //#define SIEVE_PRIMES 19
#define QS_NUM_PRIMES 53 //42 //53
#define QS_MAX_PRIME 251 //191 //251
     //#define QS_MAX_PRIME_EVEN 256
#define QS_DEFAULT_HEIGHT 10

typedef unsigned long bit_array;
#define bit_zero ((bit_array)0)
#define all_ones ((bit_array)(-1L))
typedef struct {long p; double r; long n;} entry;
typedef struct {double low; double up;} interval;
typedef struct {long p; bit_array *ptr;} sieve_spec;


class qsieve {
  bigint c[QS_MAX_DEGREE+1];
  point_processor* curve;
  int degree;
  int verbose;
  
  //"global variables"
  char init_made;
  static long prime[QS_NUM_PRIMES];
  bit_array* bits;
  char** squares;
 //squares[pn][x] = 1 if x is square mod prime[pn], 0 if not 

  char** is_f_square;//[QS_NUM_PRIMES][QS_MAX_PRIME];
  // is_f_square[pn][x] = 1 if f(x) is square mod prime[pn], 0 if not 
  long* pnn;//[QS_NUM_PRIMES];
  // holds the number of the primes used for the sieve
  long* kpa;//[QS_NUM_PRIMES];    
  // stores suitable multiples of sieving primes, depending on tech

  // ratios2[d] gives the speed ratio of testing/sieving for degree d.
  // The following values are reasonable for polynomials with coefficients
  // the size as in the examples, on a Pentium. 
  static double ratio1_def;
  static double ratios2[QS_MAX_DEGREE+1];
  static double ratio1;  /* This will be set if the -r option is given */
  static double ratio2;  /* This will be set if the -R option is given */

  bit_array*** sieve;//[QS_NUM_PRIMES][QS_MAX_PRIME_EVEN][QS_MAX_PRIME_EVEN];
  //bit k of sieve[pn][b1][a1] = 0 iff a/b is excluded mod p = prime[pnn[pn]], when  a = a1*QS_LONG_LENGTH + k mod p and b = b1 mod p. 
  bit_array*** sieve2;//[QS_NUM_PRIMES][QS_MAX_PRIME_EVEN][QS_MAX_PRIME_EVEN];
  // bit k of sieve2[pn][b1][a1] = 0 iff (2*a + 1)/b is excluded mod p = prime[pnn[pn]], when  a = a1*QS_LONG_LENGTH + k mod p and b = b1 mod p and b is even.

  sieve_spec sieves[QS_NUM_PRIMES];

  long coeffs_mod_p[QS_NUM_PRIMES][QS_MAX_DEGREE+1];
  bigint bc[QS_MAX_DEGREE+1]; //The coefficients of f, multiplied by powers of b
  bigint fff, tmp, tmp2;   //Some multi-precision integer variables 
  long sieve_primes1;  //The number of primes used for the first sieving stage
  long sieve_primes2; //The number of primes used for the both sieving stages
  bit_array *survivors;
  entry prec[QS_NUM_PRIMES];  //This array is used for sorting in order to determine the `best' sieving primes.

  long height;         //The height bound
  long w_height;       //The height bound divided by the word length in bits 
  int points_at_infty; //A flag saying if we should look for points at infty
  long b_low, b_high;  //Lower and upper bounds for the denominator
  int halt_flag;       // gets set when searching can stop
  int use_opt;         //A flag saying if to use optimised sieving for even denominators
  int no_check;        //A flag that says whether to omit the final check
  long array_size;     //The size of the survivors array (in longs: converted from kbytes to longs at beginning of search)
  long num_inter;      //The number of intervals in the search region 
  interval domain[QS_MAX_DEGREE+1];  //This contains the intervals representing the search region
  static long tech;
  int check_denom; //A flag saying if we should check if the denom is divisible by a `forbidden divisor'
  long forbidden[QS_NUM_PRIMES+2]; //The forbidden divisors, a zero-terminated array.
  int use_squares; //A flag that is set when we can use only squares as denominators (odd degree & polynomial monic)
  int odd_nums;    //A flag that says that we need only consider odd numerators
  int compute_bc;

  long num_surv1;  //Used to count the survivors of the first stage
  long num_surv2;  //Used to count the survivors of the second stage
  

  long w_floor(long a,long b)
    { if (a < 0) return -(1 + (-a-1)/b); else return a/b;}; 
  long w_ceil(long a,long b)
    { if (a <= 0) return -((-a)/b); else return 1+(a-1)/b;};

  void init_data();
  void init_all();//init bits and squares
  void init_f();//init is_f_square and sieve

  long sift(long b);
  long sift0(long b, long w_low, long w_high, int use_odd_nums);
  void check_point(bit_array nums, long b, long i, long* total, int use_odd_nums);
  //  void a_search(const long& amin, const long& amax);
  //  void a_simple_search(const long& amin, const long& amax);
  
  void dealloc_sieves();

 public:
  qsieve() {;}
  qsieve(point_processor* acurve, int deg, vector<bigint> coeff, int verb=0);
  qsieve(point_processor* acurve, int deg, vector<bigint> coeff, bigfloat h_limx, 
	 int verb=0); 
//qsieve(point_processor* acurve, int deg, vector<bigint> coeff, double h_lim, 
//       double up, double low, int verb=1);
  ~qsieve();

  void set_intervals(vector<double> interv,int nb_bnd,int start_low,int pos_x_only=0);
  //the search will be made on the intervals thus defined:
  //nb_bnd bounds are given in interv, the first one being a lower bound 
  //for the first interval if start_low=1, and the upper bound for the first
  // interval (the lower bound being -height) otherwise.
  //we must have interv[i]<interv[i+1], otherwise, behaviour undefined
  //there can be at most QS_MAX_DEGREE+1 intervals.
  //pos_x_only is 1 if only positive solutions are wanted
  //IMPORTANT: If height is to be different than QS_DEFAULT_HEIGHT and
  //constructor was not called with height argument, set_height() must 
  //be called before set_intervals()

  //All the following set_ functions (and the others, up to set_b_high) 
  //should be called (if called at all) BEFORE set_intervals, and all of 
  //those, before calling a search.

  //Sets the (naive logarithmic) height to which the search is performed 
  //for the x coordinate
  void set_height(double h_lim) 
    {height=(long)(floor(exp(h_lim)));w_height=(height-1)/QS_LONG_LENGTH + 1;};

  void set_sieve_primes1(long sp1);
  //sets the number of primes used for the first stage of sieving to sp1, 
  //unless sieve_primes2 has already been set to a smaller value.  
  //Maximum is QS_NUM_PRIMES.

  void set_sieve_primes2(long sp2);
  //sets the total number of primes used for sieving to sp2, unless 
  //sieve_primes1 has already been set to a greater value.  
  //Maximum is QS_NUM_PRIMES.

  //`ratio1' is the ratio of running time of the second versus the first 
  //stage of sieving (per bit). This is used to find the optimal number 
  //of sieving primes for the first stage automatically.
  void set_ratio1(double r1)
    { if (r1<0.0) 
      cout<<"qsieve::set_ratio1():ratio must be >0"<<endl; 
    else ratio1=r1; };
  //`ratio2' is the ratio of running time needed for checking if x-coordinates
  //give rise to points versus one step of the second sieving stage. 
  //This is used to find the optimal number of sieving primes for the second 
  //stage automatically
  void set_ratio2(double r2)
    { if (r2<0.0) 
      cout<<"qsieve::set_ratiow():ratio must be >0"<<endl; 
    else ratio2=r2; };
  //`size' is the size (in kilobytes) used for the array of bits in which 
  //the sieving is done. This affects the amount of memory used by the 
  //program. A smaller value might give better performance if the array can 
  //be held in the cache in its entirety. So it is to be expected that the 
  //effect of this parameter depends heavily on the cache size and speed. 
  //Default is DEFAULT_SIZE = 10. 
  void set_array_size(long as)
    { if (as<0) 
      cout<<"qsieve::set_array_size(): size must be >0"<<endl; 
    else array_size=as;};//size of survivors array in kbytes
  void be_quiet() {verbose=0;};
  //to switch off the check if a survivng candidate is really a point. 
  //This is useful when you use another program and this program does the 
  //test for itself
  void set_no_check() {no_check=1;};
  //If we do not want it
  void no_points_at_infinity() {points_at_infty=0;};
  //lower bound for the denominator (default:1)
  void set_b_low(long b) {if (b>0) b_low=b;};
  //upper bound for the denominator (default:height)
  void set_b_high(long b) {if (b>0) b_high=b;};


  //all perform sieve-assisted fast search according to the info provided:
  //h_lim(logarithmic, default:5), b_low,b_high(bounds z, optional)
  //intervals for x/z, optional, must be set AFTER height
  //the sieve uses height, so it is set every time search is called
  long search();
  long search(double h_lim);
  
};

/* This is a comparison function needed for sorting in order to determine
   the `best' primes for sieving. */
int compare_entries(const void *a, const void *b);

#endif
