/**************************************************************************
 *    ratpoints-1.4.c                                                     *
 *                                                                        *
 *    A fast program for finding rational points on curves of the form    *
 *      y^2 = f(x)                                                        *
 *    based on an idea of Noam Elkies, with many improvements by Colin    *
 *    Stahlke and Michael Stoll.                                          *
 *                                                                        *
 *    Call it as follows:                                                 *
 *      ratpoints 'a_0 a_1 ... a_n' h [-dl low_den] [-du up_den]          *
 *                [[-l low1] -u up1 -l low2 -u up2 ... -l lown [-u upn]]  *
 *                [-n num_primes1] [-N num_primes2] [-x]                  *
 *                [-r ratio1] [-R ratio2] [-1] [-i]                       *
 *                [-o] [-s size] [-t tech] [-f format] [-q]               *
 *                                                                        *
 *    where                                                               *
 *                                                                        *
 *    + f(x) = a_n x^n + ... + a_1 x + a_0 with 1 <= n <= 10.             *
 *    + `h' is the bound on the naive height of the x-coordinate.         *
 *    + `low_den' and `up_den' are lower/upper bounds for the denominator.*
 *    + `[-l low1] -u up1 -l low2 -u up2 ... -l lown [-u upn]' restricts  *
 *      the search for points to the union of the intervals               *
 *      [low1, up1], [low2, up2], ..., [lown, upn]. `low1' defaults to    *
 *      -infinity, `upn' defaults to +infinity. Moreover, you must have   *
 *      low1 < up1 < low2 < up2 < ... < lown < upn.                       *
 *      (The `lowi' and `upi' are floating point numbers).                *
 *    + `num_primes1' gives the number of primes used for the first       *
 *      sieving stage. (Default is to find this automatically using the   *
 *      default of `ratio1' below).                                       *
 *    + `num_primes2' gives the number of primes used for the two sieving *
 *      stages together. (Default is to find this automatically using the *
 *      default of `ratio2' below).                                       *
 *    + `ratio1' is the ratio of running time of the second versus the    *
 *      first stage of sieving (per bit). This is used to find the        *
 *      optimal number of sieving primes for the first stage auto-        *
 *      matically. See below for a escription of how to determine a good  *
 *      value of `ratio1' for your machine. This option is ignored when   *
 *      the `-n' option is given.                                         *
 *    + `ratio2' is the ratio of running time needed for checking if      *
 *      x-coordinates give rise to points versus one step of the second   *
 *      sieving stage. This is used to find the optimal number of sieving *
 *      primes for the second stage automatically. See below for a        *
 *      description of how to determine a good value of `ratio2' for your *
 *      machine. This option is ignored when the `-N' option is given.    *
 *    + Use `-x' to switch off the check if a survivng candidate is       *
 *      really a point. This is useful when you use another program (e.g. *
 *      written in Magma, Mathematica or Maple) to read in the results,   *
 *      and this program does the test for itself.                        *
 *    + Use `-1' if you only want one point. The program will stop as     *
 *      soon as it has found one.                                         *
 *    + Use `-i' if you don't want to see points at infinity.             *
 *      Use `-i -1' if you want only one point that is not at infinity.   *
 *    + `-o' switches on the optimisation at even denominators (see       *
 *      below).                                                           *
 *    + `size' is the size (in kilobytes) used for the array of bits in   *
 *      which the sieving is done. This affects the amount of memory used *
 *      by the program. A smaller value might give better performance     *
 *      if the array can be held in the cache in its entirety.            *
 *      So it is to be expected that the effect of this parameter depends *
 *      heavily on the cache size and speed.                              *
 *      Default is DEFAULT_SIZE = 10.                                     *
 *    + `tech' is some technical parameter. On a Pentium, it is best left *
 *      at its default value (= 1), but on a machine with a larger/faster *
 *      cache, you might gain a little bit by setting it to 256 (which is *
 *      NUM_RPIMES_EVEN below). Please tell me waht you get.              *
 *    + `format' is a format string to be given to printf to print points.*
 *      It should use two long integers (first the numerator, then the    *
 *      denominator). Default is "(%ld : %ld)\n".                         *
 *    + Use `-q' to suppress messages other than printing the points      *
 *      found.                                                            *
 *                                                                        *
 *    The order of optional arguments is arbitrary, except for the `-l'   *
 *    and `-u' options that must alternate.                               *
 *                                                                        *
 **************************************************************************
 *                                                                        *
 *    Noam Elkies,   pre-1995 (?)                                         *
 *    Michael Stoll, 04-Jul-1995                                          *
 *      + Introduced bit arrays instead of char arrays to be used in the  *
 *        sieving procedure.                                              *
 *    Colin Stahlke, 22-Apr-1998 (ratpgmp.c)                              *
 *      + Introduced use of the gmp library. It's now possible to check   *
 *        the putative solutions.                                         *
 *      + First attempt to select good primes for sieving.                *
 *    Michael Stoll, 23-Jul-1998: Version 1.0                             *
 *      + Sorting of sieving primes to obtain the best set of them.       *
 *      + Use speed ratio of the two stages (sieving/checking) to find    *
 *        optimal number of sieving primes.                               *
 *      + Lots of options.                                                *
 *      + Optimised innermost loop in sift() on Intel processors.         *
 *    Michael Stoll, 24-Jul-1998: Version 1.1                             *
 *      + Introduced some optimisation (use only odd numerators if        *
 *        possible) and several minor bugfixes. Thanks to John Cremona.   *
 *    Michael Stoll, 18-Sep-1998: Version 1.2                             *
 *      + Fixed a bug introduced in 1.1 (using only odd numerators wasn't *
 *        treated  correctly in sift0(). Thanks to Colin Stahlke.         *
 *    Michael Stoll, 09-Nov-1998: Version 1.3                             *
 *      + Introduced a second sieving stage.                              *
 *    Michael Stoll, 10-Nov-1998: Version 1.4                             *
 *      + Made second stage more efficient.                               *
 *                                                                        *
 *    email: stoll@math.uni-duesseldorf.de                                *
 **************************************************************************/

/*------------------------------------------------------------------------+ 
 | To compile the program, use                                            |
 | `gcc ratpoints-1.4.c -Wall -o ratpoints -O2 -fomit-frame-pointer       |
 |   -lgmp -lgcc -lc -lm' ;                                               | 
 | add  `-m386 -Dx86' if on an ix86 machine with sufficiently new gcc;    |
 | this will produce much better code in the sift0 function, resulting in |
 | a considerable speedup (> 30% on my machine -- MS).                    |
 |                                                                        |
 | You need the GNU gmp library to compile this program.                  |
 +------------------------------------------------------------------------*/

/*------------------------------------------------------------------------+
 | The program works as follows. Basically it implements a loop over all  |
 | denominators b = 1,...,h (where h is the height bound given). Within   |
 | each execution of the loop body, there is a loop over all numerators   |
 | a = -h,...,h (or, more precisely, over a in  Z cap b*I cap [-h, h],    |
 | where I is the union of intervals representing the search region).     |
 | For each numerator a it is tested if  f(a/b)  is a square.             |
 |                                                                        |
 | This test is split into two stages. First the numerators are checked   |
 | whether the result is a square mod p with a list of suitable primes p. |
 | Only the `surviving' numerators are then used to compute f(a/b) and    |
 | to check if we get a square.                                           |
 |                                                                        |
 | There are a number of improvements to this basic scheme.               |
 |                                                                        |
 | (1a) We use bits to represent the individual numerators. In this way,  |
 |     we can sieve as many numerators as bits fit into a long word       |
 |     (usually 32) at the same time, using bit-wise and operations.      |
 |                                                                        |
 | (1b) When this kind of sieving has reduced the candidates considerably,|
 |     we continue the sieving with more primes for each candidate        |
 |     separately.                                                        |
 |                                                                        |
 | (2) We take more primes than necessary for the sieving procedure       |
 |     and determine in a first step which are the `best' ones. This is   |
 |     measured by the ratio of numbers surviving the corresponding step  |
 |     of the sieve (essentially the number of points mod p, divided by   |
 |     2*(p+1)).                                                          |
 |                                                                        |
 | (3) If the polynomial is monic and of odd degree, we need only take    |
 |     perfect squares as denominators.                                   |
 |                                                                        |
 | (4) We can reduce the number of numerators that have to be considered  |
 |     with an even denominator by half by only looking at the odd ones.  |
 |     This is done by giving the `-o' option to the program.             |
 |                                                                        |
 | (5) We can exclude denominators divisible by certain prime powers      |
 |     if we can prove beforehand that there can't be a point with such   |
 |     a denominator. (This essentially means that there is no point at   |
 |     infinity modulo this prime power.) This is only possible when the  |
 |     curve has no rational points at infinity (i.e. even degree and     |
 |     non-square leading coefficient). It might be a good idea to first  |
 |     transform your curve to have this form, if possible with a leading |
 |     coefficient that is = 3 mod 4 and = 2 mod 3.                       |
 |                                                                        |
 | (6) In the same way, we can exclude even numerators if we can prove    |
 |     that f(x) can't be a square for such x (this amounts to showing    |
 |     that f(0), f(2), f(4), f(6) are non-squares mod 8). This has the   |
 |     same effect as the `-o' option has, but for all denominators,      |
 |     cutting the work that has to be done in half.                      |
 |     Thanks to John Cremona for pointing this out.                      |
 |                                                                        |
 | If you want to get even better performance, read the comments further  |
 | down before the `check_point' and `sift' functions. But before trying  |
 | to follow the suggestions there, you might want to play around with    |
 | the `-n' option to find the best number of sieving primes for your     |
 | curve. By default, the program tries to determine the optimal number   |
 | of sieving primes using the ratio1 = ratio1_def of running time of     |
 | the second stage versus the first stage of sieving.                    |
 | The default value works well on a Pentium-S processor. If above        |
 | speed ratio differs sensibly on your machine from that used here, you  |
 | might be able to gain something by adjusting ratio1_def.               |
 |                                                                        |
 | The second ratio, ratio2, should be less critical for the performance. |
 | Ideally, the last step that checks if a candidate really gives rise to |
 | a point is only done for the `real points', and the last few steps of  |
 | the second sieving stage should also deal with only few candidates.    |
 | This implies that the running time of these last steps has no great    |
 | ifluence on the overall performance. However, as a rule, larger        |
 | coefficients will need a larger ratio2, since then the testing         |
 | requires computing with larger numbers. Moreover, polynomials of       |
 | higher degrees need a higher ratio2. This is reflected in the default  |
 | value, which taken as ratios2[degree].                                 |
 |                                                                        |
 | To determine a good value for ratio1, you can proceed as follows.      |
 |  (i) Time the call `ratpoints <coeffs> 100000000 -n 40 -du 1'.         |
 |      Call t1 the time it takes.                                        |
 | (ii) Time the call `ratpoints <coeffs> 10000000 -n 0 -N 1-du 1'.       |
 |      Call t2 the time this one takes.                                  |
 | In both cases, <coeffs> stands for the coefficients of your curve.     |
 | Then a good value for ratio1 is given by                               |
 |   400 * t2/t1 .                                                        |
 | The value need not be very precise, since it is essentially its        |
 | logarithm that enters into the calculation. Check if your value is     |
 | good by giving it to the program via the `-r' option. Look at the      |
 | number of sieving primes you get, and see if the program runs faster   |
 | when you specify one sieving prime more or less (use the `-n' option). |
 | Please tell me when you obtain a value that differs from the default   |
 | one.                                                                   |
 |                                                                        |
 | CHANGE:                                                                |
 | Here is some more detailed explanation of how the program determines   |
 | the number of sieving primes, given such a ratio r.                    |
 | First, the primes are ordered: p1, p2, ..., according to the fractions |
 |  f1 <= f2 <= ...  of numbers that are not excluded modulo the          |
 | corresponding prime. We assume that the sieving at different primes is |
 | independent, so that we expect a proportion of f1*f2*...*fn of the     |
 | numerators to survive the sieve when using the first n primes. If we   |
 | take the time one sieving step (i.e. with one prime) needs as our unit |
 | of measurement, then the expected running time with n sieving primes   |
 | is  n + f1*f2*...*fn*r . As a function of n, this will first decrease, |
 | until we reach a minimum, and then increase again. The minimum is      |
 | attained at the smallest n such that  f1*f2*...*fn*(1-f{n+1})*r < 1    |
 | (since when increasing n by 1, we lose more (1) than we gain (above    |
 | expression)), and this is how the program finds it.                    |
 |                                                                        |
 | Michael Stoll, November 1998.                                          |
 +------------------------------------------------------------------------*/

/*
Examples (timings (user time) on a Pentium-S 200 with -Dx86 under Linux)
 ratpoints '793881 -1073610 181107 125476 -19824 -2704 574' 10000
  time 1.6 s, 97 points, 12+35 primes / 0.8 s auf tim
 ratpoints '313383220164 241518308040 57659086152 -4944866082 -724533076 \
           34187102 4037229' 10000
  time 2.7 s, 71 points, 12+35 primes / 1.2 s auf tim
 ratpoints '900150784 0 -956392943 0 712615534 0 10418625' 10000
  time 3.1 s, 208 points, 13+33 primes / 1.5 s auf tim
 same with -l 0
  time 1.9 s, 104 points / 1.0 s auf tim
 ratpoints '6748448400 5062056960 582010296 -150517656 -15420327 1242798 \
           115929' 10000
  time 3.2 s, 148 points, 13+35 primes / 1.5 s auf tim
 same with -o
  time 3.0 s / 1.5 s auf tim
 ratpoints '6748448400 5062056960 582010296 -150517656 -15420327 1242798 \
           115929' 20000
  time 11.5 s, 161 points, 13+35 primes / 4.7 s auf tim
  (without `-Dx86': 19.6 s, 13 primes)
 same with -o
  time 9.2 s / 3.9 s auf tim
 same with -u -3.599 -l -1.9225
  time 10.9 s, 161 points / 4.4 s auf tim
 same with -u -3.599 -l -1.9225 -o
  time 8.8 s / 3.7 s auf tim
 ratpoints '-260128816 -313104 0 1' 100000 -i
  time 1.2 s, 0 points, 10+42 primes / 0.5 s auf tim
 ratpoints '-260128816 -313104 0 1' 1000000 -i
  time 35.8 s, 1 point / 9.8 s auf tim
 ratpoints '-260128816 -313104 0 1' 1000000 -i -o
  time 27.2 s, 1 point / 7.6 s auf tim
 ratpoints '-260128816 -313104 0 1' 1000000 -i -1
  time 0.3 s, 1 point
 ratpoints '-260128816 -313104 0 1' 1000000 -l 799 -i
  time 0.5 s, 1 point
 ratpoints '-260128816 -313104 0 1' 1000000 -l 799 -i -1
  time 0.2 s, 1 point
 ratpoints '-260128816 -313104 0 1' 5000000 -i -o
  time 307.8 s, 1 point / 81.1 s auf tim
 ratpoints '-260128816 -313104 0 1' 5000000 -l 799 -i -o      
  time 3.6 s, 1 point
 ratpoints '-260128816 -313104 0 1' 10000000 -l 799 -i -o      
  time 10.0 s, 1 point
 ratpoints '-260128816 -313104 0 1' 10000000 -l 799 -i -o -1
  time 1.2 s, 1 point
 (this is curve 873C1 from Cremona's book (2nd edition))
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gmp.h>

/**************************************************************************
 * define                                                                 *
 **************************************************************************/

#define MAX_DEGREE 10       /* max. degree of f(x) */
#define DEFAULT_SIZE 10     /* Default value for the -s option */
#define DEFAULT_TECH 1      /* Default value for the -t option */

#define RATPGMP_VERSION \
  "This is ratpoints-1.4 by Michael Stoll (1998-11-10).\n\n" \
  "Please acknowledge use of the program in published work.\n"

#ifndef DEBUG
#define DEBUG 0
#endif

#if (DEBUG > 0)
# define NUM_PRIMES 17
# define MAX_PRIME 61
# define MAX_PRIME_EVEN 64
#else
# define NUM_PRIMES 53
# define MAX_PRIME 251
# define MAX_PRIME_EVEN 256
#endif /* DEBUG > 0 */

#define FLOOR(a,b) (((a) < 0) ? -(1 + (-(a)-1) / (b)) : (a) / (b))
#define CEIL(a,b) (((a) <= 0) ? -(-(a) / (b)) : 1 + ((a)-1) / (b))

#define LONG_LENGTH (8 * sizeof(unsigned long))
   /* number of bits in an unsigned long */
#define LONG_SHIFT ((LONG_LENGTH == 16) ? 4 : \
                    (LONG_LENGTH == 32) ? 5 : \
		    (LONG_LENGTH == 64) ? 6 : 0)
#define LONG_MASK (~(-1L<<LONG_SHIFT))

#define HALF_MASK ((bit_array)(~((unsigned long)(-1L) / 3)))

/**************************************************************************
 * global variables                                                       *
 **************************************************************************/

#if (DEBUG > 0)
long prime[NUM_PRIMES+1] =  
{3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61};
#else
long prime[NUM_PRIMES+1] =  
{3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,
 103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,
 199,211,223,227,229,233,239,241,251};
#endif
long pnn[NUM_PRIMES+1];  /* This array holds the numbers (= index in prime[])
                            of the sieving primes */
long kpa[NUM_PRIMES];    /* stores suitable multiples of sieving primes,
                            depending on tech */

#ifdef x86
double ratio1_def = 1800;
#else
double ratio1_def = 1000;
#endif

/* ratios2[d] gives the speed ratio of testing/sieving for degree d.
   The following values are reasonable for polynomials with coefficients
   the size as in the examples, on a Pentium. */
double ratios2[MAX_DEGREE+1] =
{ 0, 2.5, 4, 5, 6.5, 7.5, 9, 10.5, 12, 14, 15.5 };

double ratio1 = 0.0;  /* This will be set if the -r option is given */
double ratio2 = 0.0;  /* This will be set if the -R option is given */

typedef unsigned long bit_array;
#define zero ((bit_array)0);
bit_array bits[LONG_LENGTH]; /* An array of bit masks */

typedef struct {long p; double r; long n;} entry;
typedef struct {double low; double up;} interval;

char squares[NUM_PRIMES][MAX_PRIME];
/* squares[pn][x] = 1 if x is a square mod prime[pn], 0 if not */
char is_f_square[NUM_PRIMES][MAX_PRIME_EVEN];
/* is_f_square[pn][x] = 1 if f(x) is a square mod prime[pn], 0 if not */

bit_array sieve[NUM_PRIMES][MAX_PRIME_EVEN][MAX_PRIME_EVEN];
/* bit k of sieve[pn][b1][a1] = 0
   iff a/b is excluded mod p = prime[pnn[pn]],
   when  a = a1*LONG_LENGTH + k mod p and b = b1 mod p. */
bit_array sieve2[NUM_PRIMES][MAX_PRIME_EVEN][MAX_PRIME_EVEN];
/* bit k of sieve2[pn][b1][a1] = 0
   iff (2*a + 1)/b is excluded mod p = prime[pnn[pn]],
   when  a = a1*LONG_LENGTH + k mod p and b = b1 mod p and b is even. */

typedef struct {long p; bit_array *ptr;} sieve_spec;

sieve_spec sieves[NUM_PRIMES];

MP_INT c[MAX_DEGREE+1];  /* The coefficients of f */
MP_INT bc[MAX_DEGREE+1]; /* The coefficients of f, multiplied by powers of b */
MP_INT fff, tmp, tmp2;   /* Some multi-precision integer variables */
long coeffs_mod_p[NUM_PRIMES][MAX_DEGREE+1];
                         /* The coefficients of f reduced modulo the various
                            primes */
entry prec[NUM_PRIMES];  /* This array is used for sorting in order to
                            determine the `best' sieving primes. */

long degree;         /* The degree of f(x) */
long height;         /* The height bound */
long w_height;       /* The height bound divided by the word length in bits */
long sieve_primes1;  /* The number of primes used for the first sieving stage */
long sieve_primes2;  /* The number of primes used for the both sieving stages */
int quiet;           /* A flag saying whether to suppress messages */
int one_point;       /* A flag saying if one point is enough */
int points_at_infty; /* A flag saying if we should look for points at infty */
long b_low, b_high;  /* Lower and upper bounds for the denominator */
int use_opt;         /* A flag saying if to use optimised sieving for
                        even denominators */
int no_check;        /* A flag that says whether to omit the final check */
char *print_format;  /* The printf format for printing points */
long array_size;     /* The size of the survivors array (in longs) */
long num_inter;      /* The number of intervals in the search region */
interval domain[MAX_DEGREE];
                     /* This contains the intervals representing the
                        search region */

long tech = DEFAULT_TECH;

int check_denom = 1; /* A flag saying if we should check if the denom
                        is divisible by a `forbidden divisor' */
long forbidden[NUM_PRIMES+2];
                     /* The forbidden divisors, a zero-terminated array. */
int use_squares = 0; /* A flag that is set when we can use only squares
                        as denominators (odd degree & polynomial monic) */
int odd_nums = 0;    /* A flag that says that we need only consider
                        odd numerators */
int compute_bc;

long num_surv1 = 0;  /* Used to count the survivors of the first stage */
long num_surv2 = 0;  /* Used to count the survivors of the second stage */

bit_array *survivors; /* In this array the sieving takes place */

/**************************************************************************
 * prototypes                                                             *
 **************************************************************************/

void init_main(void);
long find_points(void);
void read_input(long, char *argv[]);
char *scan_mpz(char*, MP_INT*);
static inline int relprime(long, long);
void init_squares(void);
void init_fmodpsquare(void);
void init_sieve(void);
void test_sieve(void);
static inline int check_point(bit_array, long, long, long *, int);
long sift(long);
long sift0(long, long, long, int);
void print_poly(MP_INT*, long);
void message(long, long);
void error(long);

/**************************************************************************
 * main                                                                   *
 **************************************************************************/

int main(int argc, char *argv[])
{
  long total, s;

  if(LONG_SHIFT == 0) error(1);
  init_main();
  /* read input */
  if(argc < 3) error(2);
  read_input(argc-1, &argv[0]);
  if(!quiet)
  { message(0, 0);
    message(5, degree);
    message(6, height);
    message(3, 0);
  }
  s = 2*CEIL(height, LONG_LENGTH);
  array_size <<= 13 - LONG_SHIFT; /* from kbytes to longs */
  if(s < array_size) array_size = s;
  
  /* find and count points */
  total = find_points();
  if(!quiet) { message(12, 0); message(2, total); }
  return(0);
}

/**************************************************************************
 * procedures                                                             *
 **************************************************************************/

/**************************************************************************
 * get at the input                                                       *
 **************************************************************************/

void read_input(long argc, char *argv[])
{
  int l_seen = 0; /* flag for dealing with -l -u */
  { char *s = argv[1];
    degree = 0;
    while((degree <= MAX_DEGREE) && (s = scan_mpz(s, &c[degree]))) degree++;
    degree--;
    if(scan_mpz(s, &fff)) error(3);
  }
  if(degree == 0) error(5);
  if(sscanf(argv[2], " %ld", &height) != 1 || height < 1)
  { error(4); }
  w_height = ((height - 1) / LONG_LENGTH) + 1;
  /* Set global variables to their default values */
  num_inter = 0;       /* No interval up to now */
  sieve_primes1 = -1;  /* automatic determination of this number */
  sieve_primes2 = -1;  /* automatic determination of this number */
  no_check = 0;        /* do the check by default */
  quiet = 0;           /* don't be quiet */
  one_point = 0;       /* look for all points */
  points_at_infty = 1; /* also find points at infinity */
  b_low = 1;           /* denominators go from 1 to h */
  b_high = height;
  use_opt = 0;         /* don't use even denominator optimisation */
  print_format = "(%ld : %ld)\n";
  array_size = DEFAULT_SIZE;
  /* recognise optional args */
  { long i = 3;
    while(i <= argc)
    {
      if(*(argv[i]) != '-') error(6);
      switch(argv[i][1])
      { case 'l': /* lower intevral endpoint */
          if(argc == i) error(7);
          i++;
          if(l_seen) error(7); /* -l -l */
          if(num_inter == MAX_DEGREE) error(7);
          if(sscanf(argv[i], " %lf", &domain[num_inter].low) != 1) error(7);
          if(num_inter > 0 && domain[num_inter-1].up >= domain[num_inter].low)
            error(7);
          i++;
          l_seen = 1;
          break;
        case 'u': /* upper interval endpoint */
          if(argc == i) error(7);
          i++;
          if(!l_seen)
          { if(num_inter == 0) domain[0].low = -height-0.5;
            else error(7); /* -u -u */
          }
          if(sscanf(argv[i], " %lf", &domain[num_inter].up) != 1) error(7);
          if(domain[num_inter].low >= domain[num_inter].up) error(7);
          i++;
          l_seen = 0;
          num_inter++;
          break;
        case 'n': /* number of primes used for first stage of sieving */
          if(argc == i) error(6);
          i++;
          if(sscanf(argv[i], " %ld", &sieve_primes1) != 1) error(6);
          if(sieve_primes1 < 0) sieve_primes1 = 0;
          else
          { if(sieve_primes2 >= 0 && sieve_primes1 > sieve_primes2)
              sieve_primes1 = sieve_primes2;
            else { if(sieve_primes1 > NUM_PRIMES) sieve_primes1 = NUM_PRIMES; }
          }
          i++;
          break;
        case 'N': /* number of primes used for sieving altogether */
          if(argc == i) error(6);
          i++;
          if(sscanf(argv[i], " %ld", &sieve_primes2) != 1) error(6);
          if(sieve_primes2 < 0) sieve_primes2 = 0;
          else
          { if(sieve_primes1 >= 0 && sieve_primes1 > sieve_primes2)
              sieve_primes2 = sieve_primes1;
            else { if(sieve_primes2 > NUM_PRIMES) sieve_primes2 = NUM_PRIMES; }
          }
          i++;
          break;
        case 'r': /* speed ratio 1 */
          if(argc == i) error(6);
          i++;
          if(sscanf(argv[i], " %lf", &ratio1) != 1) error(6);
          if(ratio1 <= 0.0) error(6);
          i++;
          break;
        case 'R': /* speed ratio 2 */
          if(argc == i) error(6);
          i++;
          if(sscanf(argv[i], " %lf", &ratio2) != 1) error(6);
          if(ratio2 <= 0.0) error(6);
          i++;
          break;
        case 's': /* size of survivors array in kbytes */
          if(argc == i) error(6);
          i++;
          if(sscanf(argv[i], " %ld", &array_size) != 1) error(6);
          if(array_size <= 0) error(6);
          i++;
          break;
        case 'f': /* printing format */
          if(argc == i) error(6);
          i++;
          { long l = strlen(argv[i]);
            print_format = malloc((l+1)*sizeof(char));
            strcpy(print_format, argv[i]);
          }
          i++;
          break;
        case 'q': /* quiet */
          quiet = 1;
          i++;
          break;
        case 'x': /* no check */
          no_check = 1;
          i++;
          break;
        case '1': /* only one point */
          one_point = 1;
          i++;
          break;
        case 'i': /* no points at infty */
          points_at_infty = 0;
          i++;
          break;
        case 'o': /* optimisation for even denoms */
          use_opt = 1;
          i++;
          break;
        case 'd': /* Bounds for denom */
          switch(argv[i][2])
          { case 'l': /* lower bound */
              if(argc == i) error(6);
              i++;
              if(sscanf(argv[i], " %ld", &b_low) != 1) error(6);
              if(b_low <= 0) b_low = 1;
              i++;
              break;
            case 'u': /* upper bound */
              if(argc == i) error(6);
              i++;
              if(sscanf(argv[i], " %ld", &b_high) != 1) error(6);
              if(b_high <= 0) b_high = 1;
              i++;
              break;
            default: error(6);
          }
          break;
        case 't':
          if(argc == i) error(6);
          i++;
          if(sscanf(argv[i], " %ld", &tech) != 1) error(6);
          if(tech <= 0) tech = 1;
          else { if(tech > MAX_PRIME_EVEN) tech = MAX_PRIME_EVEN; }
          i++;
          break;
        default: error(6);
  } } }
  if(l_seen)
  /* complete last interval */
  { domain[num_inter].up = height+0.5; num_inter++; }
  if(num_inter == 0)
  /* default interval (effectively ]-infty,infty[) if none is given */
  { domain[0].low = -height-0.5; domain[0].up = height+0.5; num_inter = 1; }
}

/* Read in a long long long integer. Should really be in the library. */
char *scan_mpz(char *s, MP_INT *x)
{
  long neg = 0;
  if(s == NULL || *s == 0) return NULL;
  while(*s == ' ') s++;
  if(*s == 0) return NULL;
  if(*s == '-') {neg = 1; s++;}
  else if(*s == '+') s++;
  mpz_set_si(&tmp2, 0);
  while('0' <= *s && *s <= '9')
  { mpz_mul_ui(&tmp2, &tmp2, 10);
    mpz_add_ui(&tmp2, &tmp2, (long)(*s - '0'));
    s++; }
  if(neg) mpz_neg(&tmp2, &tmp2);
  mpz_set(x, &tmp2);
  return s;
}

/**************************************************************************
 * initialisations                                                        *
 **************************************************************************/

void init_main(void)
{
  bit_array bit;
  long n;
  /* initialise multi-precision integer variables */
  for(n=0; n<=MAX_DEGREE; n++)
  { mpz_init(&c[n]); mpz_init(&bc[n]); }
  mpz_init(&fff);
  mpz_init(&tmp);
  mpz_init(&tmp2);
  /* intialise bits[] */
  bit = (bit_array)1;
  for(n = 0; n < LONG_LENGTH; n++) { bits[n] = bit; bit <<= 1; }
  /* initialise squares[][] */
  init_squares();
}

/**************************************************************************
 * initialise the tables                                                  *
 **************************************************************************/

void init_squares(void)
/* initialise squares[][] */
{
  long a, pn, p;
  for(pn = 0; pn < NUM_PRIMES; pn++)
  {
    p = prime[pn];
    for(a = 0 ; a < p ; a++) squares[pn][a] = 0;
    for(a = 0 ; a < p; a += 2) squares[pn][(a*a)%p] = 1;
         /* if p==2 then this does not work! */
  }
  return;
}

/* This is a comparison function needed for sorting in order to determine
   the `best' primes for sieving. */
int compare_entries(const void *a, const void *b)
{
  double diff = (((entry *)a)->r - ((entry *)b)->r);
  return (diff > 0) ? 1 : (diff < 0) ? -1 : 0;
}

void init_fmodpsquare(void)
/* initialise is_f_square[][] */
{
  long a, pn, p, n, s, np, fdc = 0;

  if(!(degree&1))
  { unsigned long t1, t2;
    /* check if we can exclude even numerators */
    t1 = mpz_fdiv_r_ui(&tmp, &c[0], 8);
    t2 = mpz_fdiv_r_ui(&tmp, &c[1], 4) + 2*mpz_fdiv_r_ui(&tmp, &c[2], 2);
    switch(t1)
    { case 2:
      case 3:
      case 6:
      case 7: if((t2&1) == 0) odd_nums = 1; break;
      case 5: if((t2&3) == 0) odd_nums = 1; break;
  } }
  if(check_denom)
  { unsigned long t, t1, t2, t3;
    /* check if we can exclude denominators divisible by low powers of
       2,3,5, or 7 */
    /* check powers of 2 first */
    t1 = mpz_fdiv_r_ui(&tmp, &c[degree], 8);
    t2 = mpz_fdiv_r_ui(&tmp, &c[degree-1], 4);
    t3 = mpz_fdiv_r_ui(&tmp, &c[degree-2], 2);
    switch(t1)
    { case 2:
      case 3:
      case 6:
      case 7:
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
    /* check 3, 5, 7 */
    t = mpz_fdiv_r_ui(&tmp, &c[degree], 9);
    t1 = t%3;
    t2 = mpz_fdiv_r_ui(&tmp, &c[degree-1], 3);
    if(t1 == 2)
    { forbidden[fdc] = 3; fdc++; }
    else
    { if(t == 3 || t == 6)
      { forbidden[fdc] = (t2 == 0) ? 3 : 9; fdc++; }
    }
        
    t = mpz_fdiv_r_ui(&tmp, &c[degree], 25);
    t1 = t%5;
    t2 = mpz_fdiv_r_ui(&tmp, &c[degree-1], 5);
    if((t1 == 2) || (t1 == 3))
    { forbidden[fdc] = 5; fdc++; }
    else
    { if(t != 0 && t1 == 0)
      { forbidden[fdc] = (t2 == 0) ? 5 : 25; fdc++; }
    }
    
    t = mpz_fdiv_r_ui(&tmp, &c[degree], 49);
    t1 = t%7;
    t2 = mpz_fdiv_r_ui(&tmp, &c[degree-1], 7);
    if((t1 == 3) || (t1 == 5) || (t1 == 6))
    { forbidden[fdc] = 7; fdc++; }
    else
    { if(t != 0 && t1 == 0)
      { forbidden[fdc] = (t2 == 0) ? 7 : 49; fdc++; }
    }
  }
  for(pn = 0; pn < NUM_PRIMES; pn++)
  {
    p = prime[pn];
    /* compute coefficients mod p */
    for(n = 0; n <= degree; n++)
      coeffs_mod_p[pn][n] = mpz_fdiv_r_ui(&tmp, &c[n], p);
    np = squares[pn][coeffs_mod_p[pn][0]];
    is_f_square[pn][0] = np;
    for(a = 1 ; a < p; a++) {
      s = coeffs_mod_p[pn][degree];
      for(n = degree - 1 ; n >= 0 ; n--)
      { s *= a;
        s += coeffs_mod_p[pn][n];
        s %= p;
      }
      if((is_f_square[pn][a] = squares[pn][s])) np++;
    }
    /* Fill array with info for p */
    prec[pn].p = p;
    prec[pn].n = pn;
    if(degree&1 || squares[pn][coeffs_mod_p[pn][degree]])
    { prec[pn].r = ((double)(np*(p-1) + p))/((double)(p*p)); }
    else
    { prec[pn].r = (double)np/(double)p;
         /* denominator divisible by p is excluded */
      if(check_denom && p > 7) { forbidden[fdc] = p; fdc++; } 
    }
  }
  /* sort the array to get at the best primes */
  qsort(prec, NUM_PRIMES, sizeof(entry), compare_entries);
  /* Determine the number of sieving primes:
     First stage: take smallest n such that prod(j<=n) prob(j) * ratio1 < 1.
     Second stage: take smallest N such that (1-prob(N+1)) * ratio2 < 1. */
  if(sieve_primes1 < 0)
  { long n = 0, m = (sieve_primes2 >= 0) ? sieve_primes2-1 : NUM_PRIMES-1;
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
    for(n = sieve_primes1; n < NUM_PRIMES; n++)
      if(ratio2*(1.0 - prec[n].r) < 1.0) break;
    sieve_primes2 = n;
  }
  for(n = 0; n < sieve_primes2; n++)
  { pnn[n] = prec[n].n;
    sieves[n].p = prime[pnn[n]];
  }
  if(!quiet)
  { message(9, 0);
    if(odd_nums) message(11, 0);
    message(4, 0);
    if(sieve_primes2 > 0)
    { message(7, 0); message(8, 0); }
  }
  /* terminate array of forbidden divisors */
  if(check_denom)
  { forbidden[fdc] = 0;
    if(fdc == 0) check_denom = 0;
  }
  if(!quiet && check_denom) message(10, 0);
  return;
}

/* help is an array for temporarily storing the sieving information.
   Bit j of help[b][a0] says whether f(a/b) is a square mod p, where
   a = a0*LONG_LENGTH + j. a runs from 0 to k*p - 1, where k*p is the
   least multiple of p exceeding LONG_LENGTH. */
bit_array help[MAX_PRIME][MAX_PRIME / LONG_LENGTH + 2];
/* help2 plays the same role when using the optimisation for even denoms */
bit_array help2[MAX_PRIME][MAX_PRIME / LONG_LENGTH + 2];

/* initalise sieve[][][] */
void init_sieve(void)
{
  long a, b, i, pn, p, kp, bb, aa;
  
#if (DEBUG >= 2)
  printf("\n sieve:\n");
#endif
  for(pn = 0; pn < sieve_primes2; pn++)
  {
    p = prime[pnn[pn]];
    /* kp = k*p such that k = 1 if p > tech and k*p maximal <= tech otherwise */
    if(p > tech)
      kp = p;
    else
      kp = (tech/p) * p;
    kpa[pn] = kp;
    /* determine which quotients a/b are excluded mod p */
    bb = (b_high >= p) ? p : b_high + 1;
    aa = p>>LONG_SHIFT;
    if(aa == 0) aa = 1;
#if (DEBUG >= 2)
    printf("  p = %ld, kp = %ld, bb = %ld\n", p, kp, bb);
#endif
    if(!odd_nums)
    { for(b = 0; b < bb; b++)
	for(a = 0; a <= aa; a++) help[b][a] = zero;
      if(p < LONG_LENGTH)
	for(a = 0; a < p; a++)
	{ if(is_f_square[pnn[pn]][a])
	    for(b = 1; b < bb; b++)
	    { long ab = (a * b) % p;
	      for(i = 0; i <= LONG_LENGTH / p; i++)
	      {
		help[b][ab>>LONG_SHIFT] |= bits[ab & LONG_MASK];
		ab += p;
	    } }
	}
      else
	for(a = 0; a < p; a++)
	{ if(is_f_square[pnn[pn]][a])
	    for(b = 1; b < bb; b++)
	    { long ab = (a * b) % p;
	      help[b][ab>>LONG_SHIFT] |= bits[ab & LONG_MASK];
	    }
	}
#if (DEBUG >= 3)
      printf(" help(%ld):\n", p);
      for(b = 1; b < bb; b++)
      { printf(" b = %3ld: ", b);
	for(a = 0; a <= aa; a++)
	  printf(" %8.8lx", help[b][a]);
	printf("\n");
      }
#endif
    }
    if(odd_nums || use_opt)
    { /* prepare help2 */
      for(b = 0; b < bb; b++)
        for(a = 0; a <= aa; a++) help2[b][a] = zero;
      if(p < LONG_LENGTH)
        for(a = 0; a < p; a++)
        { if(is_f_square[pnn[pn]][a])
	    for(b = 1; b < bb; b++)
	    { long ab = (a * b) % p;
              ab = (ab&1) ? (ab-1)>>1 : (ab-1+p)>>1;
	      for(i = 0; i <= LONG_LENGTH / p; i++)
	      {
	        help2[b][ab>>LONG_SHIFT] |= bits[ab & LONG_MASK];
	        ab += p;
	    } }
        }
      else
        for(a = 0; a < p; a++)
        { if(is_f_square[pnn[pn]][a])
	    for(b = 1; b < bb; b++)
	    { long ab = (a * b) % p;
              ab = (ab&1) ? (ab-1)>>1 : (ab-1+p)>>1;
	      help2[b][ab>>LONG_SHIFT] |= bits[ab & LONG_MASK];
	    }
        }
    }
    if(!odd_nums)
    { /* fill the bit pattern from help[][] into sieve[pn][][].
	  sieve[pn][b][a0] has the same semantics as help[b][a0],
	  but here, a0 runs from 0 to kp-1 and all bits are filled. */
      for(b = 1; b < bb; b++)
      { bit_array *si = sieve[pn][b];
	bit_array *he = help[b];
	long p1 = (LONG_LENGTH/p + 1) * p;
	long diff_shift = p1 & LONG_MASK;
	long diff = LONG_LENGTH - diff_shift;
	bit_array diff_mask = ~(bit_array)(-1L<<diff);
	long a1;
	long wp = p1>>LONG_SHIFT;
	/* copy the first chunk from help[b][] into sieve[pn][b][] */
	for(a = 0; a < wp; a++) si[a] = he[a];
	/* now keep repeating the bit pattern, rotating it in help */
	for(a1 = a ; a < kp; a++)
	{ 
	  he[a1] |= (he[(a1 == wp) ? 0 : a1 + 1] & diff_mask)<<diff_shift;
	  si[a] = he[a1];
	  if(a1 == wp) a1 = 0; else a1++;
	  he[a1] >>= diff;
      } }
#if (DEBUG >= 3)
      printf(" sieve(%ld):\n", p);
      for(b = 0; b < bb; b++)
      { printf(" b = %3ld: ", b);
	for(a = 0; a < kp; a++)
	  printf(" %8.8lx", sieve[pn][b][a]);
	printf("\n");
      }
#endif
    }
    if(odd_nums || use_opt)
    { /* fill in sieve2 for the case of only odd numerators */
      for(b = 1; b < bb; b++)
      { bit_array *si = sieve2[pn][b];
        bit_array *he = help2[b];
        long p1 = (LONG_LENGTH/p + 1) * p;
        long diff_shift = p1 & LONG_MASK;
        long diff = LONG_LENGTH - diff_shift;
        bit_array diff_mask = ~(bit_array)(-1L<<diff);
        long a1;
        long wp = p1>>LONG_SHIFT;
        /* copy the first chunk from help2[b][] into sieve2[pn][b][] */
        for(a = 0; a < wp; a++) si[a] = he[a];
        /* now keep repeating the bit pattern, rotating it in help2 */
        for(a1 = a ; a < kp; a++)
        { 
          he[a1] |= (he[(a1 == wp) ? 0 : a1 + 1] & diff_mask)<<diff_shift;
	  si[a] = he[a1];
	  if(a1 == wp) a1 = 0; else a1++;
	  he[a1] >>= diff;
      } }      
    }
    /* implement check for denominator divisble by p */
    { bit_array pattern = zero;
      if((degree & 1L) || squares[pnn[pn]][coeffs_mod_p[pnn[pn]][degree]])
        pattern = ~zero;
      if(!odd_nums)
        for(a = 0; a < kp; a++) sieve[pn][0][a] = pattern;
      if(odd_nums || use_opt)
        for(a = 0; a < kp; a++) sieve2[pn][0][a] = pattern;
    }
  }
  return;
}

/**************************************************************************
 * find points by looping over the denominators and sieving numerators    *
 **************************************************************************/

long find_points()
{
  long total = 0, b;
  int print_infty = 0;
  
  if(degree&1 || mpz_perfect_square_p(&c[degree]))
  { if(points_at_infty) print_infty = 1;
    check_denom = 0;
  }
  /* Can use only squares as denoms if degree is odd and poly is monic */
  if(degree&1 && mpz_cmp_si(&c[degree], 1) == 0) use_squares = 1;
  /* initialise is_f_square[][] */
  init_fmodpsquare();
  /* initalise sieve[][][] */
  init_sieve();
  if(print_infty)
  { printf(print_format, 1, 0);
    fflush(stdout);
    total++;
    if(one_point) return(total);
  }
  if(sieve_primes2 > 0 && prec[0].r == 0.0)
  { if(!quiet) message(1,0); return(0); }
  survivors = (bit_array *)malloc(array_size*sizeof(bit_array));
  if(use_squares)
  /* need only take squares as denoms */
  { long limit = floor(sqrt((double)b_high));
    for(b = b_low; b <= limit; b++)
    { long num = sift(b*b);
      total += num;
      if(one_point && num) return(total);
  } }
  else
  { if(check_denom)
    { long *forb;
      for(b = b_low; b <= b_high; b++)
      { /* check if denominator is excluded */
        for(forb = &forbidden[0] ; *forb && (b % (*forb)); forb++) {};
        if(*forb == 0)
        { long num = sift(b);
          total += num;
          if(one_point && num) return(total);
        }
    } }
    else
      for(b = b_low; b <= b_high; b++)
      { long num = sift(b);
        total += num;
        if(one_point && num) return(total);
      }
  }
  return(total);
}

/**************************************************************************
 * check if a and b are relatively prime, used in check_point             *
 **************************************************************************/

static inline int relprime(long m, long n)
{
  if(m < 0) m = -m;
  /* n (the denominator) is always positive here */
  while(1)
  { m %= n;
    if(m == 0) return(n == 1);
    n %= m;
    if(n == 0) return(m == 1);
  }
}

/**************************************************************************
 * check a `survivor' of the sieve if it really gives a point             *
 **************************************************************************/

/*------------------------------------------------------------------------+
 | Apart from the innermost loops of sift0, the speed of the following    |
 | procedure influences the overall speed of the program to some extent.  |
 | If it is faster, then it is possible to use fewer primes for sieving.  |
 | Since it is possible to determine beforehand a good bound for the size |
 | of the multi-precision integers needed for the computation (the result |
 | and all intermediate results are at most of absolute value             |
 |  (l^1-norm of the coeffs of f)*(height bound)^(degree[+1])  ),         |
 | it might be possible to gain a little bit by using one's own           |
 | fixed-length multi-precision integers here.                            |
 +------------------------------------------------------------------------*/

static inline int check_point(bit_array nums, long b, long i, long *total, 
                       int use_odd_nums) 
{ /* to be more precise, this checks a complete bit_array of possible
     survivors */   
  sieve_spec *ssp = &sieves[sieve_primes1];
  num_surv1++;
  if(i < 0)
  { long n;
    for(n = sieve_primes2 - sieve_primes1; n && nums; n--)
    { long p = ssp->p;
      long ip = i%p;
      if(ip) ip += p;
      nums &= ssp->ptr[ip];
      ssp++;
    }
  }
  else
  { long n;
    for(n = sieve_primes2 - sieve_primes1; n && nums; n--)
    { long p = ssp->p;
      long ip = i%p;
      nums &= ssp->ptr[ip];
      ssp++;
    }
  }
  if(nums)
  { long a, k, d;
    bit_array test; /* test has one bit set and is used to isolate
                       the corresponding bit in nums */
    /* a will be the numerator corresponding to the selected bit */
    i <<= LONG_SHIFT;
    if(!use_odd_nums)
    { d = 1; a = i; }
    else
    { d = 2; a = 2*i + 1; }
    for(test = 1 ; test; a += d, test <<= 1)
    { /* test one bit */
      if((nums & test) && relprime(a, b))
      { /* a/b is a surviving fraction so far. Now check if it really
           gives a point. */
	num_surv2++;
        if(!no_check)
        { /* Compute F(a, b), where F is the homogenized version of f
             of smallest possible even degree  */
          if(compute_bc)
          { /* compute entries bc[k] = c[k] * b^(degree-k), k < degree */
            mpz_set_si(&tmp, 1);
            for(k = degree-1; k >= 0; k--)
            { mpz_mul_ui(&tmp, &tmp, b);
              mpz_mul(&bc[k], &c[k], &tmp);
            }
            compute_bc = 0;
          }
          mpz_set(&fff, &c[degree]);
          for(k = degree-1; k >= 0; k--)
          { /* Why isn't there a `mpz_mul_si' ? */
            if(a > 0)
              mpz_mul_ui(&fff, &fff, (unsigned long)a);
            else
            { mpz_mul_ui(&fff, &fff, (unsigned long)(-a));
              mpz_neg(&fff, &fff);
            }
            mpz_add(&fff, &fff, &bc[k]);
          }
          if(degree&1 && !use_squares) mpz_mul_ui(&fff, &fff, b);
            /* If degree is odd and lcf is 1, b is a square anyway... */
#if (DEBUG >= 2)
          printf("(%ld : %ld) -> f = ", a, b);
          mpz_out_str((FILE *)NULL, 10, &fff);
          printf("\n");
#endif
          if(mpz_cmp_si(&fff, 0) == 0 || mpz_perfect_square_p(&fff))
          { printf(print_format, a, b);
            fflush(stdout);
            (*total)++;
            if(one_point) return(1);
          }
        } /* if(!no_check) */
        else
        { printf(print_format, a, b);
          fflush(stdout);
          (*total)++;
          if(one_point) return(1);
        } 
      } /* if((nums & test) && relprime(a, b)) */
    } /* for(test = 1 ; test; a += d, test <<= 1) */
  }
  return(0);
}

/**************************************************************************
 * The sieving procedure itself                                           *
 **************************************************************************/

/*------------------------------------------------------------------------+ 
 | The following procedure is the heart of the matter.                    |
 | The overall speed of the program highly depends on the quality         |
 | of the code the compiler produces for the innermost loops in sift0.    |
 | It is advisable to force the compiler to put the variables             |
 | surv, siv0, siv1 into registers.                                       |
 | With gcc-2.7.2.1, I got a speedup of more than 30% ! -- MS             |
 +------------------------------------------------------------------------*/

long sift(long b)
/* print fractions a/b surviving sieve; return their number. */
{
  long total = 0, low, high, range;
  long w_low, w_high;
  long w_low0, w_high0;
  long k;
  int use_odd_nums; /* flag that tells whether to use only odd numerators */
  
#if (DEBUG >= 1)
  printf("\n sift(b = %ld)\n", b);
#endif
  /* check if using only odd numerators */
  use_odd_nums = odd_nums || (!(b&1) && use_opt);
  { long n;
    if(use_odd_nums)
      for(n = 0; n < sieve_primes2; n++)
      { long pn = pnn[n], p = prime[pn];
        sieves[n].ptr = &sieve2[n][b%p][0];
      }
    else
      for(n = 0; n < sieve_primes2; n++)
      { long pn = pnn[n], p = prime[pn];
        sieves[n].ptr = &sieve[n][b%p][0];
  }   }
#if (DEBUG >= 2)
  { long n;
    for(n = 0; n < sieve_primes2; n++)
      printf(" sieves[%ld] = { %ld, %p } -> %8.8lx\n", 
             n, sieves[n].p, sieves[n].ptr, *(sieves[n].ptr));
    printf("\n");
  }
#endif
  compute_bc = 1;
  for(k = 0; k < num_inter; k++)
  { /* Determine relevant interval [low, high[ of numerators. */
    { double hb = (double)height/(double)b;
      interval inter = domain[k];
      if(inter.low <= -hb)
        low = -height;
      else
      { if(inter.low > hb) return(total); 
        low = ceil(b*inter.low);
      }
      if(inter.up >= hb)
        high = height;
      else
      { if(inter.up < -hb) high = -height-1;
        else high = floor(b*inter.up);
      }
      high++;
    }
    if(use_odd_nums)
    { low >>= 1;
      high--; high >>= 1;
    }
#if (DEBUG >= 1)
    printf(" low = %ld, high = %ld\n", low, high);
#endif
    if(low < high)
    { /* Now the range of longwords (= bit_arrays) */
      w_low = FLOOR(low, LONG_LENGTH);
      w_high = CEIL(high, LONG_LENGTH);
      for(w_low0 = w_low; w_low0 < w_high; w_low0 += array_size)
      { w_high0 = w_low0 + array_size;
        if(w_high0 > w_high) w_high0 = w_high;
        range = w_high0 - w_low0;
        /* initialise the bits */
        {
#ifdef x86
          register bit_array *surv asm("%ecx");
          register long i asm("%edi");
#else
          register bit_array *surv;
          register long i;
#endif
          surv = &survivors[0];
          if(!use_odd_nums && !b&1)
	    for(i = range; i; i--) *surv++ = HALF_MASK;
	  else
	    for(i = range; i; i--) *surv++ = ~zero;
        }
        if(w_low0 == w_low)
          survivors[0] &= (~0UL)<<(low - LONG_LENGTH * w_low);
        if(w_high0 == w_high)
          survivors[range-1] &= (~0UL)>>(LONG_LENGTH * w_high - high);
#if (DEBUG >= 3)
        { long i;
          for(i = 0; i < range; i++) printf(" %8.8lx",survivors[i]);
        }
#endif
#if (DEBUG >= 1)
        printf("\n sift0(%ld, %ld)\n", w_low0, w_high0);
#endif
        { long num = sift0(b, w_low0, w_high0, use_odd_nums);
          total += num;
          if(num && one_point) return(total);
    } } }
  }
  return(total);
}

long sift0(long b, long w_low, long w_high, int use_odd_nums)
{
   /* now do the sieving (fast!) */
  long total = 0, n, i, range;
  range = w_high - w_low;
  for(n = 0; n < sieve_primes1; n++)
  {
    bit_array *sieve_n = sieves[n].ptr;
    long kp = kpa[n];
    long p_low = CEIL(w_low, kp), p_high = FLOOR(w_high, kp);
#ifdef x86
    register bit_array *surv asm("%ecx");
#else
    register bit_array *surv;
#endif
#if (DEBUG >= 4)
    { long i, p = sieves[n].p, r = (range < p) ? range : p;
      printf(" p = %ld:", p);
      printf(" sieve_n = %p (%p); *sieve_n =", sieve_n, &sieve[pnn[n]][b%p][0]);
      for(i = 0; i < r; i++)
        printf(" %8.8lx", sieve_n[i]);
      printf("\n");
      printf(" sieve[n][b%%p][] =");
      for(i = 0; i < r; i++)
        printf(" %8.8lx", sieve[pnn[n]][b%p][i]);
      printf("\n");
      for(i = 0; i < range; i++) printf(" %8.8lx",survivors[i]);
      printf("\n");
    }
#endif
    surv = survivors;
    
    if(p_high < p_low)
    { 
#ifdef x86
      register bit_array *siv1 asm("%ebx");
      register long i asm("%edi");
#else
      register bit_array *siv1;
      register long i;
#endif
      siv1 = &sieve_n[w_low - kp * p_high];
      for(i = range; i ; i--) *surv++ &= *siv1++;
    }
    else
    {
#ifdef x86
      register bit_array *siv1 asm("%ebx");
      register long j asm("%esi");
#else
      register bit_array *siv1;
      register long j;
#endif
      j = kp * p_low - w_low;
      siv1 = &sieve_n[kp-j];
      if(j)
      { 
#ifdef x86
        register long i asm("%edi");
#else
        register long i;
#endif
        for(i = j; i; i--) *surv++ &= *siv1++;
      }
      for(i = p_high - p_low; i; i--)
      {
#ifdef x86
        register bit_array *siv0 asm("%edi");
#else
        register bit_array *siv0;
#endif
        siv0 = siv1;
        siv1 -= kp;
        while(siv1 != siv0) *surv++ &= *siv1++;
      }
      j = w_high - kp * p_high;
      if(j)
      {
#ifdef x86
        register long i asm("%edi");
#else
        register long i;
#endif
        siv1 -= kp;
        for(i = j; i; i--) *surv++ &= *siv1++;
  } } }
#if (DEBUG >= 3)
  for(i = 0; i < range; i++) printf(" %8.8lx",survivors[i]);
  printf("\n");
#endif
  /* Check the points that have survived the sieve if they really are points */
  { bit_array *surv0 = &survivors[0];
    bit_array nums;
    for(i = w_low; i < w_high; i++)
      if((nums = *surv0++))
        if(check_point(nums, b, i, &total, use_odd_nums)
            && one_point)
          return(total);
  }
  return(total);
}

/**************************************************************************
 * output routines                                                        *
 **************************************************************************/

void print_poly(MP_INT *coeffs, long degree)
{
  int flag = 0;
  int i;
  for(i = degree; i >= 0; i--)
  { mpz_set(&tmp, &coeffs[i]);
    if(mpz_cmp_si(&tmp, 0) != 0)
    { if(mpz_cmp_si(&tmp, 0) > 0)
      { printf(flag ? " + " : ""); }
      else
      { printf(flag ? " - " : "- ");
        mpz_neg(&tmp, &tmp);
      }
      flag = 1;
      switch(i)
      { case 0: printf("%s", mpz_get_str((char *) 0, 10, &tmp)); break;
        case 1: if(mpz_cmp_si(&tmp, 1) == 0)
                  printf("x");
                else
                  printf("%s x", mpz_get_str((char *) 0, 10, &tmp));
                break;
        default: if(mpz_cmp_si(&tmp, 1) == 0)
                   printf("x^%d", i);
                 else
                   printf("%s x^%d", mpz_get_str((char *) 0, 10, &tmp), i);
                 break;
  } } }
  printf("\n\n");
  fflush(stdout);
}

void message(long n, long total)
{
  switch(n)
  { case 0: printf("\n%s\n", RATPGMP_VERSION); break;
    case 1: printf("\nprob = 0, hence no solutions.\n"); break;
    case 2: printf("\n%ld rational point pairs found.\n", total); break;
    case 3: printf("Search region:\n  ");
            { long i;
              for(i = 0; i < num_inter; i++)
              { if(i) printf(" U ");
                printf("[%f, %f]", domain[i].low, domain[i].up);
            } }
            printf("\n");
            break;
    case 4: printf("%ld primes used for first stage of sieving,\n", 
                   sieve_primes1);
            printf("%ld primes used for both stages of sieving together.\n",
                   sieve_primes2);
            break;
    case 5: printf("\ny^2 = "); print_poly(c, total); break;
    case 6: printf("max. Height = %ld\n", total); break;
    case 7: { long i;
              printf("Sieving primes:\n First stage: ");
              for(i = 0; i < sieve_primes1; i++)
              { printf("%ld", prime[pnn[i]]);
                if(i < sieve_primes1 - 1) printf(", ");
              }
              printf("\n Second stage: ");
              for( ; i < sieve_primes2; i++)
              { printf("%ld", prime[pnn[i]]);
                if(i < sieve_primes2 - 1) printf(", ");
              }
              printf("\n");
              break;
            }
    case 8:
      printf("Probabilities: Min(%ld) = %f, Cut1(%ld) = %f, ",
              prec[0].p, prec[0].r, 
              prec[sieve_primes1-1].p, prec[sieve_primes1-1].r);
      printf("Cut2(%ld) = %f, Max(%ld) = %f\n\n",
              prec[sieve_primes2-1].p, prec[sieve_primes2-1].r,
              prec[NUM_PRIMES-1].p, prec[NUM_PRIMES-1].r);  break;
    case 9: printf("Using speed ratios %f and %f\n", ratio1, ratio2); break;
    case 10:
    { long n;
      printf("Forbidden divisors of the denominator:\n  ");
      for(n = 0; forbidden[n]; n++)
      { printf("%ld", forbidden[n]);
        if(forbidden[n+1]) printf(", ");
      }
      printf("\n\n");
      break;
    }
    case 11: printf("Even numerators are excluded.\n"); break;
    case 12: printf("\n%ld candidates survived the first stage,\n", num_surv1);
             printf("%ld candidates survived the second stage.\n", num_surv2);
             break;
  }
  fflush(stdout);
}

void error(long errno)
{
  switch(errno)
  { case 1:
      printf("\nUnusual size of `unsigned long' type: %d\n\n", 
             (int)LONG_LENGTH);
      break;
    case 3: printf("\nToo many coefficients.\n\n"); break;
    case 4: printf("\nIncorrect height argument.\n\n"); break;
    case 5: printf("\nThe polynomial must have degree at least 1.\n\n"); break; 
    case 6: printf("\nWrong syntax for optional arguments:\n\n");
    case 2:
      printf("\n");
      printf("Usage: ratpoints 'a_0 a_1 ... a_n' max_height\n");
      printf("                 [-dl low_den] [-du up_den] [-f format]\n");
      printf("                 [[-l low1] -u up1 ... -l lown [-u upn]]\n");
      printf("                 [-n num_primes1] [-N num_primes2] [-x]\n");
      printf("                 [-r ratio1] [-R ratio2]\n");
      printf("                 [-q] [-1] [-i] [-o] [-s size] [-t tech]\n\n");
      break;
    case 7: printf("\nIncorrect interval arguments (not alternating, ");
            printf("too many, or not ordered).\n\n");
            break;
  }
  fflush(stdout);
  exit(errno);
}
