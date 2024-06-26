// periods.h : class for integrating newforms
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

#ifndef _ECLIB_PERIODS_H
#define _ECLIB_PERIODS_H      1

#include <eclib/newforms.h>

class character {
private:
  long modul;
  vector<int> chartable;
  void init();
public:
  character(long m=1);
  void reset(long m);
  long modulus(void) {return modul;}
  int operator()(long n) {return chartable[n%modul];}
};

class summer {
protected:
  bigfloat sum1, sum2;  // sum2 not necessarily used
  long limit, limit1, limit2;
  bigfloat rootlimit, rootmod, factor, factor1, factor2, rp, ip;
  long type;
  long N;  vector<long> aplist;  vector<long> primelist;
  vector<long> an_cache;  // holds a_n for n up to rootlimit
  vector<long> a2p_cache;  // holds a_n for n=2^e up to rootlimit
  vector<long> a3p_cache;  // holds a_n for n=3^e up to rootlimit
  vector<long> a5p_cache;  // holds a_n for n=5^e up to rootlimit
  vector<long> a7p_cache;  // holds a_n for n=7^e up to rootlimit
  long n2p,n3p,n5p,n7p;
  void initaplist(const level* iN, const vector<long>& apl);
  virtual bigfloat func1(long n) {return to_bigfloat(0);}
  virtual bigfloat func2(long n) {return to_bigfloat(0);}
  virtual void use(long n, long an)=0;
  void use1(long n, long an);
  void use2(long n, long an);
  void use2357(long n, long an);
  void add(long n, long pindex, long y, long z);
  void add2357(long n, long pindex, long y, long z);
  void sumit(void);             // do the sum
public:
  virtual ~summer() {;}
  virtual void compute(void)=0;         // ditto with pre and post-processing
  bigfloat rper(void) {return rp;}
  bigfloat iper(void) {return ip;}
  Cperiods getperiods() {Cperiods per(rp,ip,type);  return per;}
};

class periods_via_lfchi :public summer {
private:
  character chi1, chi2;
  long mplus, mminus, dp0;
  void use(long n, long an) override {use2(n,an);}
  bigfloat func1(long n) override { return to_bigfloat(chi1(n)) * pow(factor1,to_bigfloat(n)); }
  bigfloat func2(long n) override { return to_bigfloat(chi2(n)) * pow(factor2,to_bigfloat(n)); }
public:
  periods_via_lfchi (const level* iN, const newform* f);
  void compute(void) override;
};

class periods_direct :public summer {
private:
  long eps_N, a, b, c, d;
  bigfloat theta1,theta2;
  vector<bigfloat> ctab, stab;  // array of cos(j*2*pi/d), sin(j*2*pi/d), for j mod d
  void use(long n, long an) override;

public:
  periods_direct (const level* iN, const newform* f);
  void compute(void) override;
  void compute(long ta, long tb, long tc, long td);
                               // period of (a,b;Nc,d) in Gamma_0(N)
};

class part_period :public summer {
private:
  bigfloat efactor,x0,y0,xn;
  bigfloat func1(long n)  override
  { xn=to_bigfloat(n); efactor = exp(-xn*y0);
    return efactor*cos(xn*x0); }
  bigfloat func2(long n)  override
  { return efactor*sin(xn*x0); }
  void use(long n, long an) override
  { use2(n,an);}

public:
  part_period (const level* iN, const newform* f);
  ~part_period () {;}
  void compute(const bigcomplex& z0);
  void compute() override;
  bigcomplex getperiod() {return bigcomplex(rp,ip);}
};

class ldash1 : public summer {
private:
  long r;
  int computed;
  bigfloat ld1;
  bigfloat G(bigfloat x);  // G_r(x)
  void init(const level* N, const vector<long>& f_aplist, long f_sfe, const rational& f_loverp);
  void use(long n, long an) override {use1(n,an);}
  bigfloat func1(long n) override;
public:
  ldash1 (const level* iN, const newform* f);
  ldash1 (const newforms* nf, long i);  // the i'th newform
  void compute(void) override;
  long rank() {compute(); return r;}
  bigfloat value() {compute(); return ld1;}
//
// NB this value is equal to r!*L^{(r)}(f,1) -- note the r! factor!
//
};

class lfchi : public summer {
private:
  long limit0;
  bigfloat val;
  character chi;
  bigfloat func1(long n) override
  { return chi(n)*pow(factor1,to_bigfloat(n));}
  void use(long n, long an) override
  {use1(n,an);}
public:
  lfchi (const level* iN, const newform* f);
  void compute(long ell);
  void compute(void) override {} // not called but has to exist;
  bigfloat value(void) {return val;}
  bigfloat scaled_value(void) {return sqrt(to_bigfloat(chi.modulus()))*val;}
};

vector<long> resort_aplist(const level* iN,
			   const vector<long>& primelist,
			   const vector<long>& apl);

bigfloat myg1(bigfloat x);

#endif
