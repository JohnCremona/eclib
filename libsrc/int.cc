// FILE INT.CC -- implementation of some of the interface functions in int.h

#include <assert.h>
#include "eclib/int.h"

INT::INT(const std::string& s)
{
  fmpz_init(z);
  fmpz_set_str(z, s.c_str(), 10);
}

std::ostream& operator<<(std::ostream& s, const INT& a)
{
  char* st = fmpz_get_str(NULL, 10, a.z);
  s << std::string(st);
  flint_free(st);
  return s;
}

std::istream& operator>>(std::istream& s, INT& x)
{
  // std::cout << "Reading into an INT..." << std::endl;
  x = 0;
  int digit = 0, neg = 0;
  char c;
  s >> std::ws;
  c = s.peek();
  if (c=='-')
    {
      neg = 1;
    }
  else if (c=='+')
    {
      ;
    }
  else
    {
      digit = c - '0';
      // std::cout << "First character read was " << c << ", digit = " <<digit << std::endl;
      if ((digit<0) || (digit>9))
        {
          std::cerr << "Invalid input character to INT: " << c << std::endl;
          return s;
        }
      x = digit;
      // std::cout << "Setting x to " << x << std::endl;
    }
  // std::cout << "First character read was " << c << ", neg = " << neg << ", digit = "
  //           << digit << ", so far x = " << x << std::endl;
  s.get();
  c = s.peek();
  digit = c - '0';
  while ((digit>=0) && (digit<=9))
    {
      x *= 10;
      x += digit;
      // std::cout << "Next character read was " << c << ", digit = " << digit
      //           << ", so far x = " << x << std::endl;
      s.get();
      c = s.peek();
      digit = c - '0';
    }
  if (neg)
    x = -x;
  // std::cout << "Final character read was " << c << ", final x = " << x << std::endl;
  return s;
}

int divrem(const INT& a, const INT& b, INT& quo, INT& rem)
{
  fmpz_fdiv_qr(quo.z, rem.z, a.z, b.z);
  return is_zero(rem);
}

INT rounded_division(const INT& a, const INT& b, int round_down)
{
  INT q, r;
  divrem(a,b,q,r);
  INT r2 = 2*r;
  if (round_down)
    {
      // We want -b <= r2 < +b, so q=round(a.b) with halves going *down*
      if (r2<-b)
        q-=1;
      else
        if (r2>=b)
          q+=1;
    }
  else
    {
      // We want -b < r2 <= +b, so q=round(a.b) with halves going *up*
      if (r2<=-b)
        q-=1;
      else
        if (r2>b)
          q+=1;
    }
  return q;
}

std::vector<INT> pdivs(const INT& a)
{
  fmpz_factor_t f;
  fmpz_factor_init(f);
  fmpz_factor(f, a.z);
  std::vector<INT> ans;
  for (int i =0; i< f->num; i++)
    {
      INT p(f->p + i);
      ans.push_back(p);
    }
  fmpz_factor_clear(f);
  return ans;
}

std::vector<INT> sqdivs(const INT& a)
{
  fmpz_factor_t f;
  fmpz_factor_init(f);
  fmpz_factor(f, a.z);
  std::vector<INT> plist;
  std::vector<int> elist;
  int nd = 1;
  for (slong i =0; i< f->num; i++)
    {
      INT p;
      fmpz_set(p.z, f->p + i);
      plist.push_back(p);
      int e = (f->exp[i])/2;
      elist.push_back(e);
      nd *= (1+e);
    }
  fmpz_factor_clear(f);
  std::vector<INT> dlist(1, INT(1));
  dlist.resize(nd);
  nd = 1;
  auto pr = plist.begin();
  auto ei = elist.begin();
  while(pr!=plist.end())
    {
      INT p=*pr++;
      int e=*ei++;
      for (int j=0; j<e; j++)
        for (int k=0; k<nd; k++)
          dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
      nd*=(e+1);
    }
  return dlist;
}

INT invmod(const INT&a, const INT& p)
{
  INT b;
  int r = fmpz_invmod(b.z, a.z, p.z);
  assert (r);
  return b;
}

long invmod(const INT&a, long p)
{
  fmpz_t b, P;
  fmpz_init_set_si(P,p);
  fmpz_init(b);
  int r = fmpz_invmod(b, a.z, P);
  assert (r);
  long ans(fmpz_get_si(b));
  fmpz_clear(b);
  fmpz_clear(P);
  return ans;
}

long I2long(const INT& a)
{
  if (a.is_long())
    return (long)fmpz_get_si(a.z);
  assert(0 && "INT does not fit into a long");
}

INT sqrt_mod_p(const INT& a, const INT& p)
{
  INT b;
  sqrt_mod_p(b, a, p);
  return b;
}
