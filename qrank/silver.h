double silverman_bound(const Curvedata& CD);
double siksek_bound(const Curvedata& CD);
inline double height_constant(const Curvedata& CD)
{
  //  return silverman_bound(CD);
  double b1 = silverman_bound(CD), b2 = siksek_bound(CD);
  return min(b1,b2);
}
