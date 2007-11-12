#include <iostream.h>
#include <LiDIA/quadratic_order.h>
int main ()
{
   LiDIA::quadratic_order O (-108708);
   cout << O << O.class_number () << endl;
}
