#include <LiDIA/gf_element.h>
#include <LiDIA/eco_prime.h>

using namespace LiDIA; 

main(int, char**)
{
    int p=103;
    //    cout << " p = " ;
    //    cin >> p;
    cout << " p = " << p << endl;
    bigmod::set_modulus(p);
    galois_field F(p);
    gf_element a(F), b(F);
    a.randomize();
    b.randomize();
 
    cout << " a : " << a <<endl;
    cout << " b : " << b <<endl;
 
    bigint order = compute_group_order(a,b);
 
    cout <<" group order : "<< order <<endl;
} 
