//   h1clist.cc  --  Outputs table of curves computing a_i, r, |T|
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
//

#include <iomanip>
#include <eclib/periods.h>
#include <eclib/points.h>

#define AUTOLOOP

#include <eclib/curvesort.h>
#define LMFDB_ORDER       // if defined, sorts newforms into LMFDB order before output
                          // otherwise, sorts newforms into Cremona order before output

#define CURVE_IS_ONE_FIELD // outputs      [a1,a2,a3,a4,a6]
                           // else outputs a1 a2 a3 a4 a6

const scalar modulus(default_modulus<scalar>());

int main(void)
{
  int prec0 = 100;
  int maxprec = 300;
  int prec = prec0;
  int delta_prec = 50;
  set_precision(prec);
 int limit,n=1;
 string code;
#ifdef AUTOLOOP
 cerr<<"Enter first and last N: ";cin>>n>>limit;
 n--; cerr<<endl;
#ifdef CURVE_IS_ONE_FIELD
 cerr << "    N\tC\t#\t[a1,a2,a3,a4,a6]\tr\t|T|\tdeg(phi)" << endl;
#else
 cerr << "    N\tC\t#\ta1 a2  a3 a4          a6\tr\t|T|\tdeg(phi)" << endl;
#endif
 cerr << "----------------------------------------------------------------------" << endl;

 while (n<limit) { n++;
#else
 while (n>1) { cout<<"Enter level: "; cin>>n;
#endif
 if (n>1)
{
  newforms nf(n, modulus, 0);
  int noldap=25;
  nf.createfromdata(1,noldap,0); // do not create from scratch if data absent
  int nnf = nf.n1ds;
  //  int inf = 1;
  if(nnf>0)
    {
      // cout<<"****************************"<<endl;
      // cout<<"N="<<n<<" before sorting: "<<endl;
      // cout<<"****************************"<<endl;
      // for (inf=0; inf<nnf; inf++)
      //   {
      //     code = codeletter(inf);
      //     cout<<code<<": ";
      //     vec_out(cout,nf.nflist[inf].aplist,20);  // outputs at most 20 eigs.
      //     cout<<endl;
      //   }
#ifdef LMFDB_ORDER
  nf.sort_into_LMFDB_label_order();
#else
  nf.sort_into_Cremona_label_order();
#endif
      // cout<<"****************************"<<endl;
      // cout<<"N="<<n<<" after  sorting: "<<endl;
      // cout<<"****************************"<<endl;
      // for (inf=0; inf<nnf; inf++)
      //   {
      //     code = codeletter(inf);
      //     cout<<code<<": ";
      //     vec_out(cout,nf.nflist[inf].aplist,20);  // outputs at most 20 eigs.
      //     cout<<endl;
      //   }
      // cout<<"****************************"<<endl;
    }
 for(int xi=0; xi<nnf; xi++)
   { int i=xi;
     code = codeletter(xi);
     newform& nfi = nf.nflist[i];
     int degphi = nfi.degphi;
     bigfloat rperiod;
     int r = nfi.rank();
     Curve C;
     Curvedata CD;
     CurveRed CR;
     bigint nc; int nt;
     bigint a1,a2,a3,a4,a6;
     prec = prec0-delta_prec;
     while (C.isnull() && (prec<maxprec))
       {
	 prec += delta_prec;
	 set_precision(prec);
	 C = nf.getcurve(i, -1, rperiod);
         CD = Curvedata(C,1);  // The 1 causes minimalization
         CR = CurveRed(CD);
         nc = getconductor(CR);
         if(n!=nc) C = Curve();  // wrong conductor; reset to null curve
       }
     if (C.isnull())
       {
	 cout<<setw(6)<<n<<code<<": bad curve"<<endl;
	 continue;
       }

     CD.getai(a1,a2,a3,a4,a6);
     nt = ntorsion(CD);
     //     cout.form("%4d\t%s\t1\t",n,code);
     cout<<setw(6)<<n<<"\t"<<code<<"\t1\t";
#ifdef CURVE_IS_ONE_FIELD
     cout<<"["<<a1<<","<<a2<<","<<a3<<","<<a4<<","<<a6<<"]";
#else
     cout<<a1<<" "<<a2<<" "<<a3<<" "<<a4<<" "<<a6;
#endif
     cout << "\t" << r << "\t" << nt << "\t" << degphi << endl;
   }
}       // end of if(n)
}       // end of while()
}       // end of main()
