//  Needs updating to use Cperiods class etc
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2007 John Cremona
// 
// This file is part of the mwrank/g0n package.
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
//
// FILE H1DEGPHITEX.CC: (not upgraded) Program to output tex table of deg(phi)

#include <time.h>
#include <fstream.h>
#include "marith.h"
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "oldforms.h"
#include "newforms.h"
#include "h1newforms.h"
#include "periods.h"
#include "getai.h"

#define AUTOLOOP
#define BOOKORDER       // if defined, sorts newforms/curves into order
                        // in the Book (relevant up to 500 only)
#ifdef BOOKORDER
#include "curvesort.cc"
#endif

char* antwerpletters[201][17]=

{{},{},{},{},{},{},{},{},{},{},{},{"BCA"},{},{},{"CDEAFB"},{"CEBFHGDA"},{},
/*17:*/{"CBDA"},{},{"BCA"},{"BADC"},{"BDCAFE"},{},{},
/*24:*/{"BCDAFE"},{},{"BCA","DE"},{"BDAC"},{},{},{"ABCDEFGH"},{},{"BACD"},{"BADC"},
/*34:*/{"ABCD"},{"BCA"},{"ABCD"},{"A","CDB"},{"DEC","AB"},{"BCDA"},
/*40:*/{"BDAC"},{},{"ABCDFE"},{"A"},{"AB"},{"ABDCEFHG"},{"AB"},{},{"BDCAFE"},
/*49:*/{"ABCD"},{"EFGH","ABCD"},{"AB"},{"BA"},{"A"},{"EFD","ACB"},
/*55:*/{"BDCA"},{"CDEF","AB"},{"E","BACD","FG"},{"A","BC"},{},{},{"A"},{"ABCD"},
/*63:*/{"ABCDFE"},{"BCDA"},{"AB"},{"ABCD","EFHG","IJLK"},{"A"},{},{"AB"},
/*70:*/{"ABDC"},{},{"ABDCFE"},{"BA"},{},{"AB","EFGHIJLK","CD"},
/*76:*/{"A"},{"F","DEC","AB"},{"ABCD"},{"A"},{"FEHG","BADC"},
/*81:*/{},{"AB"},{"A"},{"CDEF","AB"},{"AB"},{},{},{"A"},{"C","AB"},
/*90:*/{"MNOP","ABCD","EFGIHJLK"},{"A","BCD"},{"AB","C"},{},
/*94:*/{"AB"},{},{"EFHG","ADBC"},{},{"BADCFE"},{"AB","HIKJ","FG","CDE"},
/*100:*/{"ABCD"},
/*101:*/{"A"},{"EF","GHJILK","ABCD"},{},{"A"},{"ABDC"},{"BC","A","EF","D"},{},
/*108:*/{"AB"},{"A"},{"CD","AB","EF"},{},{"KL","ABDC","EFGHIJ"},{"BA"},
/*114:*/{"ABCD","EF","GHJI"},{"A"},{"E","AB","DC"},{"ABDC"},{"A","BC","D","E"},
/*119:*/{},{"EFHGJI","ABCD"},{"HI","DE","FG","ABC"},{"A"},{"AB","C"},{"BC","A"},
/*125:*/{},{"ABCDEF","GHJILK"},{},{"CD","FE","AB","GH"},{"E","BACD"},
/*130:*/{"EFGH","ABDC","JI"},{"A"},{"AB","CD"},{},{},{"A","B"},{"AB","CD"},
/*137:*/{},{"EF","GHIJ","ABDC"},{"A"},{"AB","C"},{"E","GF","ABCD","I","H"},
/*142:*/{"F","E","AB","CD","G"},{"A"},{"ABCD","EFGHJI"},{"AB"},{},
/*147:*/{"CDEFHG","IJ","AB"},{"A"},{},{"ABCD","GHEF","IJKLMNOP"},{},{"A","B"},
/*153:*/{"C","AB","EFHG","D"},{"CD","EFGH","AB"},{"DE","AB","C"},{"EF","ABCD"},
/*157:*/{},{"E","D","HI","BCA","FG"},{},{"AB","DC"},{"BACD"},
/*162:*/{"KL","GHIJ","ABDC","EF"},{"A"},{},{},{"A"},{},{"BACD","EFGH"},
/*169:*/{},{"AB","HIJK","FG","DE","C"},{"DEFG","ABC","IJ","H"},{"AB"},
/*173:*/{},{"IJ","GH","F","ABCD","E"},{"BA","CDE","FG"},{"C","DEF","AB"},{},
/*178:*/{"AB","CD"},{"A"},{"ABCD"},{},{"EFGH","ABC","J","D","I"},{},
/*184:*/{"C","B","DE","A"},{"D","A","BC"},{"D","BC","A"},{"AB","C"},{},
/*189:*/{"A","CDE","FGH","B"},{"D","C","AB"},{},{"QRTS","ABDC","KLMNPO","EFHGJI"},
/*193:*/{},{"AB"},{"ABDCEFHG","I","K","J"},{"AB","CD"},{"A"},
/*198:*/{"IJLK","EFGH","MNOP","ABCD","QRST"},{},{"B","CD","GHJI","EF","A"}};
 
int main(void)
{
 cout.precision(15);
 int limit,n=1; 
#ifdef AUTOLOOP
 cout<<"Enter first and last N: ";cin>>n>>limit; 
 n--; cout<<endl;
 cout << "\n\\entry{N}{\\#}{[a1,a2,a3,a4,a6]}{\\deg(\\phi)}{p\\div\\deg(\\phi)}"<<endl;
 while (n<limit) { n++;
#else
 while (cout<<"Enter level: ", cin>>n, n>0) {
#endif
 h1newforms nf(n,0,1);
 for(int xi=0; xi<nf.n1ds; xi++)
   { int i = xi;
#ifdef BOOKORDER
     i=booknumber0(n,i);
#endif
     int degphi = nf.nflist[i].degphi;
     int sfe = nf.nflist[i].sfe;
     Complex * periods;
     if(usedirect)
       {
	 periods_direct per(&(nf.nflist[i]));
	 periods = per.getperiods();
       }
     else
       {
	 periods_via_lfchi per(&(nf.nflist[i]));
	 periods = per.getperiods();
       }
     Complex w1 = periods[0], w2=periods[1],c4,c6;
     Complex tau=normalize(w1,w2);
     getc4c6(w1,w2,c4,c6);
     Integer ic4 = Iround(real(c4));
     Integer ic6 = Iround(real(c6));
#include "fixc6.cc"
     IntegerArray ai = getai(ic4,ic6);
     char acode; if (n<201) acode=antwerpletters[n][xi][0];
     cout << "\\entry{"<<n << "}{" << "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[xi] << "1";
     if(n<201) cout << " ("<<acode<<")";
     cout << "}";
     cout << "{["<<ai[0]<<","<<ai[1]<<","<<ai[2]<<","<<ai[3]<<","<<ai[4]<<"]}";
     cout << "{" << degphi;
     if(degphi>1)
       {
	 longlist pd = pdivs(degphi);
	 if(pd[0]<degphi)
	   {cout << " = ";
	    for(longvar pvar(pd); pvar.ok(); pvar++)
	      {
		long p = pvar.value(); int e = val(p,degphi);
		if(pvar.index>0) cout <<" \\cdot ";
		cout << p;
		if(e>1) cout << "^{"<<e<<"}";
	      }
	  }
       }
     cout << "}";
//     cout << "{" << pdivs(degphi) << "}";
     cout <<endl;
   }
}       // end of while()
//time(&stoptime);
//cout << "cpu time = " << (stoptime-starttime) << " seconds" << endl;
}       // end of main()


