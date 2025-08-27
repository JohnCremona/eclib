// FILE CURVESORT.CC:  isogeny class id codes etc
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

#include <sstream>
#include <algorithm> // for reverse()
#include "eclib/curvesort.h"

int curvesort_BA[] = {1,0};
int curvesort_ACB[] = {0,2,1};
int curvesort_BAC[] = {1,0,2};
int curvesort_BCA[] = {1,2,0};
int curvesort_CAB[] = {2,0,1};
int curvesort_CBA[] = {2,1,0};
int curvesort_ABDC[] = {0,1,3,2};
int curvesort_ACBD[] = {0,2,1,3};
int curvesort_ACDB[] = {0,2,3,1};
int curvesort_ADBC[] = {0,3,1,2};
int curvesort_ADCB[] = {0,3,2,1};
int curvesort_BACD[] = {1,0,2,3};
int curvesort_BADC[] = {1,0,3,2};
int curvesort_BCAD[] = {1,2,0,3};
int curvesort_BCDA[] = {1,2,3,0};
int curvesort_BDAC[] = {1,3,0,2};
int curvesort_CABD[] = {2,0,1,3};
int curvesort_CBAD[] = {2,1,0,3};
int curvesort_CDAB[] = {2,3,0,1};
int curvesort_CDBA[] = {2,3,1,0};
int curvesort_DABC[] = {3,0,1,2};
int curvesort_DACB[] = {3,0,2,1};
int curvesort_DBAC[] = {3,1,0,2};
int curvesort_DCAB[] = {3,2,0,1};
int curvesort_AEBDC[] = {0,4,1,3,2};
int curvesort_AECBD[] = {0,4,2,1,3};
int curvesort_AEDCB[] = {0,4,3,2,1};
int curvesort_BCADE[] = {1,2,0,3,4};
int curvesort_BCAED[] = {1,2,0,4,3};
int curvesort_BEACD[] = {1,4,0,2,3};
int curvesort_BECAD[] = {1,4,2,0,3};
int curvesort_CEADB[] = {2,4,0,3,1};
int curvesort_CABDE[] = {2,0,1,3,4};
int curvesort_CEDAB[] = {2,4,3,0,1};
int curvesort_CEDBA[] = {2,4,3,1,0};
int curvesort_DAECB[] = {3,0,4,2,1};
int curvesort_DCAEB[] = {3,2,0,4,1};
int curvesort_DBCAE[] = {3,1,2,0,4};
int curvesort_EABDC[] = {4,0,1,3,2};
int curvesort_EACBD[] = {4,0,2,1,3};
int curvesort_EBCAD[] = {4,1,2,0,3};
int curvesort_EDACB[] = {4,3,0,2,1};
int curvesort_ADECFB[] = {0,3,4,2,5,1};
int curvesort_AECDBF[] = {0,4,2,3,1,5};
int curvesort_BADECF[] = {1,0,3,4,2,5};
int curvesort_BDAECF[] = {1,3,0,4,2,5};
int curvesort_BFADCE[] = {1,5,0,3,2,4};
int curvesort_CABEDF[] = {2,0,1,4,3,5};
int curvesort_CAFBDE[] = {2,0,5,1,3,4};
int curvesort_CEABDF[] = {2,4,0,1,3,5};
int curvesort_DABFEC[] = {3,0,1,5,4,2};
int curvesort_DCABEF[] = {3,2,0,1,4,5};
int curvesort_EABCDF[] = {4,0,1,2,3,5};
int curvesort_EBACDF[] = {4,1,0,2,3,5};
int curvesort_ECABDF[] = {4,2,0,1,3,5};
int curvesort_AEDBGFC[] = {0,4,3,1,6,5,2};
int curvesort_AFGDEBC[] = {0,5,6,3,4,1,2};
int curvesort_ECFBDGA[] = {4,2,5,1,3,6,0};
int curvesort_EFCAGDB[] = {4,5,2,0,6,3,1};
int curvesort_EFGCABD[] = {4,5,6,2,0,1,3};
int curvesort_FGACBED[] = {5,6,0,2,1,4,3};
int curvesort_FGBDACE[] = {5,6,1,3,0,2,4};
int curvesort_FGDBAEC[] = {5,6,3,1,0,4,2};
int curvesort_FGDBEAC[] = {5,6,3,1,4,0,2};
int curvesort_FGDEABC[] = {5,6,3,4,0,1,2};
int curvesort_AGBDEFHC[] = {0,6,1,3,4,5,7,2};
int curvesort_AFHCGDEB[] = {0,5,7,2,6,3,4,1};
int curvesort_BECADFGH[] = {1,4,2,0,3,5,6,7};
int curvesort_EGBAFHCD[] = {4,6,1,0,5,7,2,3};
int curvesort_GBEAHDFC[] = {6,1,4,0,7,3,5,2};

int booknumber0(int level, int form)  // permutes numbers starting from 0
{
  if(level<56) return form;
  if(level>450) return form;
switch (level) {
case 56: return curvesort_BA[form]; break;
case 77: return curvesort_ACB[form]; break;
case 84: return curvesort_BA[form]; break;
case 99: return curvesort_ACBD[form]; break;
case 102: return curvesort_ACB[form]; break;
case 106: return curvesort_DACB[form]; break;
case 110: return curvesort_CBA[form]; break;
case 114: return curvesort_CAB[form]; break;
case 116: return curvesort_CAB[form]; break;
case 118: return curvesort_ACDB[form]; break;
case 120: return curvesort_BA[form]; break;
case 121: return curvesort_CABD[form]; break;
case 123: return curvesort_BA[form]; break;
case 124: return curvesort_BA[form]; break;
case 126: return curvesort_BA[form]; break;
case 128: return curvesort_ADCB[form]; break;
case 130: return curvesort_ACB[form]; break;
case 132: return curvesort_BA[form]; break;
case 136: return curvesort_BA[form]; break;
case 140: return curvesort_BA[form]; break;
case 141: return curvesort_EBCAD[form]; break;
case 142: return curvesort_EABDC[form]; break;
case 144: return curvesort_BA[form]; break;
case 147: return curvesort_ACB[form]; break;
case 150: return curvesort_CAB[form]; break;
case 153: return curvesort_ADCB[form]; break;
case 154: return curvesort_ACB[form]; break;
case 155: return curvesort_CBA[form]; break;
case 158: return curvesort_EACBD[form]; break;
case 162: return curvesort_ACBD[form]; break;
case 168: return curvesort_BA[form]; break;
case 170: return curvesort_DAECB[form]; break;
case 171: return curvesort_ADBC[form]; break;
case 174: return curvesort_CEDBA[form]; break;
case 175: return curvesort_CAB[form]; break;
case 178: return curvesort_BA[form]; break;
case 182: return curvesort_CEADB[form]; break;
case 184: return curvesort_DABC[form]; break;
case 185: return curvesort_BCA[form]; break;
case 186: return curvesort_ACB[form]; break;
case 187: return curvesort_BA[form]; break;
case 189: return curvesort_ADBC[form]; break;
case 190: return curvesort_BAC[form]; break;
case 192: return curvesort_ACBD[form]; break;
case 195: return curvesort_CDBA[form]; break;
case 196: return curvesort_BA[form]; break;
case 198: return curvesort_CEDAB[form]; break;
case 200: return curvesort_BECAD[form]; break;
case 201: return curvesort_BCA[form]; break;
case 203: return curvesort_BCA[form]; break;
case 205: return curvesort_CAB[form]; break;
case 208: return curvesort_DABC[form]; break;
case 210: return curvesort_DBCAE[form]; break;
case 212: return curvesort_BA[form]; break;
case 214: return curvesort_DBAC[form]; break;
case 219: return curvesort_BCA[form]; break;
case 221: return curvesort_BA[form]; break;
case 222: return curvesort_EDACB[form]; break;
case 234: return curvesort_BCADE[form]; break;
case 235: return curvesort_CAB[form]; break;
case 236: return curvesort_BA[form]; break;
case 238: return curvesort_DAECB[form]; break;
case 240: return curvesort_BCAD[form]; break;
case 242: return curvesort_BA[form]; break;
case 245: return curvesort_CAB[form]; break;
case 246: return curvesort_EFCAGDB[form]; break;
case 249: return curvesort_BA[form]; break;
case 252: return curvesort_BA[form]; break;
case 254: return curvesort_DCAB[form]; break;
case 256: return curvesort_BACD[form]; break;
case 262: return curvesort_BA[form]; break;
case 264: return curvesort_DABC[form]; break;
case 267: return curvesort_BA[form]; break;
case 270: return curvesort_BCDA[form]; break;
case 272: return curvesort_ADBC[form]; break;
case 274: return curvesort_CBA[form]; break;
case 278: return curvesort_BA[form]; break;
case 285: return curvesort_CAB[form]; break;
case 286: return curvesort_BDAECF[form]; break;
case 288: return curvesort_AECBD[form]; break;
case 291: return curvesort_DBAC[form]; break;
case 294: return curvesort_EFGCABD[form]; break;
case 297: return curvesort_DACB[form]; break;
case 298: return curvesort_BA[form]; break;
case 300: return curvesort_ACDB[form]; break;
case 302: return curvesort_BAC[form]; break;
case 304: return curvesort_ECABDF[form]; break;
case 306: return curvesort_CBAD[form]; break;
case 312: return curvesort_CAFBDE[form]; break;
case 315: return curvesort_BA[form]; break;
case 318: return curvesort_DCAEB[form]; break;
case 320: return curvesort_EABCDF[form]; break;
case 322: return curvesort_BADC[form]; break;
case 324: return curvesort_ABDC[form]; break;
case 325: return curvesort_EABDC[form]; break;
case 326: return curvesort_ACB[form]; break;
case 330: return curvesort_AEDCB[form]; break;
case 333: return curvesort_DABC[form]; break;
case 336: return curvesort_CABEDF[form]; break;
case 338: return curvesort_ADECFB[form]; break;
case 339: return curvesort_BAC[form]; break;
case 342: return curvesort_FGDEABC[form]; break;
case 345: return curvesort_AECDBF[form]; break;
case 348: return curvesort_BDAC[form]; break;
case 350: return curvesort_BFADCE[form]; break;
case 352: return curvesort_CEABDF[form]; break;
case 354: return curvesort_DCABEF[form]; break;
case 360: return curvesort_BCAED[form]; break;
case 364: return curvesort_BA[form]; break;
case 366: return curvesort_FGBDACE[form]; break;
case 368: return curvesort_AEDBGFC[form]; break;
case 369: return curvesort_BA[form]; break;
case 370: return curvesort_ACBD[form]; break;
case 372: return curvesort_ACDB[form]; break;
case 378: return curvesort_GBEAHDFC[form]; break;
case 380: return curvesort_BA[form]; break;
case 381: return curvesort_BA[form]; break;
case 384: return curvesort_BECADFGH[form]; break;
case 387: return curvesort_CABDE[form]; break;
case 390: return curvesort_AFGDEBC[form]; break;
case 392: return curvesort_EBACDF[form]; break;
case 396: return curvesort_BAC[form]; break;
case 400: return curvesort_AFHCGDEB[form]; break;
case 402: return curvesort_ABDC[form]; break;
case 404: return curvesort_BA[form]; break;
case 405: return curvesort_DABFEC[form]; break;
case 406: return curvesort_ACDB[form]; break;
case 408: return curvesort_ADCB[form]; break;
case 410: return curvesort_BDAC[form]; break;
case 414: return curvesort_CDAB[form]; break;
case 418: return curvesort_BAC[form]; break;
case 423: return curvesort_ECFBDGA[form]; break;
case 425: return curvesort_ADBC[form]; break;
case 426: return curvesort_CAB[form]; break;
case 427: return curvesort_BCA[form]; break;
case 432: return curvesort_EGBAFHCD[form]; break;
case 434: return curvesort_AEBDC[form]; break;
case 435: return curvesort_BACD[form]; break;
case 437: return curvesort_BA[form]; break;
case 438: return curvesort_FGACBED[form]; break;
case 440: return curvesort_ADBC[form]; break;
case 441: return curvesort_BADECF[form]; break;
case 442: return curvesort_BEACD[form]; break;
case 446: return curvesort_ADCB[form]; break;
case 448: return curvesort_AGBDEFHC[form]; break;
case 450: return curvesort_FGDBEAC[form]; break;

 default: return form;
}
  return form;  //default case for levels not yet sorted manually
}

int booknumber(int level, int form)  // permutes numbers starting from 1
{
  return 1+booknumber0(level,form-1);
}

// new-new codes (from 01.08.05) are:
// a,b,...,z,ba,bb,...,bz,ca,cb,... etc., i.e. straight base 26 with
// digits a=0, b=1, ..., z=25

// Function to convert new code to integer (from 0) for any length of code.
int codeletter_to_int(string code)  // i counts from 0!
{
  int n=0;
  std::for_each(code.begin(), code.end(), [&n] ( const char& c ) { n*=26; n+=(c-'a');});
  return n;
}

// Function to convert integer (from 0) to new code

string new_codeletter(int i)  // i counts from 0!
{
  if (i==0) return string("a"); // special case -- otherwise leading
                                // a's are omitted
  stringstream code;
  int n = i;
  while (n)
  {
    std::div_t x = div(n,26);
    code << alphabet[x.rem];
    n = x.quot;
  }
  string res = code.str();
  reverse(res.begin(),res.end());
  return res;
}

