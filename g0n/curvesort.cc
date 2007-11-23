// FILE CURVESORT.CC:  isogeny class id codes etc
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

#define USE_NEW_CODE_LETTERS

int BA[] = {1,0};
int ACB[] = {0,2,1};
int BAC[] = {1,0,2};
int BCA[] = {1,2,0};
int CAB[] = {2,0,1};
int CBA[] = {2,1,0};
int ABDC[] = {0,1,3,2};
int ACBD[] = {0,2,1,3};
int ACDB[] = {0,2,3,1};
int ADBC[] = {0,3,1,2};
int ADCB[] = {0,3,2,1};
int BACD[] = {1,0,2,3};
int BADC[] = {1,0,3,2};
int BCAD[] = {1,2,0,3};
int BCDA[] = {1,2,3,0};
int BDAC[] = {1,3,0,2};
int CABD[] = {2,0,1,3};
int CBAD[] = {2,1,0,3};
int CDAB[] = {2,3,0,1};
int CDBA[] = {2,3,1,0};
int DABC[] = {3,0,1,2};
int DACB[] = {3,0,2,1};
int DBAC[] = {3,1,0,2};
int DCAB[] = {3,2,0,1};
int AEBDC[] = {0,4,1,3,2};
int AECBD[] = {0,4,2,1,3};
int AEDCB[] = {0,4,3,2,1};
int BCADE[] = {1,2,0,3,4};
int BCAED[] = {1,2,0,4,3};
int BEACD[] = {1,4,0,2,3};
int BECAD[] = {1,4,2,0,3};
int CEADB[] = {2,4,0,3,1};
int CABDE[] = {2,0,1,3,4};
int CEDAB[] = {2,4,3,0,1};
int CEDBA[] = {2,4,3,1,0};
int DAECB[] = {3,0,4,2,1};
int DCAEB[] = {3,2,0,4,1};
int DBCAE[] = {3,1,2,0,4};
int EABDC[] = {4,0,1,3,2};
int EACBD[] = {4,0,2,1,3};
int EBCAD[] = {4,1,2,0,3};
int EDACB[] = {4,3,0,2,1};
int ADECFB[] = {0,3,4,2,5,1};
int AECDBF[] = {0,4,2,3,1,5};
int BADECF[] = {1,0,3,4,2,5};
int BDAECF[] = {1,3,0,4,2,5};
int BFADCE[] = {1,5,0,3,2,4};
int CABEDF[] = {2,0,1,4,3,5};
int CAFBDE[] = {2,0,5,1,3,4};
int CEABDF[] = {2,4,0,1,3,5};
int DABFEC[] = {3,0,1,5,4,2};
int DCABEF[] = {3,2,0,1,4,5};
int EABCDF[] = {4,0,1,2,3,5};
int EBACDF[] = {4,1,0,2,3,5};
int ECABDF[] = {4,2,0,1,3,5};
int AEDBGFC[] = {0,4,3,1,6,5,2};
int AFGDEBC[] = {0,5,6,3,4,1,2};
int ECFBDGA[] = {4,2,5,1,3,6,0};
int EFCAGDB[] = {4,5,2,0,6,3,1};
int EFGCABD[] = {4,5,6,2,0,1,3};
int FGACBED[] = {5,6,0,2,1,4,3};
int FGBDACE[] = {5,6,1,3,0,2,4};
int FGDBAEC[] = {5,6,3,1,0,4,2};
int FGDBEAC[] = {5,6,3,1,4,0,2};
int FGDEABC[] = {5,6,3,4,0,1,2};
int AGBDEFHC[] = {0,6,1,3,4,5,7,2};
int AFHCGDEB[] = {0,5,7,2,6,3,4,1};
int BECADFGH[] = {1,4,2,0,3,5,6,7};
int EGBAFHCD[] = {4,6,1,0,5,7,2,3};
int GBEAHDFC[] = {6,1,4,0,7,3,5,2};

int booknumber0(int level, int form)  // utes number starting from 0
{
  if(level<56) return form;
  if(level>450) return form;
switch (level) {
case 56: return BA[form]; break;
case 77: return ACB[form]; break;
case 84: return BA[form]; break;
case 99: return ACBD[form]; break;
case 102: return ACB[form]; break;
case 106: return DACB[form]; break;
case 110: return CBA[form]; break;
case 114: return CAB[form]; break;
case 116: return CAB[form]; break;
case 118: return ACDB[form]; break;
case 120: return BA[form]; break;
case 121: return CABD[form]; break;
case 123: return BA[form]; break;
case 124: return BA[form]; break;
case 126: return BA[form]; break;
case 128: return ADCB[form]; break;
case 130: return ACB[form]; break;
case 132: return BA[form]; break;
case 136: return BA[form]; break;
case 140: return BA[form]; break;
case 141: return EBCAD[form]; break;
case 142: return EABDC[form]; break;
case 144: return BA[form]; break;
case 147: return ACB[form]; break;
case 150: return CAB[form]; break;
case 153: return ADCB[form]; break;
case 154: return ACB[form]; break;
case 155: return CBA[form]; break;
case 158: return EACBD[form]; break;
case 162: return ACBD[form]; break;
case 168: return BA[form]; break;
case 170: return DAECB[form]; break;
case 171: return ADBC[form]; break;
case 174: return CEDBA[form]; break;
case 175: return CAB[form]; break;
case 178: return BA[form]; break;
case 182: return CEADB[form]; break;
case 184: return DABC[form]; break;
case 185: return BCA[form]; break;
case 186: return ACB[form]; break;
case 187: return BA[form]; break;
case 189: return ADBC[form]; break;
case 190: return BAC[form]; break;
case 192: return ACBD[form]; break;
case 195: return CDBA[form]; break;
case 196: return BA[form]; break;
case 198: return CEDAB[form]; break;
case 200: return BECAD[form]; break;
case 201: return BCA[form]; break;
case 203: return BCA[form]; break;
case 205: return CAB[form]; break;
case 208: return DABC[form]; break;
case 210: return DBCAE[form]; break;
case 212: return BA[form]; break;
case 214: return DBAC[form]; break;
case 219: return BCA[form]; break;
case 221: return BA[form]; break;
case 222: return EDACB[form]; break;
case 234: return BCADE[form]; break;
case 235: return CAB[form]; break;
case 236: return BA[form]; break;
case 238: return DAECB[form]; break;
case 240: return BCAD[form]; break;
case 242: return BA[form]; break;
case 245: return CAB[form]; break;
case 246: return EFCAGDB[form]; break;
case 249: return BA[form]; break;
case 252: return BA[form]; break;
case 254: return DCAB[form]; break;
case 256: return BACD[form]; break;
case 262: return BA[form]; break;
case 264: return DABC[form]; break;
case 267: return BA[form]; break;
case 270: return BCDA[form]; break;
case 272: return ADBC[form]; break;
case 274: return CBA[form]; break;
case 278: return BA[form]; break;
case 285: return CAB[form]; break;
case 286: return BDAECF[form]; break;
case 288: return AECBD[form]; break;
case 291: return DBAC[form]; break;
case 294: return EFGCABD[form]; break;
case 297: return DACB[form]; break;
case 298: return BA[form]; break;
case 300: return ACDB[form]; break;
case 302: return BAC[form]; break;
case 304: return ECABDF[form]; break;
case 306: return CBAD[form]; break;
case 312: return CAFBDE[form]; break;
case 315: return BA[form]; break;
case 318: return DCAEB[form]; break;
case 320: return EABCDF[form]; break;
case 322: return BADC[form]; break;
case 324: return ABDC[form]; break;
case 325: return EABDC[form]; break;
case 326: return ACB[form]; break;
case 330: return AEDCB[form]; break;
case 333: return DABC[form]; break;
case 336: return CABEDF[form]; break;
case 338: return ADECFB[form]; break;
case 339: return BAC[form]; break;
case 342: return FGDEABC[form]; break;
case 345: return AECDBF[form]; break;
case 348: return BDAC[form]; break;
case 350: return BFADCE[form]; break;
case 352: return CEABDF[form]; break;
case 354: return DCABEF[form]; break;
case 360: return BCAED[form]; break;
case 364: return BA[form]; break;
case 366: return FGBDACE[form]; break;
case 368: return AEDBGFC[form]; break;
case 369: return BA[form]; break;
case 370: return ACBD[form]; break;
case 372: return ACDB[form]; break;
case 378: return GBEAHDFC[form]; break;
case 380: return BA[form]; break;
case 381: return BA[form]; break;
case 384: return BECADFGH[form]; break;
case 387: return CABDE[form]; break;
case 390: return AFGDEBC[form]; break;
case 392: return EBACDF[form]; break;
case 396: return BAC[form]; break;
case 400: return AFHCGDEB[form]; break;
case 402: return ABDC[form]; break;
case 404: return BA[form]; break;
case 405: return DABFEC[form]; break;
case 406: return ACDB[form]; break;
case 408: return ADCB[form]; break;
case 410: return BDAC[form]; break;
case 414: return CDAB[form]; break;
case 418: return BAC[form]; break;
case 423: return ECFBDGA[form]; break;
case 425: return ADBC[form]; break;
case 426: return CAB[form]; break;
case 427: return BCA[form]; break;
case 432: return EGBAFHCD[form]; break;
case 434: return AEBDC[form]; break;
case 435: return BACD[form]; break;
case 437: return BA[form]; break;
case 438: return FGACBED[form]; break;
case 440: return ADBC[form]; break;
case 441: return BADECF[form]; break;
case 442: return BEACD[form]; break;
case 446: return ADCB[form]; break;
case 448: return AGBDEFHC[form]; break;
case 450: return FGDBEAC[form]; break;

 default: return form;
}
  return form;  //default case for levels not yet sorted manually
}


const char alphabet[] = "abcdefghijklmnopqrstuvwxyz";

int booknumber(int level, int form)  // permutes number starting from 1
{
  return 1+booknumber0(level,form-1);
}

// Old codes were:  A,B,...,Z,AA,BB,...,ZZ,AAA,BBB,... etc

// Function to convert old code to integer (from 0)

int old_codeletter_to_int(char* code)  // i counts from 0!
{
  int i = code[0]-'A';
  int n=1; while(code[n]!='\0') n++;
  return 26*(n-1)+i;
}

// Function to convert integer (from 0) to old code

void old_codeletter(int i, char* code, int width=0)
{
  int n=width;    // pads string to this width with blanks
  code[n]='\0';
  while (n) code[--n]=' ';

  int nc = i%26;
  char c = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[nc];
  n = 1 + (i-nc)/26;
  if(width==0) code[n]='\0';
  while (n) code[--n]=c;
  /*
  int j = old_codeletter_to_int(code);
  if(i==j) return;
  cout<<i<<" -> "<<code<<" -> "<<j<<endl;
  */
}

// old-new codes were:  a,b,...,z,aa,ab,...,az,ba,bb,... etc

// Function to convert old-new code to integer (from 0)

int old_new_codeletter_to_int(char* code)  // i counts from 0!
{
  int b = code[0]-'a';
  if(code[1]=='\0')  return b;
  int a = code[1]-'a';
  return 26*(b+1)+a;
}

// Function to convert integer (from 0) to old-new code

void old_new_codeletter(int i, char* code)  // i counts from 0!
{
  int b = i%26;
  int a = (i-b)/26;
  if (a==0) {code[0]=alphabet[b]; code[1]='\0';}
  else {code[0]=alphabet[a-1]; code[1]=alphabet[b]; code[2]='\0';}
  /*
  int j = codeletter_to_int(code);
  if(i==j) return;
  cout<<i<<" -> "<<code<<" -> "<<j<<endl;
  */
}

// new-new codes (from 01.08.05) are:  a,b,...,z,ba,bb,...,bz,ca,cb,... etc.,  i.e. straight base 26 with digits a=0, b=1, ..., z=25

// Function to convert new code to integer (from 0)

int codeletter_to_int(char* code)  // i counts from 0!
{
  int b = code[0]-'a';
  if(code[1]=='\0')  return b;
  int a = code[1]-'a';
  return 26*b+a;
}

// Function to convert integer (from 0) to new code

void new_codeletter(int i, char* code)  // i counts from 0!
{
  int b = i%26;
  int a = (i-b)/26;
  if (a==0) {code[0]=alphabet[b]; code[1]='\0';}
  else {code[0]=alphabet[a]; code[1]=alphabet[b]; code[2]='\0';}
  /*
  int j = codeletter_to_int(code);
  if(i==j) return;
  cout<<i<<" -> "<<code<<" -> "<<j<<endl;
  */
}

#ifdef USE_NEW_CODE_LETTERS
inline void codeletter(int i, char* code) {return new_codeletter(i,code);}
#else
inline void codeletter(int i, char* code) {return old_codeletter(i,code);}
#endif

