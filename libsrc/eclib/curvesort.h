// FILE CURVESORT.H:  isogeny class id codes etc
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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

int booknumber0(int level, int form);  // permutes numbers starting from 0

const char alphabet[] = "abcdefghijklmnopqrstuvwxyz";

int booknumber(int level, int form);  // permutes numbers starting from 1

// new-new codes (from 01.08.05) are:
// a,b,...,z,ba,bb,...,bz,ca,cb,... etc., i.e. straight base 26 with
// digits a=0, b=1, ..., z=25

// Function to convert new code to integer (from 0) for any length of code.
// NB N=176400 has 516 < 26^2 newforms, with codes from a to vt!
int codeletter_to_int(string code);  // i counts from 0!

// Function to convert integer (from 0) to new code

string new_codeletter(int i);  // i counts from 0!

#ifdef USE_NEW_CODE_LETTERS
inline string codeletter(int i) {return new_codeletter(i);}
#else
inline string codeletter(int i) {return old_codeletter(i);}
#endif

