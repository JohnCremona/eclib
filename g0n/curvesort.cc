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

#define ORDER(level,reorder) case level: return reorder(form); break;

#define BA(x) (1-x)  // interchange 0 and 1
#define ACB(x) ((int[]){0,2,1})[x]
#define BAC(x) ((int[]){1,0,2})[x]
#define BCA(x) ((int[]){1,2,0})[x]
#define CAB(x) ((int[]){2,0,1})[x]
#define CBA(x) ((int[]){2,1,0})[x]
#define ABDC(x) ((int[]){0,1,3,2})[x]
#define ACBD(x) ((int[]){0,2,1,3})[x]
#define ACDB(x) ((int[]){0,2,3,1})[x]
#define ADBC(x) ((int[]){0,3,1,2})[x]
#define ADCB(x) ((int[]){0,3,2,1})[x]
#define BACD(x) ((int[]){1,0,2,3})[x]
#define BADC(x) ((int[]){1,0,3,2})[x]
#define BCAD(x) ((int[]){1,2,0,3})[x]
#define BCDA(x) ((int[]){1,2,3,0})[x]
#define BDAC(x) ((int[]){1,3,0,2})[x]
#define CABD(x) ((int[]){2,0,1,3})[x]
#define CBAD(x) ((int[]){2,1,0,3})[x]
#define CDAB(x) ((int[]){2,3,0,1})[x]
#define CDBA(x) ((int[]){2,3,1,0})[x]
#define DABC(x) ((int[]){3,0,1,2})[x]
#define DACB(x) ((int[]){3,0,2,1})[x]
#define DBAC(x) ((int[]){3,1,0,2})[x]
#define DCAB(x) ((int[]){3,2,0,1})[x]
#define AEBDC(x) ((int[]){0,4,1,3,2})[x]
#define AECBD(x) ((int[]){0,4,2,1,3})[x]
#define AEDCB(x) ((int[]){0,4,3,2,1})[x]
#define BCADE(x) ((int[]){1,2,0,3,4})[x]
#define BCAED(x) ((int[]){1,2,0,4,3})[x]
#define BEACD(x) ((int[]){1,4,0,2,3})[x]
#define BECAD(x) ((int[]){1,4,2,0,3})[x]
#define CEADB(x) ((int[]){2,4,0,3,1})[x]
#define CABDE(x) ((int[]){2,0,1,3,4})[x]
#define CEDAB(x) ((int[]){2,4,3,0,1})[x]
#define CEDBA(x) ((int[]){2,4,3,1,0})[x]
#define DAECB(x) ((int[]){3,0,4,2,1})[x]
#define DCAEB(x) ((int[]){3,2,0,4,1})[x]
#define DBCAE(x) ((int[]){3,1,2,0,4})[x]
#define EABDC(x) ((int[]){4,0,1,3,2})[x]
#define EACBD(x) ((int[]){4,0,2,1,3})[x]
#define EBCAD(x) ((int[]){4,1,2,0,3})[x]
#define EDACB(x) ((int[]){4,3,0,2,1})[x]
#define ADECFB(x) ((int[]){0,3,4,2,5,1})[x]
#define AECDBF(x) ((int[]){0,4,2,3,1,5})[x]
#define BADECF(x) ((int[]){1,0,3,4,2,5})[x]
#define BDAECF(x) ((int[]){1,3,0,4,2,5})[x]
#define BFADCE(x) ((int[]){1,5,0,3,2,4})[x]
#define CABEDF(x) ((int[]){2,0,1,4,3,5})[x]
#define CAFBDE(x) ((int[]){2,0,5,1,3,4})[x]
#define CEABDF(x) ((int[]){2,4,0,1,3,5})[x]
#define DABFEC(x) ((int[]){3,0,1,5,4,2})[x]
#define DCABEF(x) ((int[]){3,2,0,1,4,5})[x]
#define EABCDF(x) ((int[]){4,0,1,2,3,5})[x]
#define EBACDF(x) ((int[]){4,1,0,2,3,5})[x]
#define ECABDF(x) ((int[]){4,2,0,1,3,5})[x]
#define AEDBGFC(x) ((int[]){0,4,3,1,6,5,2})[x]
#define AFGDEBC(x) ((int[]){0,5,6,3,4,1,2})[x]
#define ECFBDGA(x) ((int[]){4,2,5,1,3,6,0})[x]
#define EFCAGDB(x) ((int[]){4,5,2,0,6,3,1})[x]
#define EFGCABD(x) ((int[]){4,5,6,2,0,1,3})[x]
#define FGACBED(x) ((int[]){5,6,0,2,1,4,3})[x]
#define FGBDACE(x) ((int[]){5,6,1,3,0,2,4})[x]
#define FGDBAEC(x) ((int[]){5,6,3,1,0,4,2})[x]
#define FGDBEAC(x) ((int[]){5,6,3,1,4,0,2})[x]
#define FGDEABC(x) ((int[]){5,6,3,4,0,1,2})[x]
#define AGBDEFHC(x) ((int[]){0,6,1,3,4,5,7,2})[x]
#define AFHCGDEB(x) ((int[]){0,5,7,2,6,3,4,1})[x]
#define BECADFGH(x) ((int[]){1,4,2,0,3,5,6,7})[x]
#define EGBAFHCD(x) ((int[]){4,6,1,0,5,7,2,3})[x]
#define GBEAHDFC(x) ((int[]){6,1,4,0,7,3,5,2})[x]

int booknumber0(int level, int form)  // permutes number starting from 0
{
  if(level<56) return form;
  if(level>450) return form;
switch (level) {
ORDER(56,BA)
ORDER(77,ACB)
ORDER(84,BA)
ORDER(99,ACBD)
ORDER(102,ACB)
ORDER(106,DACB)
ORDER(110,CBA)
ORDER(114,CAB)
ORDER(116,CAB)
ORDER(118,ACDB)
ORDER(120,BA)
ORDER(121,CABD)
ORDER(123,BA)
ORDER(124,BA)
ORDER(126,BA)
ORDER(128,ADCB)
ORDER(130,ACB)
ORDER(132,BA)
ORDER(136,BA)
ORDER(140,BA)
ORDER(141,EBCAD)
ORDER(142,EABDC)
ORDER(144,BA)
ORDER(147,ACB)
ORDER(150,CAB)
ORDER(153,ADCB)
ORDER(154,ACB)
ORDER(155,CBA)
ORDER(158,EACBD)
ORDER(162,ACBD)
ORDER(168,BA)
ORDER(170,DAECB) 
ORDER(171,ADBC)
ORDER(174,CEDBA)
ORDER(175,CAB)
ORDER(178,BA)
ORDER(182,CEADB)
ORDER(184,DABC)
ORDER(185,BCA)
ORDER(186,ACB)
ORDER(187,BA)
ORDER(189,ADBC)
ORDER(190,BAC)
ORDER(192,ACBD)
ORDER(195,CDBA)
ORDER(196,BA)
ORDER(198,CEDAB)
ORDER(200,BECAD)
ORDER(201,BCA)
ORDER(203,BCA)
ORDER(205,CAB)
ORDER(208,DABC)
ORDER(210,DBCAE)
ORDER(212,BA)
ORDER(214,DBAC)
ORDER(219,BCA)
ORDER(221,BA)
ORDER(222,EDACB)
ORDER(234,BCADE)
ORDER(235,CAB)
ORDER(236,BA)
ORDER(238,DAECB)
ORDER(240,BCAD)
ORDER(242,BA)
ORDER(245,CAB)
ORDER(246,EFCAGDB)
ORDER(249,BA)
ORDER(252,BA)
ORDER(254,DCAB)
ORDER(256,BACD)
ORDER(262,BA)
ORDER(264,DABC)
ORDER(267,BA)
ORDER(270,BCDA)
ORDER(272,ADBC)
ORDER(274,CBA)
ORDER(278,BA)
ORDER(285,CAB)
ORDER(286,BDAECF)
ORDER(288,AECBD)
ORDER(291,DBAC)
ORDER(294,EFGCABD)
ORDER(297,DACB)
ORDER(298,BA)
ORDER(300,ACDB)
ORDER(302,BAC)
ORDER(304,ECABDF)
ORDER(306,CBAD)
ORDER(312,CAFBDE)
ORDER(315,BA)
ORDER(318,DCAEB)
ORDER(320,EABCDF)
ORDER(322,BADC)
ORDER(324,ABDC)
ORDER(325,EABDC)
ORDER(326,ACB)
ORDER(330,AEDCB)
ORDER(333,DABC)
ORDER(336,CABEDF)
ORDER(338,ADECFB)
ORDER(339,BAC)
ORDER(342,FGDEABC)
ORDER(345,AECDBF)
ORDER(348,BDAC)
ORDER(350,BFADCE)
ORDER(352,CEABDF)
ORDER(354,DCABEF)
ORDER(360,BCAED)
ORDER(364,BA)
ORDER(366,FGBDACE)
ORDER(368,AEDBGFC)
ORDER(369,BA)
ORDER(370,ACBD)
ORDER(372,ACDB)
ORDER(378,GBEAHDFC)
ORDER(380,BA)
ORDER(381,BA)
ORDER(384,BECADFGH)
ORDER(387,CABDE)
ORDER(390,AFGDEBC)
ORDER(392,EBACDF)
ORDER(396,BAC)
ORDER(400,AFHCGDEB)
ORDER(402,ABDC)
ORDER(404,BA)
ORDER(405,DABFEC)
ORDER(406,ACDB)
ORDER(408,ADCB)
ORDER(410,BDAC)
ORDER(414,CDAB)
ORDER(418,BAC)
ORDER(423,ECFBDGA)
ORDER(425,ADBC)
ORDER(426,CAB)
ORDER(427,BCA)
ORDER(432,EGBAFHCD)
ORDER(434,AEBDC)
ORDER(435,BACD)
ORDER(437,BA)
ORDER(438,FGACBED)
ORDER(440,ADBC)
ORDER(441,BADECF)
ORDER(442,BEACD)
ORDER(446,ADCB)
ORDER(448,AGBDEFHC)
ORDER(450,FGDBEAC)

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

char* antwerpletters[201][17]=

{{},{},{},{},{},{},{},{},{},{},{},
/*11:*/{"BCA"},
/*12:*/{},
/*13:*/{},
/*14:*/{"CDEAFB"},
/*15:*/{"CEBFHGDA"},
/*16:*/{},
/*17:*/{"CBDA"},
/*18:*/{},
/*19:*/{"BCA"},
/*20:*/{"BADC"},
/*21:*/{"BDCAFE"},
/*22:*/{},
/*23:*/{},
/*24:*/{"BCDAFE"},
/*25:*/{},
/*26:*/{"BCA","DE"},
/*27:*/{"BDAC"},
/*28:*/{},
/*29:*/{},
/*30:*/{"ABCDEFGH"},
/*31:*/{},
/*32:*/{"BACD"},
/*33:*/{"BADC"},
/*34:*/{"ABCD"},
/*35:*/{"BCA"},
/*36:*/{"ABCD"},
/*37:*/{"A","CDB"},
/*38:*/{"DEC","AB"},
/*39:*/{"BCDA"},
/*40:*/{"BDAC"},
/*41:*/{},
/*42:*/{"ABCDFE"},
/*43:*/{"A"},
/*44:*/{"AB"},
/*45:*/{"ABDCEFHG"},
/*46:*/{"AB"},
/*47:*/{},
/*48:*/{"BDCAFE"},
/*49:*/{"ABCD"},
/*50:*/{"EFGH","ABCD"},
/*51:*/{"AB"},
/*52:*/{"BA"},
/*53:*/{"A"},
/*54:*/{"EFD","ACB"},
/*55:*/{"BDCA"},
/*56:*/{"CDEF","AB"},
/*57:*/{"E","BACD","FG"},
/*58:*/{"A","BC"},
/*59:*/{},
/*60:*/{},
/*61:*/{"A"},
/*62:*/{"ABCD"},
/*63:*/{"ABCDFE"},
/*64:*/{"BCDA"},
/*65:*/{"AB"},
/*66:*/{"ABCD","EFHG","IJLK"},
/*67:*/{"A"},
/*68:*/{},
/*69:*/{"AB"},
/*70:*/{"ABDC"},
/*71:*/{},
/*72:*/{"ABDCFE"},
/*73:*/{"BA"},
/*74:*/{},
/*75:*/{"AB","EFGHIJLK","CD"},
/*76:*/{"A"},
/*77:*/{"F","DEC","AB"},
/*78:*/{"ABCD"},
/*79:*/{"A"},
/*80:*/{"FEHG","BADC"},
/*81:*/{},
/*82:*/{"AB"},
/*83:*/{"A"},
/*84:*/{"CDEF","AB"},
/*85:*/{"AB"},
/*86:*/{},
/*87:*/{},
/*88:*/{"A"},
/*89:*/{"C","AB"},
/*90:*/{"MNOP","ABCD","EFGIHJLK"},
/*91:*/{"A","BCD"},
/*92:*/{"AB","C"},
/*93:*/{},
/*94:*/{"AB"},
/*95:*/{},
/*96:*/{"EFHG","ADBC"},
/*97:*/{},
/*98:*/{"BADCFE"},
/*99:*/{"AB","HIKJ","FG","CDE"},
/*100:*/{"ABCD"},
/*101:*/{"A"},
/*102:*/{"EF","GHJILK","ABCD"},
/*103:*/{},
/*104:*/{"A"},
/*105:*/{"ABDC"},
/*106:*/{"BC","A","EF","D"},
/*107:*/{},
/*108:*/{"AB"},
/*109:*/{"A"},
/*110:*/{"CD","AB","EF"},
/*111:*/{},
/*112:*/{"KL","ABDC","EFGHIJ"},
/*113:*/{"BA"},
/*114:*/{"ABCD","EF","GHJI"},
/*115:*/{"A"},
/*116:*/{"E","AB","DC"},
/*117:*/{"ABDC"},
/*118:*/{"A","BC","D","E"},
/*119:*/{},
/*120:*/{"EFHGJI","ABCD"},
/*121:*/{"HI","DE","FG","ABC"},
/*122:*/{"A"},
/*123:*/{"AB","C"},
/*124:*/{"BC","A"},
/*125:*/{},
/*126:*/{"ABCDEF","GHJILK"},
/*127:*/{},
/*128:*/{"CD","FE","AB","GH"},
/*129:*/{"E","BACD"},
/*130:*/{"EFGH","ABDC","JI"},
/*131:*/{"A"},
/*132:*/{"AB","CD"},
/*133:*/{},
/*134:*/{},
/*135:*/{"A","B"},
/*136:*/{"AB","CD"},
/*137:*/{},
/*138:*/{"EF","GHIJ","ABDC"},
/*139:*/{"A"},
/*140:*/{"AB","C"},
/*141:*/{"E","GF","ABCD","I","H"},
/*142:*/{"F","E","AB","CD","G"},
/*143:*/{"A"},
/*144:*/{"ABCD","EFGHJI"},
/*145:*/{"AB"},
/*146:*/{},
/*147:*/{"CDEFHG","IJ","AB"},
/*148:*/{"A"},
/*149:*/{},
/*150:*/{"ABCD","GHEF","IJKLMNOP"},
/*151:*/{},
/*152:*/{"A","B"},
/*153:*/{"C","AB","EFHG","D"},
/*154:*/{"CD","EFGH","AB"},
/*155:*/{"DE","AB","C"},
/*156:*/{"EF","ABCD"},
/*157:*/{},
/*158:*/{"E","D","HI","BCA","FG"},
/*159:*/{},
/*160:*/{"DC","AB"},
/*161:*/{"BACD"},
/*162:*/{"KL","GHIJ","ABDC","EF"},
/*163:*/{"A"},
/*164:*/{},
/*165:*/{},
/*166:*/{"A"},
/*167:*/{},
/*168:*/{"BACD","EFGH"},
/*169:*/{},
/*170:*/{"AB","HIJK","FG","DE","C"},
/*171:*/{"DEFG","ABC","IJ","H"},
/*172:*/{"AB"},
/*173:*/{},
/*174:*/{"IJ","GH","F","ABCD","E"},
/*175:*/{"BA","CDE","FG"},
/*176:*/{"C","DEF","AB"},
/*177:*/{},
/*178:*/{"AB","CD"},
/*179:*/{"A"},
/*180:*/{"ABCD"},
/*181:*/{},
/*182:*/{"EFGH","ABC","J","D","I"},
/*183:*/{},
/*184:*/{"C","B","DE","A"},
/*185:*/{"D","A","BC"},
/*186:*/{"D","BC","A"},
/*187:*/{"AB","C"},
/*188:*/{},
/*189:*/{"A","CDE","FGH","B"},
/*190:*/{"D","C","AB"},
/*191:*/{},
/*192:*/{"QRTS","ABDC","KLMNPO","EFHGJI"},
/*193:*/{},
/*194:*/{"AB"},
/*195:*/{"ABDCEFHG","I","K","J"},
/*196:*/{"AB","CD"},
/*197:*/{"A"},
/*198:*/{"IJLK","EFGH","MNOP","ABCD","QRST"},
/*199:*/{},
/*200:*/{"B","CD","GHJI","EF","A"}};

