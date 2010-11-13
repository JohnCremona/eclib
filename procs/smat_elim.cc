// smat_elim.cc: implementation of class smat_elim
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2005 John Cremona
// 
// This file is part of the mwrank package.
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
 
//  implements structured modular elimination
//  original written by Luiz Figueiredo
//  this version by JEC December 2005

#include <time.h>
#include <sys/times.h>

//#define TRACE_ELIM

void smat_elim::init_elim( )
{
#ifdef TRACE_ELIM
  cout<<"smat_elim::init()"<<endl; //cout<<" smat =\n"<<(smat)(*this)<<endl;
#endif
  reduce_mod_p();
  rank = 0;
  position = vector<int>(nro+1,-1);  // row numbers start at 1
  elim_col = vector<int>(nco+1,-1);  // col numbers start at 1
  elim_row = vector<int>(nro+1,0);   //  ditto
  column = vector<std::set<int> >(nco+1); //  ditto
  light_col_flag = vector<int>(nco+1,0);
  ech_form = red_ech_form = 0;
  int row;
  vector<svec>::const_iterator ai;
  map<int,scalar>::const_iterator aij;
  for( (ai = rows.begin())++, row=1; ai!=rows.end(); ai++, row++ ) 
    for( aij = ai->begin(); aij!=ai->end(); aij++)
      {
	if((aij->first<1)||(aij->first>nco)) 	
	  cout<<"inserting "<<row<<" into column "<<aij->first<<endl;
	column[aij->first].insert(row);
      }
  //  cout<<"At end of init(), columns are: \n";
  //  for(int l = 1; l <= nco; l++ ) 
  //    cout<<l<<": "<<column[l]<<"\n";

  nrows_left=nro;
  for( int row=1; row<=nro; row++ ) 
    if( rows[row].size() ==0 ) 
      {position[row]=0; nrows_left--;}
  ncols_left=nco;
  vector<std::set<int> >::iterator ci;
  for((ci=column.begin())++; ci!=column.end(); ci++)
    if(ci->size()==0) ncols_left--;
#ifdef TRACE_ELIM
  cout<<"After init:"; show_progress();
#endif
}

// This just does steps 0,1,2,3,4,5,6 in turn:
void smat_elim::echelon_form( )
{
  init_elim();
#ifdef TRACE_ELIM
  long starttime,stoptime;
  starttime=clock();
  int flag=0, flag2=1;
  int pop=get_population(*this);
  double density = (double)pop/(nco*nro);
  cout<<"Starting sparse elimination: "<<nro<<" rows, "<<nco<<" columns, "<<pop<<" entries (density = "<<density<<")\n";
  cout<<"Starting steps 1-2-3..."<<flush;
#endif
  int prog=nrows_left+ncols_left+1, pass=0;
  while(nrows_left+ncols_left<prog)
    {
      prog=nrows_left+ncols_left;
      pass++;
#ifdef TRACE_ELIM
      cout<<"Pass "<<pass<<endl;
#endif
  for(int fc=1; fc<=2; fc++)
    {
      elim_light_rows(fc);
#ifdef TRACE_ELIM
      cout << "after elim_light_rows("<<fc<<")" << endl; 
      cout << "# of rows eliminated so far: " << get_rank( ) << endl;
      if( flag ) cout << " matrix: " << (*this) << endl;
      if(flag>1) cout << " = " << this->as_mat();
      if( flag2 ) cout<<"("<<get_population( *this )<<" entries)"<<endl;
#endif
      //if(fc==1)    
  elim_light_cols(fc);
#ifdef TRACE_ELIM
      cout << "after elim_light_cols("<<fc<<")" << endl; 
      cout << "# of rows eliminated so far: " << get_rank( ) << endl;
      if( flag ) cout << " matrix: " << (*this) << endl;
      if(flag>1) cout << " = " << this->as_mat();
      if( flag2 ) cout<<"("<<get_population( *this )<<" entries)"<<endl;
#endif
    }
    }
#ifdef TRACE_ELIM
  stoptime=clock();
  cout<<"finished: ";
  cout << "cpu time (steps 1-3) = " 
       << ((double)(stoptime-starttime)/CLOCKS_PER_SEC) 
       << " seconds"<<endl;
  pop=get_population(*this);
  density = (double)pop/(nco*nro);
  cout<<"Rank so far = "<<rank<<"; "
      <<pop<<" entries (density = "<<density<<")\n";
  cout<<"Starting (new) step 4..."<<flush;
  //cout<<"Starting (old) step 4..."<<flush;
  starttime=clock();
#endif
  step4new();
  //step4();
#ifdef TRACE_ELIM
  stoptime=clock();
  cout<<"finished: ";
  cout << "cpu time (steps 4) = " 
       << ((double)(stoptime-starttime)/CLOCKS_PER_SEC) 
       << " seconds"<<endl;
  pop=get_population(*this);
  density = (double)pop/(nco*nro);
  cout<<"Rank so far = "<<rank<<"; "
      <<pop<<" entries (density = "<<density<<")\n";
  cout<<"Starting step 5 (all remaining elimination)..."<<flush;
  starttime=clock();
#endif
  step5( );
#ifdef TRACE_ELIM
  stoptime=clock();
  cout<<"finished: ";
  cout << "cpu time (step 5) = " 
       << ((double)(stoptime-starttime)/CLOCKS_PER_SEC) 
       << " seconds"<<endl;
  //  if(!check_echelon())
  //    cout<<"Echelon check fails!"<<endl;
  pop=get_population(*this);
  density = (double)pop/(nco*nro);
  cout<<"Rank so far = "<<rank<<"; "
      <<pop<<" entries (density = "<<density<<")\n";
#endif
  ech_form=1;
}

//#define TRACE_ELIM
void smat_elim::reduced_echelon_form( )
{
#ifdef TRACE_ELIM
  cout<<"Starting step 6 (back-substitution)..."<<flush;
  long starttime,stoptime;
  starttime=clock();
#endif
  step6( );
#ifdef TRACE_ELIM
  stoptime=clock();
  cout<<"finished: ";
  cout << "cpu time (back-substitution) = " 
       << ((double)(stoptime-starttime)/CLOCKS_PER_SEC) 
       << " seconds"<<endl;
  //  if(!check_red_echelon())
  //    cout<<"Reduced echelon check fails!"<<endl;
  int pop=get_population(*this);
  double density = (double)pop/(nco*nro);
  cout<<"Rank so far = "<<rank<<"; "
      <<pop<<" entries (density = "<<density<<")\n";
#endif
  red_ech_form=1;
}

#undef TRACE_ELIM
// function to call when the (row,col) entry is the unique nonzero
// entry in its column: the entry in elim_col[] shows which row
// eliminated this col; the entry in position[] shows that this row
// has been eliminated, and which is its pivotal column; the entry in
// elim_row[] records which row is the (cumulative) rank'th o be
// eliminated; the set column[col] is cleared (emptied) to show that
// this col has been eliminated
void smat_elim::eliminate( int& row, int& col ) //1<=col<=nco;
{
  rank++;
#ifdef TRACE_ELIM
  //  cout<<"Eliminating (row,col) = ("<<row<<","<<col<<"); rank now "<<rank<<endl;
#else
  //  cout<<"."<<flush;
#endif
  elim_col[ col ]  = row;
  position[ row ]  = col;
  elim_row[ rank ] = row;
  column[col].clear();
  nrows_left--;
  ncols_left--;
}

// eliminate (row,col) for rows of weight at most fr
void smat_elim::elim_light_rows(int fr)
{
  if(finished()) 
    {
#ifdef TRACE_ELIM
      cout<<"elim_light_rows("<<fr<<") quitting, already finished"<<endl;
#endif
      return;
    }
  int row, col, wt;

  // Find the rows not yet eliminated with at most fr entries,:

  //  light_rows.clear();  
  for( row=1; row<=nro; row++ ) 
    if( position[row]==-1 )
      {
	wt=rows[row].size();
	if( (wt>0) && (wt <= fr) )
	  light_rows.push(row);
      }
#ifdef TRACE_ELIM
  //  cout<<"In elim_light_rows("<<fr<<"), light_rows = "<<light_rows<<endl;
#endif
  // Loop through these (but we may add to the list as we go along):

  while(!light_rows.empty())
    {
      row=light_rows.front(); light_rows.pop();
      // In case this row has been eliminated since we started this loop:
      if( position[row] != -1 ) continue;
      if( rows[row].size() == 0 ) {position[row]=0; continue;}
      col = rows[row].first_index();
      clear_col ( row, col, fr );      
      eliminate( row, col );
    }
}

// eliminate (row,col) for columns of weight <= fc
void smat_elim::elim_light_cols(int fc)
{
  if(finished()) 
    {
#ifdef TRACE_ELIM
      cout<<"elim_light_cols("<<fc<<") quitting, already finished"<<endl;
#endif
      return;
    }
  int wt, row, col;

  // Find the cols with at most fc entries:

  //  light_cols.clear();  
  for( col = 1; col <= nco; col++ )
    {
      wt=column[col].size();
      if((  wt>0 ) && ( wt <= fc ) ) light_cols.push(col);
    }
    
  // Loop through these (but we may add to the list as we go along):
  
  std::set<int>::const_iterator ri;
  while(!light_cols.empty())
    {
      col=light_cols.front(); light_cols.pop();
      if(column[col].size() ==0 ) continue;
      row = *(column[col].begin());
      wt=rows[row].size();
      for((ri=column[col].begin())++; ri!=column[col].end(); ri++)
	if(rows[*ri].size()<wt) 
	  wt=rows[row=*ri].size();
      clear_col( row, col, 0, fc );
      eliminate( row, col );
    }
}

//#define TRACE_ELIM
// eliminate (row,col) for "light" columns 
void smat_elim::step4 ( )
{
  if(finished()) 
    {
#ifdef TRACE_ELIM
      cout<<"step4() quitting, already finished"<<endl;
#endif
      return;
    }
  vector<int>::iterator l;
  vector<std::set<int> >::iterator c;
  map<int,scalar>::const_iterator rij;
  int M, maxcolwt, wt, r, row, col;

#undef TRACE_ELIM
  // Find max weight of an active col:  
  maxcolwt=0;
  for( col = 1; col <= nco; col++ )
    {
      wt=column[col].size();
      if(wt>maxcolwt)  maxcolwt=wt;
    }
#ifdef TRACE_ELIM
  cout<<"Max column weight = "<<maxcolwt<<endl;
#endif
  // Find heavy cols to mark as inactive, successively marking more

  for( M = maxcolwt; M >= 3; M--)
    {
#ifdef TRACE_ELIM
      cout<<"(r="<<rank<<",M="<<M<<")"<<flush;
#endif
      // tag columns which are `light' : not yet eliminated but no
      // more than M entries
      for( (l = light_col_flag.begin())++, (c=column.begin())++; 
	   l!= light_col_flag.end(); 
	   l++, c++ ) 
	{
	  wt = c->size();
	  *l = ( 0 < wt && wt <= M );  // 0 or 1
	}
      
      // eliminate all (row,col) where col is the unique light column in row

      int fr;
      for(fr=1; fr<=2; fr++) {
      do {	  // find such a row not yet eliminated whose weight
		  // in the flagged columns is 1 or 2:
	  for( r = 1, row = 0; (r<=nro)&&(row==0); r++ ) 
	    if(position[r] == -1)
	      {wt=get_weight(r);
	      if(wt==fr)
		{
		  row = r;  // we found one...
		  // find # of a light col cutting this row
		  col=0;
		  for(rij=rows[row].begin(), col=0; 
		      (col==0)&&(rij!=rows[row].end());rij++)
		    if(light_col_flag[rij->first]) col=(rij->first);
		  if(col==0) cout<<"Problem"<<endl;
		  clear_col(row,col, 0, 0, M);
		  eliminate( row, col );
		}
	      }
      }
      while (row!=0);
      }
    } // end of M loop
}
//#undef TRACE_ELIM

// eliminate (row,col) for "light" columns 
int smat_elim::step4finished()
{
  for( int col = 1; col <= nco; col++ )
    if((light_col_flag[col])&&(column[col].size()>0)) 
      return 0;
  return 1;
}

//#define TRACE_ELIM
void smat_elim::step4new ( )
{
  if(finished()) 
    {
#ifdef TRACE_ELIM
      cout<<"step4() quitting, already finished"<<endl;
#endif
      return;
    }
#undef TRACE_ELIM
  vector<int>::iterator l;
  vector<std::set<int> >::iterator c;
  map<int,scalar>::const_iterator rij;
  int maxcolwt, wt, row, col;

  // Find max weight of an active col, and mark all cols as light
  // (later some will by unmarked):
  maxcolwt=0;
  for( col = 1; col <= nco; col++ )
    {
      wt=column[col].size();
      if(wt>maxcolwt)  maxcolwt=wt;
      light_col_flag[col]=1;
    }
#ifdef TRACE_ELIM
  cout<<"Max column weight = "<<maxcolwt<<endl;
#endif

  // Follow strategy of Pomerance & Smith: mark the 5% heaviest
  // columns, eliminate (row,col) which only intersect 1 light row or
  // col, then successivle mark another 0.1% columns and continue.

  int fivepercent = (long)(0.05*nco); if(fivepercent==0) fivepercent=1;
  int point1percent = (long)(0.001*nco); if(point1percent==0) point1percent=1;
  int nhcol=fivepercent, nhcol0=0;
#ifdef TRACE_ELIM
  cout<<"#Heavy columns intialliy = "<<fivepercent
      <<" with increment "<<point1percent<<endl;
#endif

  // Parameter tweaking place here:  having (nhcol<nco) as the
  // stopping condition gets the most out of this stage, though the
  // last few rounds can be slow.  But having (nhcol<nco/10) makes the
  // final elimination stage take a lot longer too.  With
  // (nhcol<nco/8) and an example of size 30000, initially 8 entries
  // per row, I found that the final elimination stage did not take
  // any longer.

  int nco8=nco; //nco>>3;
  while((nhcol<nco8)&&(!step4finished())) {

#ifdef TRACE_ELIM
    //  cout<<"Using "<<nhcol<<" heavy columns..."<<flush;
    //  cout<<"."<<flush;
#endif

  for(wt=maxcolwt; (nhcol0<nhcol)&&(wt>0); wt--)
    for(col=1; (nhcol0<nhcol)&&(col<=nco); col++)
      if(light_col_flag[col])
	if((int)column[col].size()==wt)
	  {
	    light_col_flag[col]=0;
	    nhcol0++;
	  }
  
  // eliminate all (row,col) where col is the unique light column in row

  int count=0;
  
  for( row=1; row<=nro; row++ ) 
    if( position[row]==-1 )
      if(get_weight(row)==1) 
	light_rows.push(row);
  
  while(!light_rows.empty())
    {
      row=light_rows.front(); light_rows.pop();
      // In case this row has been eliminated since we started this loop:
      if( position[row] != -1 ) continue;
      if( rows[row].size() == 0 ) {position[row]=0; continue;}
      if( get_weight(row) != 1 ) continue;
      // find unique light pivotal column in this row:
      for(rij=rows[row].begin(), col=0; 
	  (col==0)&&(rij!=rows[row].end());rij++)
	if(light_col_flag[rij->first]) col=(rij->first);
      if(col==0) cout<<"Problem"<<endl;
      clear_col ( row, col, 0,0,0,1 );      
      eliminate( row, col );
      count++;
    }
  
#ifdef TRACE_ELIM
  if(count)
    {
      cout<<"Eliminated "<<count<<" rows using "<<nhcol<<" heavy columns";
  int p=0;
  for( int col = 1; col <= nco; col++ )
    if((light_col_flag[col])) p+=column[col].size();
  cout<<", remaining live population =  "<<p<<endl; 
    }
#endif 
  nhcol += point1percent;
  }
}

//#undef TRACE_ELIM
//#define TRACE_ELIM
// all remaining elimination
void smat_elim::step5 ( )
{
  if(finished()) 
    {
#ifdef TRACE_ELIM
      cout<<"step5() quitting, already finished"<<endl;
#endif
      return;
    }
  int row, col, wt;
  // This structure, which will be sorted, allows us to loop over the
  // columns in order of their weight:
  multimap<int,int> col_wts;
  for( col = 1; col <= nco; col++ )
    {
      wt=column[col].size();
      if( wt>0 ) col_wts.insert(pair<int,int>(wt,col));
    }
  multimap<int,int>::iterator cols;
  std::set<int>::const_iterator ci;
  for( cols=col_wts.begin(); cols!=col_wts.end(); cols++ )
    {
      col=cols->second;
      if(column[col].size()>0)
	{
	  // find row of least weight in this col:
	  ci = column[col].begin();
	  row = *ci; 
          /*
          wt=rows[row].size();
	  for(; (wt>1)&&(ci!=column[col].end()); ci++)
	    if(rows[*ci].size() < wt) 
	      {
		row = *ci; wt=rows[*ci].size();
	      }
	  //	  cout<<wt<<")"<<flush;
          */
	  clear_col( row, col );
	  eliminate( row, col );
	}
    }
}
  
#undef TRACE_ELIM
// Back substitution
void smat_elim::step6 ( )
{
#ifdef TRACE_ELIM
  //  cout<<"In step6(): m=\n"<<(*this).as_mat()<<endl;
#endif
  int i, col, col2, row, row2;
  scalar a;
  map<int,scalar>::const_iterator aij;
  map<int, pair<int,scalar> > elim_list;
  map<int, pair<int,scalar> >::iterator ei;
  int ci,cc=0; // count entries to be eliminated
  for( i=rank; i>0; i--)
    { 
      row = elim_row[i];
      for(aij=rows[row].begin(); aij!=rows[row].end(); aij++)
	{
	  ci=elim_col[aij->first];
	  if((ci!=-1)&&(ci!=row))
	    cc++;
	}
    }
#ifdef TRACE_ELIM
  cout<<": "<<cc<<" entries to be eliminated..."<<endl;
#endif
  // traverse rows in reverse order of their elimination
  for( i=rank; cc&&(i>0); i--)
    { 
#ifdef TRACE_ELIM
      cout<<"Back-subst: "<<i<<" rows remaining...";
      cout<<"Population = "<<get_population(*this)<<endl;
      //      cout<<"m=\n"<<(*this).as_mat()<<endl;
#endif
      row = elim_row[i];
      col  = position[row];
//Go along the row looking for entries in pivotal columns which we
//then eliminate; when we do so we restart at the beginning of the row
//since the pivotal columns are not in natural order:
      elim_list.clear();
      for(aij=rows[row].begin(); aij!=rows[row].end(); )
	{
	  col2=aij->first;  a=aij->second; aij++;
	  row2 = elim_col[col2]; // = which row will eliminate this entry
	  if((row2==-1) // entry not in a pivotal column -- ignore
	     ||(row2==row)) // entry is this row's pivot -- ignore
	    continue;
	  elim_list[col2]=pair<int,scalar>(row2,-a);
	}
      ci=elim_list.size(); cc-=ci;
      if(ci)
	{
#ifdef TRACE_ELIM
	  cout<<"--eliminating "<<ci<<" entries, "<<cc<<" more to go..."<<endl;
#endif
	  for(ei=elim_list.begin(); ei!=elim_list.end(); ei++)
	    rows[row].add_scalar_times_mod_p(rows[(ei->second).first],(ei->second).second);
	}
    }
}

// Given row# row with entry 1 in column# col, clears the rest of this
// column by subtracting suitable mutliples of the row

// fr ("row-flag", default 0): this row has weight 1 or 2; light_rows
// is a list of such rows, which should be updated if necessary

// fc ("col-flag", default 0): this col has weight 1 or 2; light_cols
// is a list of such cols, which should be updated if necessary

// M (weight threshold, default 0): this is a "light" col with at most
// M entries; these are indicated by flag array light_cols (global to
// class), which should be updated if necessary

void 
smat_elim::clear_col( int row, int col, int fr, int fc, int M, int frl, int fcl)
{
#ifdef TRACE_ELIM
  //  cout<<"Entering clear_col ("<<row<<","<<col<<") "<<endl;
  //  cout<<"column[col] = "<<column[col]<<endl;
#endif
  map<int,scalar>::const_iterator ri;
  // for each changed row, these sets contain column #s added/deleted
  std::set<int> cols_out;
  std::set<int> cols_in;
  std::set<int>::iterator coli, ci;
  int wt, col2, row2;

// rescale row #row so its (non-zero) col'th entry is 1 mod p
  scalar s = rows[row].elem(col);
  s = xmod0(s);
  if(s==0) 
    {
      cout<<"Error in smat_elim::clear_col()!\nEntry #"<<col
	  <<" in row "<<row<<" = "<<rows[row]<<" is zero"<<endl;
      abort();
    }
  if(s!=1)
    {
      s = invmod0(s);
      rows[row].mult_by_scalar_mod_p(s);
    }

  // eliminate col from other rows cutting col 
  for( coli=column[col].begin(); coli!=column[col].end(); coli++ ) 
    {
      row2 = *coli;
      if( row2 == row ) continue;
      cols_in.clear();
      cols_out.clear();
      //      cout<<rows[row2]<<endl;
      rows[row2].add_scalar_times_mod_p(rows[row], -rows[row2].elem(col), 
					cols_in, cols_out);
      //      cout<<rows[row2]<<endl;
      // update column lists for OTHER columns (this col will be
      // cleared anyway, and deleting the current entry will
      // invalidate the iterator coli)
      for( ci=cols_in.begin(); ci!=cols_in.end(); ci++ ) 
	{if((*ci)!=col) {column[*ci].insert(row2);}}
      for( ci=cols_out.begin(); ci!=cols_out.end(); ci++ ) 
	{if((*ci)!=col) 
	  {
	    column[*ci].erase(row2); 
	    if(column[*ci].size()==0) {ncols_left--;}

	  }
	}

      if( fr ) //  check condition for rows
	{
	  wt = rows[row2].size();
	  if( wt == 0 )   {position[row2] = 0; nrows_left--;}
	  else if(wt <= fr) light_rows.push(row2);  
	}      
      if( frl ) //  check condition for rows
	if( get_weight(row2) == 1 ) 
	  light_rows.push(row2);  
    }
#ifdef TRACE_ELIM
  //  if(fr) cout<<"light_rows now = "<<light_rows<<endl;
#endif
  
  // update column lists for the cols in this row
  //  cout<<"updating column lists for cols in row "<<row<<endl;
  for( ri=rows[row].begin(); ri!=rows[row].end(); ri++ ) 
    {
      col2 = ri->first; 
      if(col2!=col) column[col2].erase(row);
      if(column[col2].size()==0) {ncols_left--;}

      if( fc ) // add column# col2 to list if its weight is <=fc
	{
	  wt = column[col2].size();
	  if( 0<wt && wt <= fc ) light_cols.push(col2);
	}        
#if(0)
      if( M ) // set flag of column# col2 if its weight is <=M
	{
	  wt = column[col2].size();
	  light_col_flag[col2] = ( (0<wt) && (wt<=M) );  // 0 or 1
	}
#endif
    }
#ifdef TRACE_ELIM
  //  if(fc) cout<<"light_cols now = "<<light_cols<<endl;
  //  if(M) cout<<"light_col_flag now = "<<light_col_flag<<endl;
#endif
}

// compute the number of "light"  columns intersecting row# row
int smat_elim::get_weight( int row ) 
{
  int wt;
  map<int,scalar>::const_iterator ri;
  for(wt=0, ri=rows[row].begin(); ri!=rows[row].end(); ri++)
    if(light_col_flag[ri->first]) wt++;
  return wt;
}

// check that we do have echelon form
int smat_elim::check_echelon()
{
  //  cout<<"Checking echelon condition..."<<endl;
  //  cout<<"position = "<<position<<endl;
  //  cout<<"elim_row = "<<elim_row<<endl;
  //  cout<<"elim_col = "<<elim_col<<endl;
  int i, r;
  for(i=1; i<=nro; i++)
    {
      if(position[i]==-1) return 0;
      if((position[i]==0)!=(rows[i].size()==0)) return 0;
    }
  //  cout<<"ech: ok so far"<<endl;
  vector<int> elim_row_inv(nro+1,-1);
  for(i=1; i<=rank; i++)
    {
      elim_row_inv[elim_row[i]]=i;
    }
  map<int,scalar>::const_iterator aij;
  for( i=rank; i>0; i--)
    { 
      r = elim_row[i];
      for(aij=rows[r].begin(); aij!=rows[r].end(); aij++)
	{
	  int ci=elim_col[aij->first];
	  if(ci!=-1)
	    if(elim_row_inv[ci]<i) return 0;
	}
    }
  return 1;
}

// check that we do have reduced echelon form
int smat_elim::check_red_echelon()
{
  //  cout<<"Checking reduced echelon condition..."<<endl;
  //  cout<<"position = "<<position<<endl;
  //  cout<<"elim_row = "<<elim_row<<endl;
  //  cout<<"elim_col = "<<elim_col<<endl;
  int i, r;
  for(i=1; i<=nro; i++) 
    {
      if(position[i]==-1) return 0;
      if((position[i]==0)!=(rows[i].size()==0)) return 0;
    }
  //  cout<<"red-ech: ok so far"<<endl;
  map<int,scalar>::const_iterator aij;
  for( i=rank; i>0; i--)
    { 
      r = elim_row[i];
      for(aij=rows[r].begin(); aij!=rows[r].end(); aij++)
	{
	  int ci=elim_col[aij->first];
	  if((ci!=-1)&&(ci!=r))
	    return 0;
	}
    }
  return 1;
}

//#define TRACE_ELIM
// extract kernel from an smat after reducing it to reduced echelon form:
smat smat_elim::oldkernel( vec& pc, vec& npc)
{
  int i,j;
#ifdef TRACE_ELIM
  cout<<"Starting echelon_form()..."<<flush;
  long starttime,stoptime;
#endif
  echelon_form();
  reduced_echelon_form();
  reduce_mod_p();
#ifdef TRACE_ELIM
  cout<<"finished echelon_form()"<<endl;
  //  cout<<"Now A = \n"<<(*this)<<endl;
  //  cout<<"After elimination, A = \n"<<(this->as_mat())<<endl;
  //  cout<<"Checking reduced echelon form: "<<check_red_echelon()<<endl;
  //  cout<<"position = "<<position<<endl;
  //  cout<<"elim_row = "<<elim_row<<endl;
  //  cout<<"elim_col = "<<elim_col<<endl;
#endif

  int nullity = nco - rank;

  // set-up vecs pc & npc 
  pc.init( rank );
  npc.init( nullity );
  vector<int> ind(nco+1,0);
#ifdef TRACE_ELIM
  //  cout<<"Setting up pc and npc..."<<flush;
#endif
  int ny = 0, k = 0;
  for( i = 1; i <= nco; i++ )
    if( elim_col[i] != -1 ) // then this is a pivotal column
      {pc.set(++k,i); ind[i]=k;}
    else 
      {npc.set(++ny,i); ind[i]=ny;}

  if(ny!=nullity) 
    cout<<"Error: nullity = "<<nullity<<" but "
	<<ny<<" non-pivotal columns"<<endl;
  if(k!=rank) 
    cout<<"Error: rank = "<<rank<<" but "
	<<k<<" pivotal columns"<<endl;
  
#ifdef TRACE_ELIM
  //  cout<<"pc = "<<pc<<endl;
  //  cout<<"npc = "<<npc<<endl;
#endif


  // create kernel basis 

#ifdef TRACE_ELIM
  starttime=clock();
  cout<<"creating kernel basis of dimension "<<ny<<"..."<<flush;
#endif
  smat kerbasis( nco,nullity );
  int ri,rk,ci,ck;
  map<int,scalar>::const_iterator rik;

  // First enter the identity matrix in rows indexed by npc:

  for(j=1; j<=nullity; j++)
    kerbasis.rows[npc[j]].entries.insert(pair<int,scalar>(j,1));

  // Next enter the rest, by rows (because of the smat data structure)
  // though by columns would be easier to think about:

  for(i=1; i<=rank; i++)
    {
      ri=elim_row[i];
      ci=position[ri];
      rk=pc[ind[ci]]; // this is the row number we'll be entering data into
      insert_iterator<map<int,scalar> > rr(kerbasis.rows[rk].entries,kerbasis.rows[rk].entries.begin());
      for(rik=rows[ri].entries.begin(); rik!=rows[ri].entries.end(); rik++)
	{
	  ck=rik->first;
	  if(ck == ci) continue;
	  *rr++=pair<int,scalar>(ind[ck],-(rik->second));
	}
    }
#ifdef TRACE_ELIM
  cout<<"done."<<endl;
  stoptime=clock();
  cout << "cpu time for kernel (old method) = " 

       << ((double)(stoptime-starttime)/CLOCKS_PER_SEC) 
       << " seconds"<<endl;
#endif
  
  return kerbasis;
}

// This version works on an echelon form which is *not* necessarily in
// reduced echelon form (avoiding the need for the back-substitution
// step6()) 
smat smat_elim::kernel( vec& pc, vec& npc)
{
  int i,j,k,l,ny;
  scalar a;
#ifdef TRACE_ELIM
  cout<<"Starting echelon_form()..."<<flush;
  long starttime,stoptime;
#endif
  echelon_form();
  reduce_mod_p();
#ifdef TRACE_ELIM
  cout<<"finished echelon_form()"<<endl;
  //  cout<<"Now A = \n"<<(*this)<<endl;
  //  cout<<"After elimination, A = \n"<<(this->as_mat())<<endl;
  //  cout<<"Checking echelon form: "<<check_echelon()<<endl;
  //  cout<<"position = "<<position<<endl;
  //  cout<<"elim_row = "<<elim_row<<endl;
  //  cout<<"elim_col = "<<elim_col<<endl;
#endif

  int nullity = nco - rank;

  // set-up vecs pc & npc 
  pc.init( rank );
  npc.init( nullity );
  vector<int> ind(nco+1,0);
#ifdef TRACE_ELIM
  //  cout<<"Setting up pc and npc..."<<flush;
#endif
  ny = 0;
  k = 0;
  for( i = 1; i <= nco; i++ )
    if( elim_col[i] != -1 ) // then this is a pivotal column
      {pc.set(++k,i); ind[i]=k;}
    else 
      {npc.set(++ny,i); ind[i]=ny;}

  if(ny!=nullity) 
    cout<<"Error: nullity = "<<nullity<<" but "
	<<ny<<" non-pivotal columns"<<endl;
  if(k!=rank) 
    cout<<"Error: rank = "<<rank<<" but "
	<<k<<" pivotal columns"<<endl;
  
#ifdef TRACE_ELIM
  //  cout<<"pc = "<<pc<<endl;
  //  cout<<"npc = "<<npc<<endl;
#endif


  // create kernel basis 


#ifdef TRACE_ELIM
  starttime=clock();
  cout<<"creating kernel basis of dimension "<<ny<<"..."<<flush;
#endif
  smat kerbasis( nco,nullity );
  int ri,ci;
  map<int,scalar>::const_iterator rik;

  // First enter the identity matrix in rows indexed by npc:

  for(k=1; k<=nullity; k++)
    kerbasis.rows[npc[k]].entries.insert(pair<int,scalar>(k,1));

  // Next enter the rest, by rows

  for(i=rank; i>0; i--)     // order *does* matter
    {
      ri=elim_row[i];
      ci=position[ri]; // this is the row of kerbasis we'll be
                       // entering data into: elim_col[ci]==ri
      svec& rowci = kerbasis.rows[ci];
      for(rik=rows[ri].entries.begin(); rik!=rows[ri].entries.end(); rik++)
	{
	  j=rik->first; a=rik->second;
	  l = elim_col[j];
	  if(l==-1) rowci.sub_mod_p(ind[j],a); 
	  else if(l!=ri) rowci.add_scalar_times_mod_p(kerbasis.rows[j],-a);
	}
      kerbasis.rows[ci].reduce_mod_p();
    }
#ifdef TRACE_ELIM
  cout<<"done."<<endl;
  stoptime=clock();
  cout << "cpu time for kernel (new method) = " 

       << ((double)(stoptime-starttime)/CLOCKS_PER_SEC) 
       << " seconds"<<endl;
#endif
  
  return kerbasis;
}

int rank(smat& sm)
{
  smat_elim sme(sm);
  sme.echelon_form();
  return sme.get_rank();
}

ssubspace::ssubspace(int n) 
:pivots(iota((scalar)n)),basis(sidmat((scalar)n))
{}

ssubspace::ssubspace(const smat& b, const vec& p)
:pivots(p),basis(b)
{}

ssubspace::ssubspace(const ssubspace& s) 
:pivots(s.pivots),basis(s.basis) 
{}

// destructor -- no need to do anything as components have their own
ssubspace::~ssubspace() 
{}

// assignment
void ssubspace::operator=(const ssubspace& s) 
{
  pivots=s.pivots; 
  basis=s.basis; 
}

// Definitions of nonmember, nonfriend operators and functions:

ssubspace combine(const ssubspace& s1, const ssubspace& s2)
{ 
  return ssubspace(mult_mod_p(s1.basis,s2.basis,BIGPRIME),s1.pivots[s2.pivots]);
}
 
smat restrict(const smat& m, const ssubspace& s)
{ 
  return mult_mod_p(m.select_rows(pivots(s)),basis(s),BIGPRIME);
}
 
ssubspace kernel(const smat& sm)
{
  vec pivs, npivs;
  //  cout<<"In kernel()"<<endl;
  smat kern = smat_elim(sm).kernel(npivs,pivs);
  return ssubspace(kern,pivs);
}
 
ssubspace eigenspace(const smat& m1, scalar lambda)
{
  //  cout<<"In eigenspace(), lambda = "<<lambda<<endl;
  smat_elim m(m1); 
  //  cout<<"Finished assigning m"<<endl;
  m.sub_mod_p(lambda); 
  //  cout<<"Finished subtracting lambda"<<endl;
  vec pivs, npivs;
  smat kern = m.kernel(npivs,pivs);
  //  cout<<"Finished finding kernel"<<endl;
  return ssubspace(kern,pivs);
}
 
ssubspace subeigenspace(const smat& m1, scalar l, const ssubspace& s)
{
  return combine(s,eigenspace(restrict(m1,s), l));
}

ssubspace make1d(const vec& bas, long&piv) 
// make a 1-D ssubspace with basis bas
{
  smat tbasis(1,dim(bas));
  svec sbas(bas);
  tbasis.setrow(1,sbas);
  vec pivs(1); // initialised to 0
  pivs[1]=sbas.first_index();
  piv=sbas.elem(pivs[1]);
  return ssubspace(transpose(tbasis),pivs);
}



#undef TRACE_ELIM
