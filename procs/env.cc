// env.cc:  test program for reading from an environment variable

#include <stdio.h>  // for sprintf()
#include <stdlib.h>  // for mkstemp()
//#include <unistd.h>  // for unlink() (not needed on linux)
#include <stream.h>
#include <streambuf.h>

int main()
{
  char* filename = new char[40];
  filename = "";
  char* tmpmatdir = getenv("TMPMATDIR");
  cout << "Value of tmpmatdir from environment = [" << tmpmatdir<<"]"<<endl;
  if (tmpmatdir==NULL) 
    {
      tmpmatdir=new char[4]; 
      sprintf(tmpmatdir,"/tmp");
    }
  cout << "Value of tmpmatdir = [" << tmpmatdir<<"]"<<endl;
  sprintf(filename,"%s%s",tmpmatdir,"/opmatXXXXXX");
  cout << "Before randomizing, *filename = [" << filename<<"]"<<endl;
  mkstemp(filename);
  cout << "After randomizing, *filename = [" << filename<<"]"<<endl;
}
