#include "kinematics.cpp"
#include <string>
#include <iostream>
using namespace std;

int main(int argc, char *argv[]) {
  
  string filepath = Form("%s",argv[1]);

  string outdir = Form("%s",argv[2]);

  string tester = Form("%s", argv[3]);
  //write "test" for the third argument to only test 100 entries for debugging

  string clas_acceptance = Form("%s", argv[4]);
  //enter acc to run with clas acceptance or anything else to run without
  cout << filepath << " " <<outdir << endl;
 
  func(filepath, outdir, tester, clas_acceptance);
  
}

