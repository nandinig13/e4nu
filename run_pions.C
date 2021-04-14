#include "nandini.cpp"
#include <string>
#include <iostream>
using namespace std;

int main(int argc, char *argv[]) {
  
  string filepath = Form("%s",argv[1]);

  string outdir = Form("%s",argv[2]);

  cout << filepath << " " <<outdir << endl;
 
  func(filepath, outdir);
  
}

