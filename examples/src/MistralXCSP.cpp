#include <XCSP3CoreParser.h>
#include "XCSP3MistralCallbacks.hpp"


using namespace XCSP3Core;


int main(int argc,char **argv) {
	
	Solver solver;
	
   XCSP3MistralCallbacks cb(solver); // my interface between the parser and the solver

   if(argc!=2) 
     throw std::runtime_error("usage: ./test xcsp3instance.xml");
  try
  {
    XCSP3CoreParser parser(&cb);
    parser.parse(argv[1]); // fileName is a string
  }
  catch (exception &e)
  {
    cout.flush();
    cerr << "\n\tUnexpected exception :\n";
    cerr << "\t" << e.what() << endl;
    exit(1);
  }

  return 0;
}

