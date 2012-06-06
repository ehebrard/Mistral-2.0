#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <iomanip>

#include "flatzinc.hpp"

using namespace std;

int main(int argc, char *argv[])
{

  cout << "c Mistral-fzn" << endl;

  list<string> args(argv+1, argv+argc);

  if(args.empty()) {
    cout << "usage: minicsp-fz [options] input.fzn";
    return 1;
  }

  Solver s;
  double cpu_time = get_run_time() ;

  FlatZinc::Printer p;
  FlatZinc::FlatZincModel *fm = 0L;
  fm = parse(args.back(), s, p);
  if( !fm ) return 0;
  double parse_time = get_run_time() - cpu_time;


  //std::cout << s << std::endl;

  fm->run(cout , p);
  delete fm;

  // std::cout << "%sParse time:" << parse_time << std::endl
  // 	    << s.statistics << std::endl;
 
  return 0;
}
