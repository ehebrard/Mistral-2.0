#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <iomanip>

#include "flatzinc.hpp"

using namespace std;

int main(int argc, char *argv[])
{

  cout << " c Mistral-fzn" << endl;

  list<string> args(argv+1, argv+argc);

  if(args.empty()) {
    cout << "usage: mistral-fz [options] input.fzn";
    return 1;
  }

  Solver s;
  double cpu_time = get_run_time() ;

  FlatZinc::Printer p;
  FlatZinc::FlatZincModel *fm = 0L;
  fm = parse(args.back(), s, p);
  if( !fm ) return 0;
  double parse_time = get_run_time() - cpu_time;


  map<string, string> options;

  // default value
  options["--var_heuristic"] = "dom/wdeg";
  options["--val_heuristic"] = "randminmax";
  options["--restart"] = "luby";
  options["--seed"] = "123456";


  string option_name = "error";
  //string option_value = "error";

  list<string>::iterator it;
  list<string>::iterator the_end = args.end();
  --the_end;
  for(it=args.begin(); it!=the_end; ++it) {
    if((*it)[0] == '-') option_name = (*it);
    else {
      //option_value = (*it);
      options[option_name] = (*it);
      //std::cout << option_name << " <- " << (*it) << std::endl;
    }
  }

  map<string,string>::iterator ito;

  // // show content:
  // for ( ito=options.begin() ; ito != options.end(); ito++ )
  //   std::cout << (*ito).first << " => " << (*ito).second << std::endl;



  fm->set_strategy(options["--var_heuristic"], options["--val_heuristic"], options["--restart"]);
  s.initialise_random_seed(atoi(options["--seed"].c_str()));



  //string restart = options["--heuristic"];


  fm->run(cout , p);
  delete fm;

  // std::cout << "%sParse time:" << parse_time << std::endl
  // 	    << s.statistics << std::endl;
 
  return 0;
}
