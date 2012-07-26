#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <iomanip>

#include "flatzinc.hpp"
#include <set>
#include <map>
using namespace std;
#include <typeinfo>

#ifdef _VERIFICATION
void write_solution(FlatZinc::FlatZincModel *fm, string filename)
{
  //	cout << "set size " << fm->verif_constraints.size() << filename <<endl;
  ofstream myfile;
  filename.insert(0,"sol_" );
  myfile.open (filename.c_str());

  for (unsigned int i = 0; i < fm->verif_constraints.size() ; i++)
    {
      if (fm->verif_constraints[i].first != "int_in")
	{
	  myfile << "constraint " << fm->verif_constraints[i].first << "(" ;
	  int size = fm->verif_constraints[i].second.size();
	  for (unsigned int j = 0; j <size ; j++)
	    {

	      myfile << fm->verif_constraints[i].second[j].get_string() ;
	      if (j< (size -1))
		myfile << " , ";
	    }
	  myfile << ");" << endl;
	}
    }
  myfile <<"solve satisfy;" << endl;
  myfile.close();
  std::cout <<" c DONE" << endl;
}
#endif


int main(int argc, char *argv[])
{

  cout << " c Mistral-fzn" << endl;

  list<string> args(argv+1, argv+argc);

  if(args.empty())
    {
      cout << "usage: mistral-fz [options] input.fzn";
      return 1;
    }

  Solver s;
  double cpu_time = get_run_time() ;

#ifdef _VERBOSE_PARSER
  std::cout << " c Parse: ";
#endif

  FlatZinc::Printer p;
  FlatZinc::FlatZincModel *fm = 0L;
  fm = parse(args.back(), s, p);
  if( !fm ) return 0;
  double parse_time = get_run_time() - cpu_time;
#ifdef _VERBOSE_PARSER
  std::cout << std::endl;
#endif
  std::cout << " d PARSETIME " << parse_time << std::endl;

  map<string, string> options;

  // default value
  options["--var_heuristic"] = "dom/wdeg";
  options["--val_heuristic"] = "randminmax";
  options["--restart"] = "luby";
  options["--seed"] = "123456";
  options["--limit"] = "10";
  options["--verbose"] = "1";
  options["--rewrite"] = "0";


  string option_name = "error";
  //string option_value = "error";

  list<string>::iterator it;
  list<string>::iterator the_end = args.end();
  --the_end;
  for(it=args.begin(); it!=the_end; ++it)
    {
      if((*it)[0] == '-') option_name = (*it);
      else
	{
	  //option_value = (*it);
	  options[option_name] = (*it);
	  //std::cout << option_name << " <- " << (*it) << std::endl;
	}
    }

  map<string,string>::iterator ito;

  double cutoff = atof(options["--limit"].c_str()) - parse_time;
  std::cout << " d CUTOFF " << cutoff << std::endl;

  fm->set_strategy(options["--var_heuristic"], options["--val_heuristic"], options["--restart"]);
  s.initialise_random_seed(atoi(options["--seed"].c_str()));
  s.set_time_limit(cutoff);
  s.parameters.verbosity = atoi(options["--verbose"].c_str());
  fm->set_rewriting(atoi(options["--rewrite"].c_str()));

  fm->run(cout , p);

#ifdef _VERIFICATION
  write_solution(fm, args.back());
#endif

  delete fm;
  return 0;
}
