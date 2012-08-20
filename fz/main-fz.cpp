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
  if (fm->finished())
    {
      //cout << fm->verif_constraints.size() << filename <<endl;
      unsigned int __size=filename.size();
      for (__size; __size>0; __size--)
	if (filename[__size-1] == '/')
	  break;
      filename.insert(__size,"sol_" );

      ofstream myfile;
      myfile.open (filename.c_str());

      //cout << filename <<endl;
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
}
#endif


int main(int argc, char *argv[])
{
#ifdef _FLATZINC_OUTPUT
	cout << "%";
#endif
  cout << " c Mistral-fzn" << endl;

  map<string, string> options;

  // default value
  options["--var_heuristic"] = "dom/wdeg";
  options["--val_heuristic"] = "guided";
  options["--restart"] = "luby";
  options["--seed"] = "123456";
  options["--limit"] = "10";
  options["--verbose"] = "0";
  options["--rewrite"] = "0";


  options["a"] = "0";
  options["p"] = "1";
  options["f"] = "1";
  



  list<string> args(argv+1, argv+argc);

  if(args.empty())
    {
	#ifdef _FLATZINC_OUTPUT
	  cout << "%";
	#endif
      cout << " usage: mistral-fz [options] input.fzn\n% options are:\n";

      map<string, string>::iterator it;
      for(it=options.begin(); it!=options.end(); ++it) {
	cout << "%  " << it->first << endl;
      }
      cout << endl;
      return 1;
    }

  Solver s;
  double cpu_time = get_run_time() ;

#ifdef _VERBOSE_PARSER
  std::cout << " c Parse: ";
#endif

  FlatZinc::Printer p;
  FlatZinc::FlatZincModel *fm = 0L;

  //FlatZinc::SolutionPrinter *sp = 

  fm = parse(args.back(), s, p);
  s.add(new FlatZinc::SolutionPrinter(&p, fm, &s));

  if( !fm ) return 0;
  double parse_time = get_run_time() - cpu_time;
#ifdef _VERBOSE_PARSER
  std::cout << std::endl;
#endif
#ifdef _FLATZINC_OUTPUT
	cout << "%";
#endif
  std::cout << " d PARSETIME " << parse_time << std::endl;

 

  string option_name = "error";
  //string option_value = "error";

  list<string>::iterator it;
  list<string>::iterator the_end = args.end();
  --the_end;
  int phase = 0;
  for(it=args.begin(); it!=the_end; ++it)
    {
      if((*it)[0] == '-') {
	if(phase) {
	  options[option_name] = "1";
	  //std::cout << "1" << std::endl;
	}
	option_name = (*it);
	//std::cout << option_name << " <- ";
	phase = 1;
      } else {
	//option_value = (*it);
	options[option_name] = (*it);
	//std::cout << (*it) << std::endl;
	phase = 0;
      }
    }
  if(phase) {
    options[option_name] = "1";
    //std::cout << "1" << std::endl;
  }

  map<string,string>::iterator ito;

  double cutoff = atof(options["--limit"].c_str()) - parse_time;
#ifdef _FLATZINC_OUTPUT
	cout << "%";
#endif
	std::cout << " d CUTOFF " << cutoff << std::endl;

  fm->set_strategy(options["--var_heuristic"], options["--val_heuristic"], options["--restart"]);
  s.initialise_random_seed(atoi(options["--seed"].c_str()));
  s.set_time_limit(cutoff);
  s.parameters.verbosity = atoi(options["--verbose"].c_str());
  fm->set_rewriting(atoi(options["--rewrite"].c_str()));

  fm->run(cout , p);
  
  //if (fm->finished())
  fm->print_final(cout , p);


#ifdef _VERIFICATION
  write_solution(fm, args.back());
#endif

  delete fm;
  return 0;
}
