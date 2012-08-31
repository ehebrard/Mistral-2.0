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
  options["--limit"] = "0";
  options["--verbose"] = "0";
  options["--rewrite"] = "0";
  options["--display_model"] = "0";


  options["a"] = "0";
  options["p"] = "1";
  options["f"] = "0";
  



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
	cout << "X";
  fm = parse(args.back(), s, p);

	cout << "Y";
  FlatZinc::SolutionPrinter *sp = new FlatZinc::SolutionPrinter(&p, fm, &s);
  s.add(sp);

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

  double total_time = atof(options["--limit"].c_str());
  double cutoff = (total_time>0 ? total_time - parse_time : 0);
#ifdef _FLATZINC_OUTPUT
  cout << "%";
#endif
  if(cutoff>0) std::cout << " d CUTOFF " << cutoff << std::endl;


  // set flatzinc model options
  fm->set_strategy(options["--var_heuristic"], options["--val_heuristic"], options["--restart"]);
  fm->set_rewriting(atoi(options["--rewrite"].c_str()));
  fm->set_enumeration(atoi(options["-a"].c_str()));
  fm->set_display_model(atoi(options["--display_model"].c_str()));
  fm->set_annotations(!(atoi(options["-f"].c_str())));

  // set solver options
  s.initialise_random_seed(atoi(options["--seed"].c_str()));
  s.parameters.verbosity = atoi(options["--verbose"].c_str());
  s.set_time_limit(cutoff);



  fm->run(cout , p);
  
  //if (fm->finished())
  fm->print_final(cout , p);


#ifdef _VERIFICATION
  write_solution(fm, args.back());
#endif

  delete fm;
  delete sp;
  return 0;
}
