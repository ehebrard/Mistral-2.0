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


#ifdef _PARALLEL
#include "omp.h"
#include <limits>
#include <math.h>
#endif

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
  
  SolverCmdLine cmd("Mistral (fzn)", ' ', "2.0");      
  
  TCLAP::SwitchArg annotationArg("","follow_annotations","Uses the annotations", false);
  cmd.add( annotationArg );

  TCLAP::ValueArg<int> parityArg("","parity","Uses parity processing", false, 0, "int");
  cmd.add( parityArg );

  TCLAP::SwitchArg simple_rewriteArg("","simple_rewrite","Uses simple rewriting", false);
  cmd.add( simple_rewriteArg );

  TCLAP::SwitchArg branchOnaux("","branch_on_aux","Branching on auxiliary variables", false);
  cmd.add( branchOnaux );

  TCLAP::ValueArg<int>  recommended_Arg("","recommended","use our recommended configurations for search", false, 0, "int");
  cmd.add( recommended_Arg );


  TCLAP::ValueArg<std::string> output_m("","output-mode","output form ",false,"c","string");
  cmd.add( output_m );

  TCLAP::SwitchArg output_o("","output-objective","output objective ",false);
  cmd.add( output_o );


//  pcommentArg = new TCLAP::ValueArg<std::string>("","prefix_comment","output comments prefix",false,"c","string");


#ifdef _PARALLEL
  //std::cout << "PARALLEL \n \n \n " << std::endl;
  TCLAP::ValueArg<int> threadsArg("p","number-threads","Use multithreading with this option.", false, 4, "int");
  cmd.add( threadsArg );
#endif

  cmd.parse(argc, argv);

#ifdef _PARALLEL
  int total = threadsArg.getValue();
  if (cmd.enumerate_solutions())
	  total=1;
  else{
	  //std::cout << "Available  threads : " << omp_get_max_threads() << std::endl;
	  //int recommended= floor((double) (omp_get_max_threads()*3) / 4.0 )   ;
	  int recommended= floor((double) (omp_get_max_threads()) / 2.0 )   ;
	  if (total >recommended){
		  //std::cout << " % " << " high value of -p. The solver will use only " << recommended << " threads (recommended) " << std::endl;
		  total=recommended;
	  }
  }
//  else
//	  std::cout << " % " << " will use " << total << " threads " << std::endl;
  omp_set_num_threads(total);
  long int global_obj =std::numeric_limits<int>::max();
  bool solution_found_elsewhere= false;

#endif



  #ifdef _PARALLEL
#pragma omp parallel
  {
	  //printf("multicore user! %d \n" , omp_get_thread_num());
	  int id = omp_get_thread_num();
#endif


  int thread_seed = cmd.get_seed();

#ifdef _PARALLEL
  thread_seed+=id;
#endif
  usrand(thread_seed);

  Solver s;

  cmd.set_parameters(s);

  std::string policy;
#ifdef _PARALLEL
  if ( id%2 != 0)
	  if (cmd.get_restart_policy()=="geom")
		  policy="luby";
	  else
		  policy="geom";
  else
#endif
	  policy =cmd.get_restart_policy();

  bool branch_on_auxilary =   branchOnaux.getValue();
  //std::cout << "branch_on_auxilary " << branch_on_auxilary << std::endl;
  //exit(1);


#ifdef _PARALLEL
//  bool branch_on_auxilary = true;
  if ( id%4 < 2 )
	  branch_on_auxilary=false;
//#pragma omp critical
//  std::cout << " " << s.parameters.prefix_statistics << "ID:  " << id  << " " << policy << " branch_on_auxilary " << branch_on_auxilary <<  "  "<< thread_seed <<  std::endl;
#endif

  double cpu_time = get_run_time() ;

#ifdef _VERBOSE_PARSER
  std::cout << " " << s.parameters.prefix_comment << " Parse: ";
#endif

  FlatZinc::Printer p;
  FlatZinc::FlatZincModel *fm = 0L;

  fm = parse(cmd.get_filename(), s, p);


  if( !fm )
#ifdef _PARALLEL
	  exit(1);
#else
  return 0;
#endif

  fm->set_enumeration(cmd.enumerate_solutions());
  double parse_time = get_run_time() - cpu_time;
#ifdef _VERBOSE_PARSER
  std::cout << std::endl;
#endif

#ifdef _PARALLEL
 //#pragma omp critical
  //cout << " " << s.parameters.prefix_statistics << " PARSETIME " << parse_time << std::endl;
#else
  if(s.parameters.verbosity >0)
	  cout << " " << s.parameters.prefix_statistics << " PARSETIME " << parse_time << std::endl;
#endif

  FlatZinc::SolutionPrinter *sp = new FlatZinc::SolutionPrinter(&p, fm, &s);
  s.add(sp);


  if(s.parameters.time_limit>0) std::cout << " " << s.parameters.prefix_statistics 
					  << " CUTOFF " << s.parameters.time_limit << std::endl;


  fm->branch_on_auxilary=branch_on_auxilary;

  if (!recommended_Arg.getValue())
	  fm->set_strategy(cmd.get_variable_ordering(), cmd.get_value_ordering(), cmd.get_randomization(), policy);
  else
	  fm->recommended_configs(recommended_Arg.getValue());

//  fm->set_strategy(cmd.get_variable_ordering(), cmd.get_value_ordering(), cmd.get_randomization(), policy);
  fm->set_display_model(cmd.print_model());
  fm->set_display_solution(cmd.print_solution());
  fm->set_annotations(annotationArg.getValue());
  //fm->set_annotations(false);
  fm->set_rewriting(cmd.use_rewrite());
  fm->set_simple_rewriting(simple_rewriteArg.getValue());
  fm->set_parity_processing(parityArg.getValue());
  fm->encode_clauses();

#ifdef _PARALLEL
  if (fm->method() == FlatZinc::FlatZincModel::MAXIMIZATION)
	  global_obj = std::numeric_limits<int>::min();
  fm->best_kown_objective = &global_obj;
  fm->set_solution_found_elsewhere(&solution_found_elsewhere);
#endif

  fm->run(cout , p);


#ifdef _PARALLEL
#pragma omp critical
  {
  if (!solution_found_elsewhere){
	  solution_found_elsewhere=true;
	  //Secure shared memory
#pragma omp flush
 // }
 // else
#endif

  if(cmd.print_solution())
    fm->print_final(cout , p);


#ifdef _PARALLEL
}
  }
#endif

  if(cmd.print_statistics())
    s.statistics.print_full(std::cout);



#ifdef _VERIFICATION
  write_solution(fm, args.back());
#endif


  delete fm;
  delete sp;
  //exit(1);

#ifdef _PARALLEL
  }
#endif

  return 0;
}
