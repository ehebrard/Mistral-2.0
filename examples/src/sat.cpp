
#include <mistral_sat.hpp>
#include <stdlib.h>
#include <fstream>
 
using namespace std;



const int nia = 9;
const char* int_ident[nia] = {
  "-h",
  "-verbose",
  "-time_limit",
  "-seed",
  "-randomize", 
  "-restart_base",
  "-checked",
  "-activity",
  "-value"
};

const char* int_man[nia] = {
  "help message",
  "verbosity     {0,1}", 
  "time limit",
  "random seed",
  "0 -> none, >0 -> shuffle vars on restart, k -> randomized heuristic", 
  "restart base",
  "whether the solution is checked",
  "initialisation method for the activity (0: none, 1: number of clauses)",
  "value selection method (0,1,2/3:most/least active,4:random,5:persistent)"
   };

			      
int int_param[nia];

const int nsa = 5;
const char* str_ident[nsa] = {
  "-restart_factor",
  "-decay",
  "-forget",
  "-policy",
  "-normalize"
};
const char* str_man[nsa] = {
  "restart factor",
  "1 - activity decay rate",
  "forgetfulness ratio",
  "restart policy {luby, geom}",
  "activity is normalized before search"
};
const char* str_param[nsa];


void outputHelpMessage()
{
  cerr << "\nUsage: \t bin/sat [path to the cnf file] <args>" << endl << endl;
  for(int i=0; i<nia; ++i)
    cerr  << setw(20) << int_ident[i] << ": \t" << int_man[i] << endl;
  for(int i=0; i<nsa; ++i)
    cerr  << setw(20) << str_ident[i] << ": \t" << str_man[i] << endl;
  cerr << endl;
}


int main(int argc, char **argv)
{

  str_param[0] = "";
  getCommandLine(int_ident, int_param, nia,
		 str_ident, str_param, nsa,
		 &(argv[1]), argc-1);
  
  if( argc < 2 || int_param[0] != NOVAL || !strcmp( argv[argc-1], "-h" ) )
    {
      outputHelpMessage();
    }
  else
    {
      SolverParameters params;

      params.verbosity       = ( int_param[1] != NOVAL ? int_param[1] : 4 );
      params.time_limit      = ( (double)(int_param[2] != NOVAL ? int_param[2] : -1 ) );
      params.seed            = ( int_param[3] != NOVAL ? int_param[3] : 11041979 );	  
      params.randomization   = ( int_param[4] != NOVAL ? abs(int_param[4]) : 2 );
      if(params.randomization == 0)
	params.randomization = 1;
      params.shuffle         = ( int_param[4] > 0 );
      params.restart_base    = ( int_param[5] != NOVAL ? int_param[5] : 200 );
      params.restart_limit   = ( int_param[5] != NOVAL ? int_param[5] : 200 );
      params.checked         = ( int_param[6] != NOVAL ? int_param[6] : 1 );
      params.init_activity   = ( int_param[7] != NOVAL ? int_param[7] : 1 );
      params.value_selection = ( int_param[8] != NOVAL ? int_param[8] : 2 );


      params.restart_factor     = ( strcmp(str_param[0],"nil") ? atof(str_param[0]) : 1.05 );
      params.activity_decay     = ( strcmp(str_param[1],"nil") ? atof(str_param[1]) : .96 );
      params.forgetfulness      = ( strcmp(str_param[2],"nil") ? atof(str_param[2]) : .75 );
      params.restart_policy     = ( ( strcmp(str_param[3],"luby") ? ( strcmp(str_param[3],"no") ? NORESTART : GEOMETRIC) : LUBY ) );
      params.normalize_activity = ( strcmp(str_param[4],"nil") ? atof(str_param[4]) : 0 );
      params.dynamic_value      = (params.value_selection>5);

      
      //SatSolver solver;

      Solver solver;

      solver.parse_dimacs(argv[1]);

      //cout << solver << endl;

      //solver.propagate();

      solver.solve();
    }
  return 0;
}



