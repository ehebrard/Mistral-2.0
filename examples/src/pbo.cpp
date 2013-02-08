
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

const int nsa = 6;
const char* str_ident[nsa] = {
  "-factor",
  "-decay",
  "-forget",
  "-policy",
  "-normalize",
  "-algo"
};
const char* str_man[nsa] = {
  "restart factor",
  "1 - activity decay rate",
  "forgetfulness ratio",
  "restart policy {luby, geom}",
  "activity is normalized before search",
  "method: pbo / cp"
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
  get_command_line(int_ident, int_param, nia,
		 str_ident, str_param, nsa,
		 &(argv[1]), argc-1);
  
  if( argc < 2 || int_param[0] != NOVAL || !strcmp( argv[argc-1], "-h" ) )
    {
      outputHelpMessage();
    }
  else
    {
      SolverParameters params;

      params.verbosity       = ( int_param[1] != NOVAL ? int_param[1] : 2 );
      params.time_limit      = ( (double)(int_param[2] != NOVAL ? int_param[2] : -1 ) );
      params.seed            = ( int_param[3] != NOVAL ? int_param[3] : 11041979 );	  
      params.randomization   = ( int_param[4] != NOVAL ? abs(int_param[4]) : 1 );
      if(params.randomization == 0)
	params.randomization = 1;
      params.shuffle         = ( int_param[4] != NOVAL ? int_param[4] : 0 );
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

      params.activity_increment = 0.012;

      int method = strcmp(str_param[5],"cp");
     
      Outcome result = UNKNOWN;

      Solver solver;
      solver.parameters.copy(params);
      solver.parse_pbo(argv[1]);
      solver.parameters.verbosity = 2;



       
      Variable *X = solver.variables.stack_;

      //solver.add(X[3] == 1);
      // solver.add(X[20] == 1);
      // solver.add(X[35] == 1);
      // solver.add(X[50] == 1);
      // solver.add(X[58] == 1);
      // solver.add(X[67] == 1);
      // solver.add(X[78] == 1);
      // solver.add(X[96] == 1);
      // solver.add(X[105] == 1);
      // solver.add(X[125] == 1);
      // solver.add(X[140] == 1);
      // solver.add(X[155] == 1);
      // solver.add(X[160] == 1);


      solver.consolidate();





      std::cout << solver << std::endl;
      
      if(method) {

	solver.parameters.backjump = 1;
	result = solver.depth_first_search(solver.variables, 
					   //new VSIDS(&solver),
					   new GenericHeuristic<
					   GenericDVO< MaxWeight, 1, LiteralActivityManager >, MinValue >(&solver),
					   new Geometric(params.restart_base,params.restart_factor)
					   );

      } else {

	result = solver.depth_first_search(solver.variables, 
					   new GenericHeuristic< WDEG, RandomMinMax >(&solver),
					   new Geometric(params.restart_base,params.restart_factor)
					   );
	
      }


      if(result == SAT) {
	Solution sol(solver.variables);
	cout << solver.statistics << endl << sol << endl;
      }
    }
  return 0;
}



