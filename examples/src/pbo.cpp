

//#include <cmath>

#include <mistral_sat.hpp>
#include <stdlib.h>
#include <fstream>
 
using namespace std;



int main(int argc, char **argv)
{

  // double x = 2;
  // double k = -7;

  // double exp = pow(x,k);

  // std::cout << exp << std::endl;

  // exit(1);

  try 
    {
      // Define the command line object.
      SolverCmdLine cmd("Mistral (PBO)", ' ', "2.0");      
      cmd.parse(argc, argv);

  
  usrand(cmd.get_seed());
  
  Solver solver;
  
  cmd.set_parameters(solver);
  
  solver.parse_pbo(cmd.get_filename());
  
  solver.consolidate();
  
  
  if(cmd.print_model())
    std::cout << solver << std::endl;
  
  
  BranchingHeuristic *strategy;
  

  if(solver.parameters.backjump)
    strategy = new GenericHeuristic< VSIDS, RandomMinMax >(&solver);
      
      //GenericDVO< MaxWeight, 2, LiteralActivityManager >, MinValue >(&solver);
  else 
    strategy = new GenericHeuristic< WDEG, RandomMinMax >(&solver);
  
  
  RestartPolicy *policy = new Geometric(solver.parameters.restart_base,
					solver.parameters.restart_factor);
  
  
  Outcome result = solver.depth_first_search(solver.variables, strategy, policy);
  
  if(cmd.print_statistics())
    cout << solver.statistics ;
  
  if(result == SAT || result == OPT) {
    Solution sol(solver.variables);
    
    if(cmd.print_solution())
      cout << " c  " << sol << endl;
  }


    } catch (TCLAP::ArgException &e)  // catch any exceptions
    { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }
  
  return 0;
}



