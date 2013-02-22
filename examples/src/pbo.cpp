

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
  
  
   
  BranchingHeuristic *strategy;
  


  if(solver.parameters.backjump) {
    if(cmd.get_value_ordering() == "min")
      strategy = new GenericHeuristic< VSIDS<2>, MinValue >(&solver);
    else if(cmd.get_value_ordering() == "max")
      strategy = new GenericHeuristic< VSIDS<2>, MaxValue >(&solver);
    else if(cmd.get_value_ordering() == "rand")
      strategy = new GenericHeuristic< VSIDS<2>, RandomMinMax >(&solver);
    else if(cmd.get_value_ordering() == "weight")
      strategy = new GenericHeuristic< VSIDS<2>, BoolMinWeightValue >(&solver);
    else if(cmd.get_value_ordering() == "guided+min")
      strategy = new GenericHeuristic< VSIDS<2>, Guided< MinValue > >(&solver);
    else if(cmd.get_value_ordering() == "guided+max")
      strategy = new GenericHeuristic< VSIDS<2>, Guided< MaxValue > >(&solver);
    else  if(cmd.get_value_ordering() == "guided+rand")
      strategy = new GenericHeuristic< VSIDS<2>, Guided< RandomMinMax > >(&solver);
    else if(cmd.get_value_ordering() == "guided+weight")
      strategy = new GenericHeuristic< VSIDS<2>, Guided< BoolMinWeightValue > >(&solver);
    else {
      std::cout << "This value ordering is not handled!\n";
      exit(1);
    }
  } else { 
    strategy = new GenericHeuristic< WDEG, RandomMinMax >(&solver);
  }


  if(cmd.print_model())
    std::cout << solver << std::endl;

  

  
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



