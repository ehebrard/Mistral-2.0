#include <XCSP3CoreParser.h>
#include "XCSP3MistralCallbacks.hpp"


using namespace XCSP3Core;


void parse(XCSP3MistralCallbacks& cb, const char* instancefile) {

	try
	{
		XCSP3CoreParser parser(&cb);
		parser.parse(instancefile); // fileName is a string
	}
	catch (exception &e)
	{
		cout.flush();
		cerr << "\n\tUnexpected exception :\n";
		cerr << "\t" << e.what() << endl;
		exit(1);
	}
	
	cb.solver.consolidate();
	
	cb.solver.set_goal(cb.goal);
}


void print_outcome(XCSP3MistralCallbacks& cb, std::ostream& os) {
	Outcome result = cb.solver.statistics.outcome;
	
  switch(result) {
  case SAT: 
    os << " s SATISFIABLE\n" ;
    break;
  case OPT: 
    os << " s OPTIMUM FOUND\n" ;
    break;
  case UNSAT: 
    os << " s UNSATISFIABLE\n" ;
    break;
  case UNKNOWN: 
    os << " s UNKNOWN\n" ;
    break;
  case LIMITOUT: 
    if(cb.solver.statistics.num_solutions > 0) 
      os << " s SATISFIABLE\n" ;
    else 
      os << " s UNKNOWN\n" ;
  }
}


void print_solution(XCSP3MistralCallbacks& cb, std::ostream& os) {
	if(cb.solver.statistics.num_solutions > 0) {
		os << " v <instantiation type=\"";
		if(cb.solver.statistics.outcome == OPT)
			os << "optimum\" cost=\"" << cb.solver.statistics.objective_value << "\">\n";
		else
			os << "solution\">\n"; 
		os << " v   <list>";
		for( auto id : cb.var_ids ) {
			os << " " << id;
		}
		os << " </list>\n v   <values>";
		
		// int linecount = 0;
		for( auto var : cb.variables ) {
			// if(linecount%50 == 0) os << endl;
			// ++linecount;
			
			os << " " << cb.solver.last_solution_lb[var.id()];
		}
		// os << endl;
		
		os << " </values>\n v </instantiation>\n";
	}
}



int main(int argc,char **argv) {
	
	SolverCmdLine cmd("Mistral (xcsp3)", ' ', "2.0"); 
	
  TCLAP::SwitchArg simple_rewriteArg("","simple_rewrite","Uses simple rewriting", false);
  cmd.add( simple_rewriteArg );

  TCLAP::SwitchArg branchOnaux("","branch_on_aux","Branching on auxiliary variables", false);
  cmd.add( branchOnaux );

	
	cmd.parse(argc, argv);
	
	
	usrand(cmd.get_seed());
	
	
	Solver solver;
	cmd.set_parameters(solver);
	
	
	double cpu_time = get_run_time() ;
	
	XCSP3MistralCallbacks cb(solver); // my interface between the parser and the solver
	parse(cb, cmd.get_filename().c_str());
	
	
	std::cout << " c parsetime " << (get_run_time() - cpu_time) << std::endl;
	
	// std::cout << solver << std::endl;
		
	BranchingHeuristic *heuristic = solver.heuristic_factory(cmd.get_variable_ordering(), cmd.get_value_ordering(), cmd.get_randomization());
	RestartPolicy *restart = solver.restart_factory(cmd.get_restart_policy());
	
	if(branchOnaux.getValue())
		solver.depth_first_search(solver.variables, heuristic, restart);
	else
		solver.depth_first_search(cb.variables, heuristic, restart);
	
	print_outcome(cb, std::cout);
	print_solution(cb, std::cout);

  return 0;
}

