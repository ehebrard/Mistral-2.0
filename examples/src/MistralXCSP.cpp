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
	
	// cb.solver.set_goal(cb.goal);
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
	
	
	// cout << solver << endl;
	
	BranchingHeuristic *heuristic = solver.heuristic_factory(cmd.get_variable_ordering(), cmd.get_value_ordering(), cmd.get_randomization());
	RestartPolicy *restart = solver.restart_factory(cmd.get_restart_policy());
	
	Outcome result = UNKNOWN;
	if(branchOnaux.getValue())
		result = solver.depth_first_search(solver.variables, heuristic, restart, cb.goal);
	else
		result = solver.depth_first_search(cb.variables, heuristic, restart, cb.goal);
	
	cout << result << endl;

  return 0;
}

