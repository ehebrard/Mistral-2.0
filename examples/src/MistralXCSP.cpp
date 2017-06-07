#include <XCSP3CoreParser.h>
#include "XCSP3MistralCallbacks.hpp"

#include <mistral_search.hpp>

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
		cout << "\n c \tUnexpected exception :\n c \t" << e.what() << endl;
		cout << " s UNSUPPORTED" << endl;
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


void print_solution(XCSP3MistralCallbacks& cb, std::ostream& os, char='v') {
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
		
		for( auto var : cb.variables ) {
			if(var.id()>=0)
				os << " " << cb.solver.last_solution_lb[var.id()];
			else
				os << " " << var.get_value();
		}
		
		// int linecount = 0;
		// for( auto var : cb.variables ) {
		// 	if(linecount%8 == 0) os << endl;
		// 	++linecount;
		//
		// 	if(var.id()>=0)
		// 		os << " " << cb.solver.last_solution_lb[var.id()];
		// 	else
		// 		os << " " << var.get_value();
		// }
		// os << endl;
		
		os << " </values>\n v </instantiation>\n";
		
	}
}


class ObjectivePrinter : public SolutionListener {

public:
	Solver *solver;
	
	ObjectivePrinter(Solver *s) : SolutionListener(), solver(s) {}
	
	void notify_solution() {
		std::cout << " o " << solver->objective->value() << std::endl;
	}

};



int main(int argc,char **argv) {
	
	SolverCmdLine cmd("Mistral (xcsp3)", ' ', "2.0"); 
	
  TCLAP::SwitchArg simple_rewriteArg("","simple_rewrite","Uses simple rewriting", false);
  cmd.add( simple_rewriteArg );

  TCLAP::ValueArg<int> branchOnaux("","branch_on_aux","Branching on auxiliary variables whose domain size is larger than value", false, 0, "int");
  cmd.add( branchOnaux );

  TCLAP::ValueArg<int> countArg("","count","count solutions (up to value)", false, 1, "int");
  cmd.add( countArg );

	
	cmd.parse(argc, argv);
	
	
	usrand(cmd.get_seed());
	
	
	Solver solver;
	cmd.set_parameters(solver);
	
	
	double cpu_time = get_run_time() ;
	
	XCSP3MistralCallbacks cb(solver); // my interface between the parser and the solver
	parse(cb, cmd.get_filename().c_str());
	
	
	std::cout << " c parsetime " << (get_run_time() - cpu_time) << std::endl;
	
	if(cmd.print_model())
		std::cout << solver << std::endl;
		
	BranchingHeuristic *heuristic = solver.heuristic_factory(cmd.get_variable_ordering(), cmd.get_value_ordering(), cmd.get_randomization());
	RestartPolicy *restart = solver.restart_factory(cmd.get_restart_policy());
	
	
	if(solver.parameters.time_limit > 0) {
		solver.parameters.time_limit -= get_run_time();
		if(solver.parameters.time_limit < 0)
			solver.set_time_limit(0.5);
	}
	
	
	if(solver.objective && solver.objective->is_optimization()) {
		solver.add( new ObjectivePrinter(&solver) );
	}
	
	// std::cout << solver.constraints[277].binary() << std::endl;
	
	if(countArg.getValue()!=1) {
		
		int num_solutions = 0;
		// Outcome res = solver.depth_first_search(cb.variables, heuristic, new NoRestart());
		solver.initialise_search(cb.variables, heuristic, new NoRestart());
		
		while( (countArg.getValue()<1 || countArg.getValue()>num_solutions) && solver.get_next_solution() == SAT ) {
			++num_solutions;
			
			cout << " numsol = " << num_solutions ;
			
			for( auto var : cb.variables ) {
				if(var.id()>=0)
					cout << " " << cb.solver.last_solution_lb[var.id()];
				else
					cout << " " << var.get_value();
			}
			cout << endl;
			
			ofstream solfile("sols/sol"+int2str(num_solutions)+".txt", ofstream::out);
			solfile << " s SATISFIABLE\n";
			print_solution(cb, solfile);
			solfile.close();
			
			
		} 
		cout << num_solutions << endl;
		
		
	} else {
	
		if(branchOnaux.getValue()>0)
		{
			Vector<Variable> search_sequence;
			BitSet search_vars(0, solver.variables.size-1, BitSet::empt);

			for(int k=0; k<cb.variables.size; ++k) {
				search_vars.add(cb.variables[k].id());
				search_sequence.add(cb.variables[k]);
			}

			for(int i=0; i<solver.variables.size; ++i) {
				int domsize = solver.variables[i].get_size();

				if(//solver.variables[i].is_boolean()
						domsize>1 && domsize<=branchOnaux.getValue()
						&& !(search_vars.contain(solver.variables[i].id()))) {
					search_vars.add(i);
					search_sequence.add(solver.variables[i]);
				}
			}
			solver.depth_first_search(search_sequence, heuristic, restart);
		}
		else {
			// if (annotationArg.getValue())
				solver.depth_first_search(cb.variables, heuristic, restart);
			// else
			// 	solver.depth_first_search(solver.variables, heuristic, restart);
	  }	
	}

	if(cmd.print_statistics())
		solver.statistics.print_full(std::cout);
	else
		print_outcome(cb, std::cout);

	print_solution(cb, std::cout);

  return 0;
}

