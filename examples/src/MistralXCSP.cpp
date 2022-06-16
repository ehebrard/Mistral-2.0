#include <signal.h> 

#include <XCSP3CoreParser.h>
#include "XCSP3MistralCallbacks.hpp"

#include <mistral_search.hpp>

using namespace XCSP3Core;


XCSP3MistralCallbacks * cb_ptr;

void parse(XCSP3MistralCallbacks& cb, const char* instancefile) {
  // try
  // {
  XCSP3CoreParser parser(&cb);

  // cout << "parse\n";

  parser.parse(instancefile); // fileName is a string
  // }
  // catch (exception &e)
  // {
  // 	cout.flush();
  // 	cout << "\nc \tUnexpected exception :\n c \t" << e.what() << endl;
  // 	cout << "s UNSUPPORTED" << endl;
  // 	exit(1);
  // }

  // cout << "consolidate\n";

  cb.solver.consolidate();

  // cout << "set goal\n";

  cb.solver.set_goal(cb.goal);

  // cout << "end parse\n";
}


void print_outcome(XCSP3MistralCallbacks& cb, std::ostream& os) {
	Outcome result = cb.solver.statistics.outcome;
	
  switch(result) {
  case SAT: 
    os << "s SATISFIABLE\n" ;
    break;
  case OPT: 
    os << "s OPTIMUM FOUND\n" ;
    break;
  case UNSAT: 
    os << "s UNSATISFIABLE\n" ;
    break;
  case UNKNOWN: 
    os << "s UNKNOWN\n" ;
    break;
  case LIMITOUT: 
    if(cb.solver.statistics.num_solutions > 0) 
      os << "s SATISFIABLE\n" ;
    else 
      os << "s UNKNOWN\n" ;
  }
}


void print_solution(XCSP3MistralCallbacks& cb, std::ostream& os, char='v') {
	if(cb.solver.statistics.num_solutions > 0) {
		os << "v <instantiation type=\"";
		if(cb.solver.statistics.outcome == OPT)
			os << "optimum\" cost=\"" << cb.solver.statistics.objective_value << "\">\n";
		else
			os << "solution\">\n"; 
		os << "v   <list>";
		for( auto id : cb.var_ids ) {
			os << " " << id;
		}
		os << " </list>\nv   <values>";
		
		
		for( size_t i=0; i<cb.variables.size; ++i ) {
			Variable var = cb.variables[i];
			int deg = cb.initial_degree[i];


			// std::cerr << var << " " << cb.declared_var_ids[i] << " " << deg << std::endl;

			if(deg==0 && var.get_initial_min() < var.get_initial_max()) {
				os << " *";
			} else if(var.id()>=0)
				os << " " << cb.solver.last_solution_lb[var.id()];
			else
				os << " " << var.get_value();
		}
		
		os << " </values>\nv </instantiation>\n";
		
	}
}


class ObjectivePrinter : public SolutionListener {

public:
	Solver *solver;
	
	ObjectivePrinter(Solver *s) : SolutionListener(), solver(s) {}
	
	void notify_solution() {
		std::cout << "o " << solver->objective->value() << std::endl;
	}

};


static void Mistral_SIGTERM_handler(int signum) {

	if(cb_ptr->solver.statistics.num_solutions > 0) {
		std::cout << "s SATISFIABLE\n" ;
		print_solution(*cb_ptr, std::cout);
	}
	else
		std::cout << "s UNKNOWN\n" ;
	exit(0);
}




int main(int argc,char **argv) {
	
	std::cout << "c Mistral XCSP3" << std::endl;
	SolverCmdLine cmd("Mistral (xcsp3)", ' ', "2.0"); 
	
  TCLAP::SwitchArg simple_rewriteArg("","simple_rewrite","Uses simple rewriting", false);
  cmd.add( simple_rewriteArg );

  TCLAP::ValueArg<int> branchOnaux("","branch_on_aux","Branching on auxiliary variables whose domain size is larger than value; if it is given -1, the solver will branch on all variables", false, 0, "int");
  cmd.add( branchOnaux );

  TCLAP::ValueArg<int> countArg("","count","count solutions (up to value)", false, 1, "int");
  cmd.add( countArg );


  TCLAP::SwitchArg minimum_outputArg("","minimum_output","use this option to avoid the output used for the competition", false);
  cmd.add( minimum_outputArg );


  TCLAP::ValueArg<int>  recommended_Arg("","recommended","use our recommended configurations for search", false, 1, "int");
  cmd.add( recommended_Arg );
	
	cmd.parse(argc, argv);
	
	
	usrand(cmd.get_seed());
	
	
	Solver solver;
	cmd.set_parameters(solver);
	
	
	// double cpu_time = get_run_time() ;
	
	XCSP3MistralCallbacks cb(solver); // my interface between the parser and the solver
	parse(cb, cmd.get_filename().c_str());
	
	cb_ptr= & cb;
	signal(SIGTERM,Mistral_SIGTERM_handler);


	//std::cout << "c parsetime " << (get_run_time() - cpu_time) << std::endl;
	
	if(cmd.print_model())
		std::cout << solver << std::endl;
		
	// BranchingHeuristic *heuristic = solver.heuristic_factory(cmd.get_variable_ordering(), cmd.get_value_ordering(), cmd.get_randomization());
	// RestartPolicy *restart = solver.restart_factory(cmd.get_restart_policy());
	
	
	if(solver.parameters.time_limit > 0) {
		solver.parameters.time_limit -= get_run_time();
		if(solver.parameters.time_limit < 0)
			solver.set_time_limit(0.5);
	}

	if(solver.objective && solver.objective->is_optimization())
		if (!minimum_outputArg.getValue())
			solver.add( new ObjectivePrinter(&solver) );

	
	BranchingHeuristic *heuristic;
	RestartPolicy *restart;


	if (!recommended_Arg.getValue()){
		heuristic = solver.heuristic_factory(cmd.get_variable_ordering(), cmd.get_value_ordering(), cmd.get_randomization());
		restart = solver.restart_factory(cmd.get_restart_policy());
	}
	else if (recommended_Arg.getValue()==1){
		//std::cout << " c search recommendation nb 1" << std::endl;
		if(solver.objective && solver.objective->is_optimization()) {
			restart	= new Luby();
			restart->base = 128;
			heuristic = new LastConflict < GenericDVO < MinDomainOverWeight, 1, ConflictCountManager >, SolutionGuided< MinValue, RandomMinMax >, SolutionGuided< MinValue, RandomMinMax >, 1 > (&solver);
		} else {
			restart = new Geometric();
			heuristic = new LastConflict < GenericDVO < MinDomainOverWeight, 1, ConflictCountManager >, SolutionGuided< MinValue, RandomMinMax >, SolutionGuided< MinValue, RandomMinMax >, 1 > (&solver);
		}
	}
	else if (recommended_Arg.getValue()==2){
		//std::cout << " c search recommendation nb 1" << std::endl;
		if(solver.objective && solver.objective->is_optimization()) {
			//restart	= new Luby();
			restart = new Geometric();
			//restart->base = 128;
			heuristic = new LastConflict < GenericDVO < MinDomainOverWeight, 1, ConflictCountManager >, Guided< MinValue >, Guided< MinValue >, 1 > (&solver);
		} else {
			//restart = new Luby();
			restart = new Geometric();
			heuristic = new LastConflict < GenericDVO < MinDomainOverWeight, 2, ConflictCountManager >,  Guided< MinValue >,  Guided< MinValue >, 1 > (&solver);
		}
	}
	else {
		std::cout << "c search recommendation not found" << std::endl;
		exit(1);
		}

	
	// std::cout << solver.constraints[277].binary() << std::endl;
	
	if(countArg.getValue()!=1) {
		
		int num_solutions = 0;
		// Outcome res = solver.depth_first_search(cb.variables, heuristic, new NoRestart());
		
		
		VarArray sequence;
		for( auto x : cb.variables ) {
			if(x.get_degree()>0) {
				sequence.add(x);
			}
		}
		
		solver.initialise_search(sequence, heuristic, new NoRestart());
		
		while( (countArg.getValue()<1 || countArg.getValue()>num_solutions) && solver.get_next_solution() == SAT ) {
			++num_solutions;
			
			if(num_solutions <= 1000 || (num_solutions <= 10000 && (num_solutions%100 == 0)) || (num_solutions%1000 == 0)) {
				cout << " numsol = " << num_solutions ;
			
				for( auto var : cb.variables ) {
					if(var.id()>=0)
						cout << " " << cb.solver.last_solution_lb[var.id()];
					else
						cout << " " << var.get_value();
				}
				cout << endl;
			
				ofstream solfile("sols/sol"+int2str(num_solutions)+".txt", ofstream::out);
				solfile << "s SATISFIABLE\n";
				print_solution(cb, solfile);
				solfile.close();
			}	
			
			
		} 
		cout << num_solutions << endl;
		
		
	} else {
	
		if(branchOnaux.getValue()>0)
		{
			Vector<Variable> search_sequence;
			BitSet search_vars(0, solver.variables.size-1, BitSet::empt);

			for(unsigned int k=0; k<cb.variables.size; ++k) {
				search_vars.add(cb.variables[k].id());
				search_sequence.add(cb.variables[k]);
			}

			for(unsigned int i=0; i<solver.variables.size; ++i) {
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
			/*VarArray sequence;
			for( auto x : cb.variables ) {
				if(x.get_degree()>0) {
					sequence.add(x);
				}
			}
			*/
			if (!branchOnaux.getValue())
				solver.depth_first_search(cb.variables, heuristic, restart);
			else
				solver.depth_first_search(solver.variables, heuristic, restart);
	  }	

		if(cmd.print_statistics())
			solver.statistics.print_full(std::cout);
		else
			print_outcome(cb, std::cout);

		if (!minimum_outputArg.getValue())
			print_solution(cb, std::cout);

	}

  return 0;
}

