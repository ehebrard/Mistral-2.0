 

#include "mistral_scheduler.hpp"


using namespace Mistral;




class DSMI {

public:
	
	SchedulingSolver *solver;
	
	std::vector<int> values_lb;
	std::vector<int> levels;

	DSMI(SchedulingSolver *s) : solver(s) {}


	
	void evaluate() {
		// std::cout << "\nhello incr backtrack at level " << solver->level << " objective = " << solver->get_objective_var().get_domain() << " \n";
		// values_lb.clear();
		// levels.clear();
		solver->get_objective_var().update_lb_history(values_lb, levels);
		
		
		

		auto v{values_lb.begin()};
		auto l{levels.begin()};
		while(v != values_lb.end()) {
			std::cout << *v << " at level " << *l << std::endl;
			++v;
			++l;
		}
		assert(l == levels.end());
	}

};

class EvalBacktrack : public BacktrackListener {

public:
	
	DSMI *criterion;

	EvalBacktrack(DSMI *c) : BacktrackListener(), criterion(c) {}
		
	void notify_backtrack() {
		std::cout << "Backtrack:\n";
		criterion->evaluate();
		std::cout << "lb >= " <<  criterion->solver->get_objective_var().get_max() << " at level " << criterion->solver->level << std::endl;
	}

};

class EvalSolution : public SolutionListener {

public:
	
	DSMI *criterion;

	EvalSolution(DSMI *c) : SolutionListener(), criterion(c) {}
		
	void notify_solution() {
		std::cout << "Solution:\n";
		criterion->evaluate();
	}

};




int main( int argc, char** argv )
{

  ParameterList params(argc, argv);
  usrand(params.Seed);

  StatisticList stats;
  stats.start();

  Instance jsp(params);
  
  std::cout << std::endl;
  
	jsp.print(std::cout);
  
	jsp.printStats(std::cout);
  params.print(std::cout);


  SchedulingSolver *solver;
  if(params.Objective == "makespan") {
    std::cout << "c Minimising Makespan" << std::endl;
    if(params.Type == "now") solver = new No_wait_Model(jsp, &params, -1, 0);
    else if(params.Type == "now2") {
      //params.Type = "now";
      solver = new No_wait_Model(jsp, &params, -1, 1);
    }
    else solver = new C_max_Model(&jsp, &params, &stats);
  } else if(params.Objective == "tardiness") {
    std::cout << "c Minimising Tardiness" << std::endl;
    solver = new L_sum_Model(&jsp, &params, &stats);
  } // else if(params.Objective == "depth") {
  //   std::cout << "c Minimising Depth" << std::endl;
  //   solver = new Depth_Model(jsp, &params, jsp.getMakespanUpperBound());
  // } else if(params.Objective == "weight") {
  //   std::cout << "c Minimising Weight" << std::endl;
  //   solver = new DTP_Model(jsp, &params, 1000);
  // }
  else {
    std::cout << "c unknown objective, exiting" << std::endl;
    exit(1);
  }
	
	solver->setup();
	

  // SchedulingSolver solver(model, &params, &stats);
  usrand(params.Seed);
  
  solver->consolidate();

  //std::cout << solver << std::endl;
  
	BranchingHeuristic *heu = new SchedulingWeightedDegree < TaskDomOverBoolWeight, Guided< MinValue >, 2 > (solver, solver->disjunct_map);
	


	/** Les variables de decision (les disjonctions) **/
	std::vector<Variable> X;
	for(auto i{0}; i<solver->disjuncts.size; ++i) {
		X.push_back(solver->disjuncts[i].get_var());
	}
	auto offset{solver->disjuncts[0].id()};
	
	// for(auto x : X) {
	// 	std::cout << x
	// 		// << " / " << solver->disjuncts[solver->sequence[i].id() - offset]
	// 		<< ": "
	// 		<< solver->disjunct_map[x.id() - offset][1] <<
	// 			" in " << solver->disjunct_map[x.id() - offset][1].get_domain() << " <> "
	// 			<< solver->disjunct_map[x.id() - offset][2] << " in "
	// 				<< solver->disjunct_map[x.id() - offset][2].get_domain() << "\n";
	// }
	
	
	
  solver->dichotomic_search(heu);

	
	

	/** Remove assigned variables **/
	X.erase(std::remove_if(X.begin(), 
                       X.end(),
                       [](const Variable& x){return x.is_ground();}),
        X.end());
	
	/** Remove assigned variables **/
	std::sort(X.begin(), X.end(), [&](const Variable& x, const Variable &y) { return heu->comp(x,y); });
	

	
	
	std::cout << std::endl << "SAVE " << solver->level << "|" << solver->decisions.size << std::endl;
	
	
	
	DSMI criterion(solver);	
	
	solver->add(new EvalBacktrack(&criterion));
	solver->add(new EvalSolution(&criterion));
	

	/** Save before making decisions **/
	solver->save();
	
	/** Beware of limits **/
	solver->parameters.restart_limit = 0;
	/** To record the max depth **/
	solver->statistics.max_depth = 0;
	/** The branch **/
	for(auto i{0}; i<3; ++i) {
		auto v{randint(2)};
		std::cout << v;
		Decision decision(X[i], Decision::ASSIGNMENT, v);
	  solver->reason_for[X[i].id()] = NULL;
		decision.make();
	}
	/** solver->stats->upper_bound is in the class SchedulingSolver, not Solver, so we can use it for the overall UB **/
	solver->set_objective(solver->stats->upper_bound);
	
	std::cout << std::endl << "DECISIONS " << solver->level << "|" << solver->decisions.size << std::endl;
	
	for(auto i{0}; i<5; ++i) {
		auto x{X[i]};
		std::cout 
			<< solver->sequence.contain(x) << " "
			<< x << " in " << x.get_domain()
			// << " / " << solver->disjuncts[solver->sequence[i].id() - offset] 
			<< ": " 
			<< solver->disjunct_map[x.id() - offset][1] << 
				" in " << solver->disjunct_map[x.id() - offset][1].get_domain() << " <> "
				<< solver->disjunct_map[x.id() - offset][2] << " in " 
					<< solver->disjunct_map[x.id() - offset][2].get_domain() << "\n";
	}
	
	/** propagate the branch and the bound **/
	solver->propagate();
	
	std::cout << "PROPAGATE " << solver->level << "|" << solver->decisions.size << std::endl;
	
	for(auto i{0}; i<5; ++i) {
		auto x{X[i]};
		std::cout 
			<< solver->sequence.contain(x) << " "
			<< x << " in " << x.get_domain()
			// << " / " << solver->disjuncts[solver->sequence[i].id() - offset] 
			<< ": " 
			<< solver->disjunct_map[x.id() - offset][1] << 
				" in " << solver->disjunct_map[x.id() - offset][1].get_domain() << " <> "
				<< solver->disjunct_map[x.id() - offset][2] << " in " 
					<< solver->disjunct_map[x.id() - offset][2].get_domain() << "\n";
	}
	
	/** set a fail limit **/
	auto numfails{solver->statistics.num_failures};
	solver->parameters.fail_limit = numfails + 1000;
	
	/** it might be possible to do restarts here **/
	auto root_level{solver->level};
	auto outcome = solver->chronological_dfs(root_level);
	
	switch(outcome) {
		case SAT: std::cout << "SOLUTION FOUND! (" << solver->get_objective() << ")\n"; break;
		case UNSAT: std::cout << "UNSAT!\n"; break;
		case LIMITOUT: std::cout << "LIMIT EXPIRED, MAX DEPTH = " << solver->statistics.max_depth << "\n"; break;
		case OPT: std::cout << "OPTIMAL\n"; break;
		default:  std::cout << "UNKNOWN\n"; break;
	}
	std::cout << "NUM FAILS = " << (solver->statistics.num_failures - numfails) << std::endl;
	
	std::cout << "AFTER SEARCH " << solver->level << "|" << solver->decisions.size << std::endl;
	
	for(auto i{0}; i<5; ++i) {
		auto x{X[i]};
		std::cout 
			<< solver->sequence.contain(x) << " "
			<< x << " in " << x.get_domain()
			// << " / " << solver->disjuncts[solver->sequence[i].id() - offset] 
			<< ": " 
			<< solver->disjunct_map[x.id() - offset][1] << 
				" in " << solver->disjunct_map[x.id() - offset][1].get_domain() << " <> "
				<< solver->disjunct_map[x.id() - offset][2] << " in " 
					<< solver->disjunct_map[x.id() - offset][2].get_domain() << "\n";
	}

	solver->restore(root_level);
	/** the method above cannot restore below the root level **/
	solver->restore();
	
	std::cout << "AFTER RESTORE " << solver->level << "|" << solver->decisions.size << std::endl;
	
	for(auto i{0}; i<5; ++i) {
		auto x{X[i]};
		std::cout 
			<< solver->sequence.contain(x) << " "
			<< x << " in " << x.get_domain()
			// << " / " << solver->disjuncts[solver->sequence[i].id() - offset] 
			<< ": " 
			<< solver->disjunct_map[x.id() - offset][1] << 
				" in " << solver->disjunct_map[x.id() - offset][1].get_domain() << " <> "
				<< solver->disjunct_map[x.id() - offset][2] << " in " 
					<< solver->disjunct_map[x.id() - offset][2].get_domain() << "\n";
	}
	


#ifdef _PROFILING
  std::cout << solver->statistics.total_propag_time << std::endl;
#endif


}
  




