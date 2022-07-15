 

#include "mistral_scheduler.hpp"


using namespace Mistral;


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
	
	auto ub{jsp.getMakespanUpperBound(10)};
	
	
	std::cout << ub << std::endl;
	
	// exit(1);


	Solver solver;
	
	VarArray start_time(jsp.nTasks(), 0, ub);
	VarArray end_time(jsp.nTasks(), 0, ub);
	
	
	for(auto i{0}; i<jsp.nTasks(); ++i) {
		solver.add(start_time[i] + jsp.getDuration(i) <= end_time[i]);
	}


	for(auto j{0}; j<jsp.nJobs(); ++j) {
		for(auto i{1}; i<jsp.nTasksInJob(j); ++i) {
			solver.add(start_time[jsp.getJobTask(j,i)] >= end_time[jsp.getJobTask(j,i-1)]);
		}
	}
	
	
	
	

	// Vector<VarArray> start_before;
	// Vector<VarArray> end_after;
	// Vector<VarArray> active_at;
	//
	//
	// for(auto i{0}; i<jsp.nTasks(); ++i) {
	//
	//
	// 	std::cout << "vars for task " << i << std::endl;
	//
	// 	VarArray si(ub);
	// 	start_before.push_back(si);
	// 	VarArray ei(ub);
	// 	end_after.push_back(ei);
	// 	VarArray ai(ub);
	// 	active_at.push_back(ai);
	//
	// }
	//
	// for(auto i{0}; i<jsp.nTasks(); ++i) {
	//
	//
	// 	std::cout << "cons for task " << i << std::endl;
	//
	// 	for(auto t{0}; t<ub; ++t) {
	// 		if(t>0) {
	// 			solver.add(start_before[i][t-1] <= start_before[i][t]);
	// 			solver.add(end_after[i][t] <= end_after[i][t-1]);
	// 		}
	// 		solver.add(active_at[i][t] <= (start_before[i][t] && end_after[i][t]));
	// 	}
	//
	// 	solver.add(Sum(active_at[i]) == jsp.getDuration(i));
	//
	//
	// 	std::cout << solver.constraints.size << std::endl;
	//
	// }
	
	VarArray end_job;
	for(auto j{0}; j<jsp.nJobs(); ++j) {
		auto i{jsp.nTasksInJob(j)};
		end_job.add(end_time[i]);
	}

	
	solver.minimize(Max(end_job));

	

  // SchedulingSolver solver(model, &params, &stats);
  usrand(params.Seed);
  
  solver.consolidate();

  std::cout << solver << std::endl;
  
	
	BranchingHeuristic *heuristic;
	RestartPolicy *restart;


	restart	= new Luby();
	restart->base = 128;
	heuristic = new LastConflict < GenericDVO < MinDomainOverWeight, 1, ConflictCountManager >, SolutionGuided< MinValue, RandomMinMax >, SolutionGuided< MinValue, RandomMinMax >, 1 > (&solver);


	solver.depth_first_search(solver.variables, heuristic, restart);
	
	// solver.depth_first_search();
	
	
	for(auto j{0}; j<jsp.nJobs(); ++j) {
		for(auto i{0}; i<jsp.nTasksInJob(j); ++i) {
			std::cout << " t" << (i+1) << ": [" << start_time[jsp.getJobTask(j,i)].get_solution_int_value() 
				<< ".." << end_time[jsp.getJobTask(j,i)].get_solution_int_value() << "]" ;
		}
		std::cout << std::endl;
	}

}
  




