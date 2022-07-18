

#include "mistral_scheduler.hpp"


using namespace Mistral;

void print_solution(Instance &jsp, Solver &solver, VarArray &start_time,
                    VarArray &end_time) {

  for (auto j{0}; j < jsp.nJobs(); ++j) {
    for (auto i{0}; i < jsp.nTasksInJob(j); ++i) {
      std::cout << " t" << (i + 1) << ": ["
                << start_time[jsp.getJobTask(j, i)].get_solution_int_value()
                << ".."
                << end_time[jsp.getJobTask(j, i)].get_solution_int_value()
                << "]";
    }
    std::cout << std::endl;
  }
}

class TaskInfo {

public:
  int index;
  int date;
  int duration;

  void setDuration(const int p) { duration = p; }

  TaskInfo() : index(0), date(0), duration(0) {}

  TaskInfo(const int t, const int d, const int p)
      : index(t), date(d), duration(p) {}

  bool operator<=(const TaskInfo &lct) const { return date <= lct.date; }

  bool operator<(const TaskInfo &lct) const { return date < lct.date; }

  std::ostream &display(std::ostream &os) const {
    os << index;
    return os;
  }
};

std::ostream &operator<<(std::ostream &os, const TaskInfo &x) {
  return x.display(os);
}

void computePreemptiveSchedule(const std::vector<int> &tasks, Instance &jsp,
                               Solver &solver, VarArray &start_time,
                               VarArray &end_time,
                               std::vector<std::pair<int, int>> &sched) {

  // for (auto a : tasks) {
  //   std::cout << " " << a;
  // }
  // std::cout << std::endl;

  // std::vector<ThetaElement> leaves;
  // int id{0};
  // for (auto a : tasks) {
  //   leaves.push_back(ThetaElement(id++, start_time[a].get_min(),
  //   end_time[a].get_min(), end_time[a].get_max(), jsp.getDuration(a)));
  // }

  // for (auto i : leaves) {
  //   std::cout << " t" << i << ": [" << start_time[i].get_min() << ".."
  //             << end_time[i].get_max() << "]\n";
  // }


																 std::cout << std::endl;
  std::vector<TaskInfo> info;
  for (auto a : tasks) {
    info.push_back(TaskInfo(a, end_time[a].get_max(), jsp.getDuration(a)));
  }

  std::sort(
      info.begin(), info.end(), [&](const TaskInfo &a, const TaskInfo &b) {
        return start_time[a.index].get_min() < start_time[b.index].get_min() or
               (start_time[a.index].get_min() ==
                    start_time[b.index].get_min() and
                end_time[a.index].get_max() < end_time[b.index].get_max());
      });

  for (auto i : info) {
    auto o{i.index};
    std::cout << " " << o << " in [" << start_time[o].get_min() << ".."
              << end_time[o].get_max() << "] (p=" << jsp.getDuration(o) << ")\n";
  }
  

  BinaryMinHeap<TaskInfo> H;

  int t{0};
  // int curdue{std::numeric_limits<int>::infinity()};
  for (auto dd : info) {

    auto next{start_time[dd.index].get_min()};

    std::cout << t << " -> " << next << "  " << dd.index << ": "
              << start_time[dd.index].get_min() << ".."
              << end_time[dd.index].get_max() << " ("
              << jsp.getDuration(dd.index) << ")" << std::endl;


    while (next > t) {

      // there is no task to run
      if (H.size() == 0) {
        std::cout << "idle from " << t << " to " << next << std::endl;
        t = next;
      } else {

        // std::cout << H << std::endl;

        // std::cout << " ==> " << (H.get_min().index) << std::endl;

        assert(H.get_min().duration > 0);

        auto r{start_time[H.get_min().index].get_min()};

        if (r > t) {

          assert(r <= next);

          std::cout << "idle from " << t << " to " << r << std::endl;
          t = r;
        } else {
          if (t + H.get_min().duration <= next) {
            std::cout << "run " << H.get_min().index
                      << " until completion (t=" << (t + H.get_min().duration)
                      << ")\n";
            t += H.get_min().duration;
            H.get_min().setDuration(0);
            H.pop_min();
          } else {
            // auto incr{next - t}
            auto p{H.get_min().duration + t - next};
            std::cout << "run " << H.get_min().index << " for " << (next - t)
                      << " unit of time (t=" << next << ", " << p
                      << " remaining)\n";
            H.get_min().setDuration(p);
            t = next;
          }
        }
      }
    }
    H.add(dd);
  }

  while (H.size() > 0) {
    auto a{H.pop_min()};

    assert(start_time[a.index].get_min() <= t);

    std::cout << "run " << a.index
              << " until completion (t=" << (t + a.duration) << ")\n";
    t += a.duration;
    assert(t <= a.date);
  }
}

int main( int argc, char** argv )
{

  ParameterList params(argc, argv);
  usrand(params.Seed);

  StatisticList stats;
  stats.start();

  Instance jsp(params);

  std::cout << std::endl;

  // jsp.print(std::cout);

  // jsp.printStats(std::cout);
  // params.print(std::cout);

  auto ub{jsp.getMakespanUpperBound(10)};

  // std::cout << ub << std::endl;

  // exit(1);

  Solver solver;
	
	solver.parameters.verbosity = 2;

  VarArray start_time(jsp.nTasks(), 0, ub);
  VarArray end_time(jsp.nTasks(), 0, ub);

  for (auto i{0}; i < jsp.nTasks(); ++i) {
    solver.add(start_time[i] + jsp.getDuration(i) <= end_time[i]);
  }

  for (auto j{0}; j < jsp.nJobs(); ++j) {
    for (auto i{1}; i < jsp.nTasksInJob(j); ++i) {
      solver.add(start_time[jsp.getJobTask(j, i)] >=
                 end_time[jsp.getJobTask(j, i - 1)]);
    }
  }
	
	
	VarArray ordering;
	for (auto k{0}; k < jsp.nMachines(); ++k) {
			for (auto i{0}; i < jsp.nTasksInMachine(k); ++i) {
				for (auto j{i+1}; j < jsp.nTasksInMachine(k); ++j) {
					ordering.add(ReifiedDisjunctive(start_time[jsp.getMachineTask(k,i)], start_time[jsp.getMachineTask(k,j)], 1, 1));
					ordering.add(ReifiedDisjunctive(start_time[jsp.getMachineTask(k,i)], end_time[jsp.getMachineTask(k,j)], 1, 1));
					ordering.add(ReifiedDisjunctive(end_time[jsp.getMachineTask(k,i)], start_time[jsp.getMachineTask(k,j)], 1, 1));
					ordering.add(ReifiedDisjunctive(end_time[jsp.getMachineTask(k,i)], end_time[jsp.getMachineTask(k,j)], 1, 1));
				}
			}
	}

	for(auto &d : ordering)
		solver.add(Free(d));



  std::vector<int> p;
  for (auto k{0}; k < jsp.nMachines(); ++k) {
    VarArray st;
    VarArray et;
    p.clear();
    for (auto i{0}; i < jsp.nTasksInMachine(k); ++i) {
      auto o{jsp.getMachineTask(k, i)};
      st.add(start_time[o]);
      et.add(end_time[o]);
      p.push_back(jsp.getDuration(o));
    }

    solver.add(NoOverlap(st, et, p));
  }

  VarArray end_job;
  for (auto j{0}; j < jsp.nJobs(); ++j) {
    // auto i{jsp.nTasksInJob(j)-1};
    end_job.add(end_time[jsp.getLastTaskofJob(j)]);
  }

  solver.minimize(Max(end_job));

  // SchedulingSolver solver(model, &params, &stats);
  usrand(params.Seed);

  solver.consolidate();

  // std::cout << solver << std::endl;

  BranchingHeuristic *heuristic;
  RestartPolicy *restart;

  restart = new Luby();
  restart->base = 128;
  heuristic =
      new LastConflict<GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
                       SolutionGuided<MinValue, RandomMinMax>,
                       SolutionGuided<MinValue, RandomMinMax>, 1>(&solver);

  // solver.initialise_search(solver.variables, heuristic, restart);

  // solver.propagate();
  //
  // // for (auto k{0}; k < 1; ++k) {
  // //   for (auto i{0}; i < jsp.nTasksInMachine(k); ++i) {
  // //     auto o{jsp.getMachineTask(k, i)};
  // //     std::cout << " " << o << " in [" << start_time[o].get_min() << ".."
  // //               << end_time[o].get_max() << "]\n";
  // //   }
  // //   std::cout << std::endl;
  // // }
  //
  // std::vector<std::pair<int, int>> intervals;
  //
  // computePreemptiveSchedule(jsp.getMachineTasks(0), jsp, solver, start_time,
  //                           end_time, intervals);

  // solver.depth_first_search(solver.variables, heuristic, restart);
	solver.depth_first_search(ordering, heuristic, restart);
	
	
  std::vector<std::pair<int, int>> intervals;

	for (auto k{0}; k < jsp.nMachines(); ++k) {
		intervals.clear();
  computePreemptiveSchedule(jsp.getMachineTasks(k), jsp, solver, start_time,
                            end_time, intervals);
													};
	
	//   for (auto k{0}; k < jsp.nMachines(); ++k) {
	//     for (auto i{0}; i < jsp.nTasksInMachine(k); ++i) {
	// 		auto x{jsp.getMachineTask(k, i)};
	// 		std::cout << "t" << x << ": [" << start_time[x].get_min() << ".." << end_time[x].get_max() << "]\n";
	// 	}
	// 	std::cout << std::endl;
	// }
	
	

  // // solver.depth_first_search();

  // print_solution(jsp, solver, start_time, end_time);
}
