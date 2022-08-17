

#include "mistral_scheduler.hpp"

#define VERBOSE

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

  // std::cout << std::endl;
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

  // for (auto i : info) {
  //   auto o{i.index};
  //   std::cout << " " << o << " in [" << start_time[o].get_min() << ".."
  //             << end_time[o].get_max() << "] (p=" << jsp.getDuration(o)
  //             << ")\n";
  // }

  BinaryMinHeap<TaskInfo> H;

  int t{0};
  // int curdue{std::numeric_limits<int>::infinity()};
  for (auto dd : info) {

    auto next{start_time[dd.index].get_min()};

    // std::cout << t << " -> " << next << "  " << dd.index << ": "
    //           << start_time[dd.index].get_min() << ".."
    //           << end_time[dd.index].get_max() << " ("
    //           << jsp.getDuration(dd.index) << ")" << std::endl;

    while (next > t) {

      // there is no task to run
      if (H.size() == 0) {

#ifdef VERBOSE
        std::cout << "idle from " << t << " to " << next << std::endl;
#endif

        t = next;
        std::pair<int, int> I{-1, next};
        sched.push_back(I);
      } else {

        // std::cout << H << std::endl;

        // std::cout << " ==> " << (H.get_min().index) << std::endl;

        assert(H.get_min().duration > 0);

        auto r{start_time[H.get_min().index].get_min()};

        if (r > t) {

          assert(r <= next);

#ifdef VERBOSE
          std::cout << "idle from " << t << " to " << r << std::endl;
#endif

          std::pair<int, int> I{-1, r};
          sched.push_back(I);

          t = r;
        } else {
          auto &urgent_task{H.get_min()};
          auto idx{urgent_task.index};

          if (t + urgent_task.duration <= next) {

#ifdef VERBOSE
            std::cout << "run " << urgent_task.index
                      << " until completion (t=" << (t + urgent_task.duration)
                      << ")\n";
#endif

            t += urgent_task.duration;
            urgent_task.setDuration(0);

            // std::pair<int, int> I{H.get_min().index, t};
            // sched.push_back(I);

            H.pop_min();
          } else {
            // auto incr{next - t}
            auto p{urgent_task.duration + t - next};

#ifdef VERBOSE
            std::cout << "run " << urgent_task.index << " for " << (next - t)
                      << " unit of time (t=" << next << ", " << p
                      << " remaining)\n";
#endif

            urgent_task.setDuration(p);
            t = next;
          }

#ifdef VERBOSE
          std::cout << "insert <" << idx << ", " << t << ">\n";
#endif

          if (sched.empty() or idx != sched.back().first) {
            std::pair<int, int> I{idx, t};
            sched.push_back(I);
          } else {
            sched.back().second = t;
          }
        }
      }
    }
    H.add(dd);
  }

  while (H.size() > 0) {
    auto a{H.pop_min()};

    assert(start_time[a.index].get_min() <= t);

    // std::cout << "run " << a.index
    //           << " until completion (t=" << (t + a.duration) << ")\n";
    t += a.duration;

    std::cout << "insert <" << a.index << ", " << t << ">\n";
    if (a.index != sched.back().first) {
      std::pair<int, int> I{a.index, t};
      sched.push_back(I);
    } else {
      sched.back().second = t;
    }

    // std::pair<int, int> I{a.index, t};
    // sched.push_back(I);
    assert(t <= a.date);
  }
}

bool checkPreemptiveScheduleSolution(const std::vector<int> &tasks,
                                     Instance &jsp, Solver &solver,
                                     VarArray &start_time, VarArray &end_time,
                                     std::vector<std::pair<int, int>> &sched) {

  std::cout << std::endl;
  std::vector<TaskInfo> info;
  for (auto a : tasks) {
    info.push_back(
        TaskInfo(a, end_time[a].get_solution_int_value(), jsp.getDuration(a)));
  }

  std::sort(info.begin(), info.end(),
            [&](const TaskInfo &a, const TaskInfo &b) {
              return start_time[a.index].get_solution_int_value() <
                         start_time[b.index].get_solution_int_value() or
                     (start_time[a.index].get_solution_int_value() ==
                          start_time[b.index].get_solution_int_value() and
                      end_time[a.index].get_solution_int_value() <
                          end_time[b.index].get_solution_int_value());
            });

  // for (auto i : info) {
  //   auto o{i.index};
  //   std::cout << " " << o << " in [" <<
  //   start_time[o].get_solution_int_value()
  //             << ".." << end_time[o].get_solution_int_value()
  //             << "] (p=" << jsp.getDuration(o) << ")\n";
  // }

  BinaryMinHeap<TaskInfo> H;

  int t{0};
  // int curdue{std::numeric_limits<int>::infinity()};
  for (auto dd : info) {

    auto next{start_time[dd.index].get_solution_int_value()};

#ifdef VERBOSE
    std::cout << t << " -> " << next << "  " << dd.index << ": "
              << start_time[dd.index].get_solution_int_value() << ".."
              << end_time[dd.index].get_solution_int_value() << " ("
              << jsp.getDuration(dd.index) << ")" << std::endl;
#endif

    while (next > t) {

      // there is no task to run
      if (H.size() == 0) {
        // std::cout << "idle from " << t << " to " << next << std::endl;
        t = next;
      } else {

        // std::cout << H << std::endl;

        // std::cout << " ==> " << (H.get_min().index) << std::endl;

        assert(H.get_min().duration > 0);

        auto r{start_time[H.get_min().index].get_solution_int_value()};

        if (r > t) {

          assert(r <= next);

#ifdef VERBOSE
          std::cout << "idle from " << t << " to " << r << std::endl;
#endif

          t = r;
        } else {
          if (t + H.get_min().duration <= next) {

#ifdef VERBOSE
            std::cout << "run " << H.get_min().index
                      << " until completion (t=" << (t + H.get_min().duration)
                      << ")\n";
#endif

            t += H.get_min().duration;

            if (t > H.get_min().date) {

#ifdef VERBOSE
              std::cout << t << " > " << H.get_min().index << "'s end ("
                        << H.get_min().date << ")\n";
#endif

              return false;
            }

            H.get_min().setDuration(0);
            H.pop_min();
          } else {
            // auto incr{next - t}
            auto p{H.get_min().duration + t - next};

#ifdef VERBOSE
            std::cout << "run " << H.get_min().index << " for " << (next - t)
                      << " unit of time (t=" << next << ", " << p
                      << " remaining)\n";
#endif

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

#ifdef VERBOSE
    std::cout << "run " << a.index
              << " until completion (t=" << (t + a.duration) << ")\n";
#endif

    t += a.duration;

    if (t > a.date) {

#ifdef VERBOSE
      std::cout << t << " > " << a.index << "'s end (" << a.date << ")\n";
#endif

      return false;
    }
  }

  return true;
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

  solver.parameters.verbosity = 3;

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

  // VarArray ordering;
  // for (auto k{0}; k < jsp.nMachines(); ++k) {
  // 		for (auto i{0}; i < jsp.nTasksInMachine(k); ++i) {
  // 			for (auto j{i+1}; j < jsp.nTasksInMachine(k); ++j) {
  // 				ordering.add(ReifiedDisjunctive(start_time[jsp.getMachineTask(k,i)],
  // start_time[jsp.getMachineTask(k,j)], 1, 1));
  // 				ordering.add(ReifiedDisjunctive(start_time[jsp.getMachineTask(k,i)],
  // end_time[jsp.getMachineTask(k,j)], 1, 1));
  // 				ordering.add(ReifiedDisjunctive(end_time[jsp.getMachineTask(k,i)],
  // start_time[jsp.getMachineTask(k,j)], 1, 1));
  // 				ordering.add(ReifiedDisjunctive(end_time[jsp.getMachineTask(k,i)],
  // end_time[jsp.getMachineTask(k,j)], 1, 1));
  // 			}
  // 		}
  // }

  // for(auto &d : ordering)
  // 	solver.add(Free(d));

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

  Variable makespan{Max(end_job)};

  solver.minimize(makespan);

  // SchedulingSolver solver(model, &params, &stats);
  usrand(params.Seed);

  std::vector<std::pair<int, int>> intervals;

  solver.consolidate();

  // std::cout << solver << std::endl;

  BranchingHeuristic *heuristic;
  RestartPolicy *restart;

  restart = new Luby();
  restart->base = 128;
  heuristic =
      new LastConflict<GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
                       SolutionGuided<MinValue, MinValue>,
                       SolutionGuided<MinValue, MinValue>, 1>(&solver);

  solver.initialise_search(solver.variables, heuristic, restart);

  solver.propagate();

  // exit(1);

  int LB = 0;
  for (auto k{0}; k < jsp.nMachines(); ++k) {
    intervals.clear();
    computePreemptiveSchedule(jsp.getMachineTasks(k), jsp, solver, start_time,
                              end_time, intervals);
    for (auto p : intervals) {
      if (p.first >= 0)
        std::cout << "run " << p.first << " until " << p.second << std::endl;
      else
        std::cout << "idle until " << p.second << std::endl;
    }

    auto trail{ub - end_time[intervals.back().first].get_max()};
    std::cout << "LB = " << intervals.back().second << " + " << trail << " = "
              << (intervals.back().second + trail) << std::endl;

    std::cout << std::endl;

    if (LB < (intervals.back().second + trail))
      LB = (intervals.back().second + trail);
  };

  solver.add(makespan >= LB);

  // solver.propagate();

  // exit(1);

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

  solver.depth_first_search();
  // solver.depth_first_search(solver.variables, heuristic, restart);
  // solver.depth_first_search(ordering, heuristic, restart);

  for (auto k{0}; k < jsp.nMachines(); ++k) {
    intervals.clear();
    if (not checkPreemptiveScheduleSolution(jsp.getMachineTasks(k), jsp, solver,
                                            start_time, end_time, intervals)) {

      std::cout << "Error on machine " << k << std::endl;
      exit(1);
    }
  };

  std::cout << "c solution checked!\n";

  //   for (auto k{0}; k < jsp.nMachines(); ++k) {
  //     for (auto i{0}; i < jsp.nTasksInMachine(k); ++i) {
  // 		auto x{jsp.getMachineTask(k, i)};
  // 		std::cout << "t" << x << ": [" << start_time[x].get_min() <<
  // ".."
  // << end_time[x].get_max() << "]\n";
  // 	}
  // 	std::cout << std::endl;
  // }

  // // solver.depth_first_search();

  // print_solution(jsp, solver, start_time, end_time);
}
