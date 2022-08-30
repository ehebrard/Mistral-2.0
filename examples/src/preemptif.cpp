

#include "mistral_scheduler.hpp"

// #define VERBOSE

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
  int release_date;
  int due_date;
  int duration;

  void setDuration(const int p) { duration = p; }

  TaskInfo() : index(0), release_date(0), due_date(0), duration(0) {}

  TaskInfo(const int t, const int rd, const int dd, const int p)
      : index(t), release_date(rd), due_date(dd), duration(p) {}

  bool operator<=(const TaskInfo &lct) const { return due_date <= lct.due_date; }

  bool operator<(const TaskInfo &lct) const { return due_date < lct.due_date; }

  std::ostream &display(std::ostream &os) const {
    os << index;
    return os;
  }
};

std::ostream &operator<<(std::ostream &os, const TaskInfo &x) {
  return x.display(os);
}

class JacksonPreemptiveScheduler {

public:
  std::vector<std::pair<int, int>> sched;
  std::vector<int> sequence;

  JacksonPreemptiveScheduler(Instance &jsp, Solver &solver,
                             VarArray &start_time, VarArray &end_time)
      : jsp(jsp), solver(solver), start_time(start_time), end_time(end_time) {
    started.resize(jsp.nTasks(), false);
  };

  bool fromDomain(const std::vector<int> &tasks);
  bool fromSolution(const std::vector<int> &tasks);

  int get_lower_bound(const int ub);

private:
  Instance &jsp;
  Solver &solver;
  VarArray &start_time;
  VarArray &end_time;

  std::vector<TaskInfo> info;
  // std::vector<int> tasks;

  std::vector<bool> started;

  bool compute();
};

int JacksonPreemptiveScheduler::get_lower_bound(const int ub) {

#ifdef VERBOSE
  std::cout << "lb\n";
#endif

  solver.propagate();

  int LB = 0;
  for (auto k{0}; k < jsp.nMachines(); ++k) {

    // for (auto a : jsp.getMachineTasks(k))
    //   std::cout << " " << a;
    // std::cout << std::endl;

    fromDomain(jsp.getMachineTasks(k));

#ifdef VERBOSE
    for (auto p : sched) {
      if (p.first >= 0)
        std::cout << "run " << p.first << " until " << p.second << std::endl;
      else
        std::cout << "idle until " << p.second << std::endl;
    }
#endif

    for (auto a : sequence)
      std::cout << " " << a;
    

    auto trail{ub - end_time[sched.back().first].get_max()};

#ifdef VERBOSE
    std::cout << "LB = " << sched.back().second << " + " << trail << " = "
              << (sched.back().second + trail) << std::endl;

    std::cout << std::endl;
#endif

    if (LB < (sched.back().second + trail))
      LB = (sched.back().second + trail);


    std::cout << " => " << LB << std::endl;
  };

  return LB;
}

bool JacksonPreemptiveScheduler::fromDomain(const std::vector<int> &tasks) {

  info.clear();
  for (auto a : tasks) {
    info.push_back(TaskInfo(a, start_time[a].get_min(), end_time[a].get_max(),
                            jsp.getDuration(a)));
  }

  return compute();
}

bool JacksonPreemptiveScheduler::fromSolution(const std::vector<int> &tasks) {

  info.clear();
  for (auto a : tasks) {
    info.push_back(TaskInfo(a, start_time[a].get_solution_int_value(),
                            end_time[a].get_solution_int_value(),
                            jsp.getDuration(a)));
  }

  return compute();
}

bool JacksonPreemptiveScheduler::compute() {

#ifdef VERBOSE
  std::cout << "JSP: sort\n";
#endif

  sched.clear();
  sequence.clear();

  std::sort(
      info.begin(), info.end(), [&](const TaskInfo &a, const TaskInfo &b) {
        return a.release_date < b.release_date or
               (a.release_date == b.release_date and a.due_date < b.due_date);
      });

#ifdef VERBOSE
  for (auto i : info) {
    std::cout << " " << i.index << " in [" << i.release_date << ".."
              << i.due_date << "] (p=" << i.duration << ")\n";
  }
#endif

  BinaryMinHeap<TaskInfo> H;

  int t{0};
  for (auto dd : info) {

    auto next{dd.release_date};

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

        assert(H.get_min().duration > 0);

        auto r{H.get_min().release_date};

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
            sequence.push_back(idx);

            H.pop_min();
          } else {
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

          if (not started[idx]) {
            sequence.push_back(idx);
            started[idx] = true;
          }
        }
      }
    }
    H.add(dd);
  }

  while (H.size() > 0) {
    auto a{H.pop_min()};

#ifdef VERBOSE
    std::cout << "run " << a.index
              << " until completion (t=" << (t + a.duration) << ")\n";
#endif

    auto idx{a.index};
    if (not started[idx]) {
      sequence.push_back(idx);
      started[a.index] = true;
    }
    sequence.push_back(idx);

    t += a.duration;

    if (t > a.due_date) {

#ifdef VERBOSE
      std::cout << t << " > " << idx << "'s end (" << a.due_date << ")\n";
#endif

      return false;
    }

    if (a.index != sched.back().first) {
      std::pair<int, int> I{idx, t};
      sched.push_back(I);
    } else {
      sched.back().second = t;
    }
  }

  return true;
}

void model(Instance& jsp, Solver& solver, VarArray& start_time, VarArray& end_time, Variable& makespan) {
#ifdef VERBOSE
  std::cout << "model\n";
#endif

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
  //    for (auto i{0}; i < jsp.nTasksInMachine(k); ++i) {
  //      for (auto j{i+1}; j < jsp.nTasksInMachine(k); ++j) {
  //        ordering.add(ReifiedDisjunctive(start_time[jsp.getMachineTask(k,i)],
  // start_time[jsp.getMachineTask(k,j)], 1, 1));
  //        ordering.add(ReifiedDisjunctive(start_time[jsp.getMachineTask(k,i)],
  // end_time[jsp.getMachineTask(k,j)], 1, 1));
  //        ordering.add(ReifiedDisjunctive(end_time[jsp.getMachineTask(k,i)],
  // start_time[jsp.getMachineTask(k,j)], 1, 1));
  //        ordering.add(ReifiedDisjunctive(end_time[jsp.getMachineTask(k,i)],
  // end_time[jsp.getMachineTask(k,j)], 1, 1));
  //      }
  //    }
  // }

  // for(auto &d : ordering)
  //  solver.add(Free(d));

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

    // for (auto i{0}; i < jsp.nTasksInMachine(k); ++i) {
    //   for (auto j{0}; j < jsp.nTasksInMachine(k); ++j) {
    //     if(i != j) {
    //       solver.add((st[i] <= st[j] && et[i] <= et[j]) <= (et[i] <= st[j]));
    //     }
    //   }
    // }

    solver.add(PreemptiveNoOverlap(st, et, p));
  }

  VarArray end_job;
  for (auto j{0}; j < jsp.nJobs(); ++j) {
    end_job.add(end_time[jsp.getLastTaskofJob(j)]);
  }

  solver.add(makespan == Max(end_job));
  

  solver.minimize(makespan);

  solver.consolidate();

}


void set_strategy(Solver& solver) {

  BranchingHeuristic *heuristic;
  RestartPolicy *restart;

  restart = new Geometric();
  restart->base = 128;
  heuristic =
      new LastConflict<GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
                       SolutionGuided<MinValue, MinValue>,
                       SolutionGuided<MinValue, MinValue>, 1>(&solver);

  solver.initialise_search(solver.variables, heuristic, restart);
}

int main( int argc, char** argv )
{

  ParameterList params(argc, argv);
  usrand(params.Seed);

  StatisticList stats;
  stats.start();

#ifdef VERBOSE
  std::cout << "read\n";
#endif

  Instance jsp(params);


#ifdef VERBOSE
  std::cout << "ub\n";
#endif

  auto ub{jsp.getMakespanUpperBound(10)};

  if(params.UBinit >= 0)
    ub = params.UBinit;


  Solver solver;
  solver.parameters.verbosity = params.Verbose;

  VarArray start_time(jsp.nTasks(), 0, ub);
  VarArray end_time(jsp.nTasks(), 0, ub);

  Variable makespan(0,ub);

  model(jsp, solver, start_time, end_time, makespan);

  JacksonPreemptiveScheduler JPS(jsp, solver, start_time, end_time);

  set_strategy(solver);

  // auto lb{get_lower_bound(jsp, solver, start_time, end_time, ub)};
  auto lb{JPS.get_lower_bound(ub)};

  solver.add(makespan >= lb);

  std::cout << " c initial bounds: [" << lb << ".." << ub << "]\n";

  solver.depth_first_search();

  // std::vector<std::pair<int, int>> intervals;
  // for (auto k{0}; k < jsp.nMachines(); ++k) {
  //   intervals.clear();
  //   if (not checkPreemptiveScheduleSolution(jsp.getMachineTasks(k), jsp,
  //   solver,
  //                                           start_time, end_time, intervals))
  //                                           {
  //     std::cout << "Error on machine " << k << std::endl;
  //     exit(1);
  //   }
  // };

  for (auto k{0}; k < jsp.nMachines(); ++k) {
    if (not JPS.fromSolution(jsp.getMachineTasks(k))) {
      std::cout << "Error on machine " << k << std::endl;
      exit(1);
    }
  };

  std::cout << " c solution checked!\n";
}
