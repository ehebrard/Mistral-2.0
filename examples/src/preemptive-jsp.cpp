

#include "mistral_scheduler.hpp"

// #define VERBOSE

enum model { basic = 0, equalities = 1, compact = 2 };

enum branching {
  minvalue = 0,
  halfsplit = 1,
  reversesplit = 2,
  randomsplit = 3,
  randompivot = 4,
  guidedsplit = 5
};

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
    os << due_date;
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
  // std::vector<bool> seen;

  JacksonPreemptiveScheduler(Instance &jsp, Solver &solver,
                             VarArray &start_time, VarArray &end_time)
      : jsp(jsp), solver(solver), start_time(start_time), end_time(end_time) {
    started.resize(jsp.nTasks(), false);
    // seen.reserve(jsp.nTasks());
  };

  bool fromDomain(const std::vector<int> &tasks);
  bool fromSolution(const std::vector<int> &tasks);
  void printSolution(const std::vector<int> &tasks);

  int get_lower_bound(const int ub);

  void print_machine() {

  for (auto i : info) {
    std::cout << " " << i.index << " in [" << i.release_date << ".."
              << i.due_date << "] (p=" << i.duration << ")\n";
  }
}

private:
  Instance &jsp;
  Solver &solver;
  VarArray &start_time;
  VarArray &end_time;

  std::vector<TaskInfo> info;
  // std::vector<int> tasks;

  std::vector<bool> started;

  bool compute(const bool print_flag=false);
  void print();
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
    

    auto trail{ub - end_time[sched.back().first].get_max()};

#ifdef VERBOSE
    std::cout << "LB = " << sched.back().second << " + " << trail << " = "
              << (sched.back().second + trail) << std::endl;

    std::cout << std::endl;
#endif

    if (LB < (sched.back().second + trail))
      LB = (sched.back().second + trail);


    // std::cout << " => " << LB << std::endl;
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
  // std::cout << std::endl;
  for (auto a : tasks) {

    // std::cout << a << " " << jsp.getDuration(a) << " "
    // << start_time[a].get_solution_int_value() << " "
    // << end_time[a].get_solution_int_value() << std::endl;

    info.push_back(TaskInfo(a, start_time[a].get_solution_int_value(),
                            end_time[a].get_solution_int_value(),
                            jsp.getDuration(a)));
  }

  return compute();
}

void JacksonPreemptiveScheduler::printSolution(const std::vector<int> &tasks) {

  info.clear();
  // std::cout << std::endl;
  for (auto a : tasks) {

    // std::cout << a << " " << jsp.getDuration(a) << " "
    // << start_time[a].get_solution_int_value() << " "
    // << end_time[a].get_solution_int_value() << std::endl;

    info.push_back(TaskInfo(a, start_time[a].get_solution_int_value(),
                            end_time[a].get_solution_int_value(),
                            jsp.getDuration(a)));
  }

  print();
}

bool JacksonPreemptiveScheduler::compute(const bool print_flag) {

#ifdef VERBOSE
  std::cout << "JSP: sort\n";
#endif

  sched.clear();
  sequence.clear();
  // seen.clear();
  // seen.resize(jsp.nTasks(), false);

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

        // assert(H.get_min().duration > 0);

        auto r{H.get_min().release_date};


        // H.display(std::cout);
        // std::cout << "seen:";
        // for(auto i{0}; i<jsp.nTasks(); ++i) {
        //   if(seen[i])
        //     std::cout << " " << i;
        // }
        // std::cout << std::endl;


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
    // assert(not seen[dd.index]);
    // std::cout << "add " << dd.index << " to heap\n";
    H.add(dd);
    // seen[dd.index] = true;
  }


  // std::cout << "here\n"; 

  while (H.size() > 0) {

    // H.display(std::cout);

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


void JacksonPreemptiveScheduler::print() {
  
}

void build_model(Instance &jsp, Solver &solver, VarArray &start_time,
                 VarArray &end_time, Variable &origin, Variable &makespan,
                 VarArray &search_vars, const int ub, const int model_choice,
                 const bool hall_pruning) {

  if (model_choice == model::compact) {

#ifdef VERBOSE
    std::cout << "jobs\n";
#endif

    search_vars.initialise(jsp.nTasks() - jsp.nJobs(), 0, ub);

    for (auto i{0}; i < jsp.nTasks(); ++i) {
      auto j{jsp.getJob(i, 0)};
      auto r{jsp.getRankInJob(i)};
      if (r == 0)
        start_time.add(origin);
      else
        start_time.add(search_vars[i - j - 1]);

      if (r == jsp.nMachines() - 1)
        end_time.add(makespan);
      else
        end_time.add(search_vars[i - j]);
    }

  } else {

    start_time.initialise(jsp.nTasks(), 0, ub);
    end_time.initialise(jsp.nTasks(), 0, ub);
  }

#ifdef VERBOSE
  std::cout << "model\n\ndurations\n";
#endif

  for (auto i{0}; i < jsp.nTasks(); ++i) {
    solver.add(start_time[i] + jsp.getDuration(i) <= end_time[i]);
  }

  if (model_choice != model::compact) {

#ifdef VERBOSE
  std::cout << "jobs\n";
#endif

  for (auto j{0}; j < jsp.nJobs(); ++j) {

    if (model_choice == model::equalities) {
      solver.add(start_time[jsp.getJobTask(j, 0)] == origin);

      for (auto i{1}; i < jsp.nTasksInJob(j); ++i) {
        solver.add(start_time[jsp.getJobTask(j, i)] ==
                   end_time[jsp.getJobTask(j, i - 1)]);
      }

      solver.add(end_time[jsp.getJobTask(j, jsp.nTasksInJob(j) - 1)] ==
                 makespan);
    } else {

      solver.add(start_time[jsp.getJobTask(j, 0)] <= 5);

      for (auto i{1}; i < jsp.nTasksInJob(j); ++i) {
        solver.add(start_time[jsp.getJobTask(j, i)] >=
                   end_time[jsp.getJobTask(j, i - 1)]);

        solver.add((start_time[jsp.getJobTask(j, i)] - 5) <=
                   end_time[jsp.getJobTask(j, i - 1)]);
      }
    }
  }
  }

#ifdef VERBOSE
  std::cout << "resources\n";
#endif

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

    // std::cout << solver.constraints.size << std::endl;
    solver.add(PreemptiveNoOverlap(st, et, p, makespan, hall_pruning));
  }

#ifdef VERBOSE
  std::cout << "makespan\n";
#endif

  // VarArray end_job;
  // for (auto j{0}; j < jsp.nJobs(); ++j) {
  //   end_job.add(end_time[jsp.getLastTaskofJob(j)]);
  // }

  // solver.add(makespan == Max(end_job));

  if (model_choice == model::basic) {
    VarArray end_job;
    for (auto j{0}; j < jsp.nJobs(); ++j) {
      end_job.add(end_time[jsp.getLastTaskofJob(j)]);
    }

    solver.add(makespan == Max(end_job));
  }

  solver.minimize(makespan);

#ifdef VERBOSE
  std::cout << "consolidate\n";
#endif

  solver.consolidate();

  if (model_choice != model::compact) {
    for (auto s : start_time)
      search_vars.add(s);

    if (model_choice == model::basic) {
      for (auto e : end_time)
        search_vars.add(e);
    }
  }
  // else {
  //   search_vars.add(makespan);
  // }

#ifdef VERBOSE
  std::cout << "end model\n";
#endif
}

void model_order(Instance& jsp, Solver& solver, VarArray& start_time, VarArray& end_time, Variable& makespan, VarArray& search_vars) {
#ifdef VERBOSE
  std::cout << "model\n";
#endif

  search_vars.clear();

  VarArray ordering;
  for (auto k{0}; k < jsp.nMachines(); ++k) {
     for (auto i{0}; i < jsp.nTasksInMachine(k); ++i) {
       for (auto j{i+1}; j < jsp.nTasksInMachine(k); ++j) {
         ordering.add(ReifiedDisjunctive(end_time[jsp.getMachineTask(k,i)],
  end_time[jsp.getMachineTask(k,j)], 1, 1));
       }
     }
  }

  for(auto &d : ordering) {
    search_vars.add(d);
    solver.add(Free(search_vars.back()));
  }
}


RestartPolicy *restart_factory(std::string rpolicy, const int b, const double f) {
  RestartPolicy *pol;
  if(rpolicy == "luby") pol = new Luby(b); 
  else if(rpolicy == "geom") pol = new Geometric(b,f); 
  else pol = new NoRestart();
  return pol;
}

BranchingHeuristic *heuristic_factory(Solver &solver,
                                  const std::string &branching_choice) {

  BranchingHeuristic *heuristic;

  if (branching_choice == "minvalue")
    heuristic = new LastConflict<
        GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>, MinValue,
        MinValue, 1>(&solver);
  else if (branching_choice == "maxvalue")
    heuristic = new LastConflict<
        GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>, MaxValue,
        MaxValue, 1>(&solver);
  else if (branching_choice == "halfsplit")
    heuristic = new LastConflict<
        GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>, HalfSplit,
        HalfSplit, 1>(&solver);
  else if (branching_choice == "randompivot")
    heuristic = new LastConflict<
        GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
        RandomPivotSplit, RandomPivotSplit, 1>(&solver);
  else if (branching_choice == "randomsplit")
    heuristic = new LastConflict<
        GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>, RandomSplit,
        RandomSplit, 1>(&solver);
  else if (branching_choice == "reversesplit")
    heuristic = new LastConflict<
        GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>, ReverseSplit,
        ReverseSplit, 1>(&solver);
  else if (branching_choice == "guidedsplit")
    heuristic = new LastConflict<
        GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
        GuidedSplit<RandomSplit>, GuidedSplit<RandomSplit>, 1>(&solver);
  else {
    std::cout << " c Warning: heuristic " << branching_choice << " is not handled, using halfsplit instead\n";
     heuristic = new LastConflict<
        GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>, HalfSplit,
        HalfSplit, 1>(&solver);
  }
 
  return heuristic;
}

int main( int argc, char** argv )
{

  std::cout << " c Mistral preemptive jobshop scheduler" << std::endl;
  SolverCmdLine cmd("Mistral preemptive scheduler", ' ', "0.0");

  // TCLAP::SwitchArg order_branching("","order","Branches on the ordering of
  // end times", false); cmd.add( order_branching );

  TCLAP::SwitchArg hall_pruning("", "hall",
                                "Applies Hall Intervals-based Pruning", false);
  cmd.add(hall_pruning);

  TCLAP::ValueArg<int> model_choice(
      "", "model", "choice of model O:default 1:reduced 2:compact", false, 2,
      "int");
  cmd.add(model_choice);

  // TCLAP::ValueArg<int> branching_choice(
  //     "", "value-order",
  //     "choice of branching strategy O:minvalue 1:halfsplit 2:reversesplit "
  //     "3:randomsplit 4:randompivot 5:guidedsplit",
  //     false, 2, "int");
  // cmd.add(branching_choice);

  TCLAP::ValueArg<std::string> format(
      "", "format",
      "Instance file format (can be: jsp/jla/osp/dyn/jet/now/sds)", false,
      "jla", "string");
  cmd.add( format );

  TCLAP::ValueArg<int> init_ub(
      "", "ub", "Sets a manual upper bound (-1 stands for no upper bound)",
      false, -1, "int");
  cmd.add( init_ub );

  cmd.branchingArg->reset("halfsplit");
  cmd.verbosityArg->reset(3);

  cmd.parse(argc, argv);

  if (hall_pruning.getValue()) {
    std::cout << " c hall pruning on (basic model enforced)\n";
    // model_choice.setValue(model::basic);
  }

  std::cout << " c instance=" << cmd.get_filename() << " model="
            << ((model_choice.getValue() == model::basic or
                 hall_pruning.getValue())
                    ? "basic"
                : (model_choice.getValue() == model::equalities) ? "reduced"
                                                                 : "compact")
            << " branching=" << cmd.get_value_ordering() << std::endl;

  usrand(cmd.get_seed());


  StatisticList stats;
  stats.start();

#ifdef VERBOSE
  std::cout << "read\n";
#endif

  Instance jsp(cmd.get_filename().c_str(), format.getValue().c_str());


#ifdef VERBOSE
  std::cout << "ub\n";
#endif

  auto ub{jsp.getMakespanUpperBound(1)};

  if(init_ub.getValue() >= 0)
    ub = init_ub.getValue();


  Solver solver;
  cmd.set_parameters(solver);

  VarArray search_vars;

  // VarArray start_time(jsp.nTasks(), 0, ub);
  // VarArray end_time(jsp.nTasks(), 0, ub);
  VarArray start_time;
  VarArray end_time;

  Variable makespan(0,ub);
  Variable origin(0, 0);

  // if(model.getValue() == 0)
  //   model(jsp, solver, start_time, end_time, makespan, search_vars);
  // else if(model.getValue() == 0)
  build_model(
      jsp, solver, start_time, end_time, origin, makespan, search_vars, ub,
      (hall_pruning.getValue() ? model::basic : model_choice.getValue()),
      hall_pruning.getValue());
  // if(order_branching.getValue())
  // else
  // model_order(jsp, solver, start_time, end_time, makespan, search_vars);

  JacksonPreemptiveScheduler JPS(jsp, solver, start_time, end_time);

  // heuristic = solver.heuristic_factory(cmd.get_variable_ordering(),
  // cmd.get_value_ordering(), cmd.get_randomization());
  RestartPolicy *restart = restart_factory(cmd.get_restart_policy(), solver.parameters.restart_base, solver.parameters.restart_factor);

  BranchingHeuristic *heuristic =
      heuristic_factory(solver, cmd.get_value_ordering());

  solver.initialise_search(search_vars, heuristic, restart);

  // auto lb{get_lower_bound(jsp, solver, start_time, end_time, ub)};
  auto lb{JPS.get_lower_bound(ub)};

  solver.add(makespan >= lb);

  std::cout << " c initial bounds: [" << lb << ".." << ub << "]\n";

  std::cout << solver << std::endl;

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
      JPS.print_machine();
      exit(1);
    }
  };

  std::cout << "c solution checked!\n";

  if (cmd.printsolArg->getValue()) {
    print_solution(jsp, solver, start_time, end_time);
  }
}
