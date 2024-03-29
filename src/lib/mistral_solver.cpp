/*
Mistral 2.0 is a constraint satisfaction and optimisation library
Copyright (C) 2009  Emmanuel Hebrard

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

The author can be contacted electronically at emmanuel.hebrard@gmail.com.
*/

#include <sstream>
#include <fstream>
#include <signal.h>
#include <assert.h>


#include <mistral_sat.hpp>
//#include <mistral_search.hpp>
#include <mistral_solver.hpp>
#include <mistral_variable.hpp>
#include <mistral_constraint.hpp>


//#define _OLD_ true
//#define _DEBUG_NOGOOD true //(statistics.num_filterings == 491)

// #define _DEBUG_SEARCH true
// #define _DEBUG_AC true //(statistics.num_propagations > 1000)
//(!decisions.empty() && decisions.back().var.id() == 137)

//((statistics.num_filterings == 48212) || (statistics.num_filterings == 46738) || (statistics.num_filterings == 44368) || (statistics.num_filterings == 43659))

//#define _DEBUG_RESTORE true
//#define _DEBUG_REWRITE true
//#define _OUTPUT_TIKZ false




Mistral::Solver* active_solver;
static void Mistral_SIGINT_handler(int signum) {
  std::cout << std::endl
            << " " << active_solver->parameters.prefix_comment
            << " *************************************** INTERRUPTED "
               "***************************************"
            << std::endl;
  std::cout << active_solver->statistics << std::endl
            << " " << active_solver->parameters.prefix_comment
            << " *************************************** INTERRUPTED "
               "***************************************"
            << std::endl;
  exit(0);
}





// #ifdef _DEBUG_PROPAG 
// std::ostringstream *o_propag = NULL;
// #endif

Mistral::Solution::Solution( Vector< Variable >& vars ) {
  min_id = INFTY;
  max_id = -INFTY;
  int aux;
  for (unsigned int i = 0; i < vars.size; ++i) {
    variables.add(vars[i].get_var());
    aux = variables[i].id();
    if (aux < min_id)
      min_id = aux;
    if (aux > max_id)
      max_id = aux;
  }

  values = new int[max_id - min_id + 1];
  values -= min_id;

  for (unsigned int i = 0; i < variables.size; ++i) {
    values[variables[i].id()] = variables[i].get_solution_int_value();
  }
}
Mistral::Solution::~Solution() {

#ifdef _DEBUG_MEMORY
  std::cout << "c delete solution" << std::endl;
#endif

  values += min_id;
  delete[] values;
}

std::ostream& Mistral::Solution::display(std::ostream& os) const {
  if (!variables.empty()) {
    os // << variables[0] << ":"
        << values[variables[0].id()];
    for (unsigned int i = 1; i < variables.size; ++i) {
      os << " " // << variables[i] << ":"
         << values[variables[i].id()];
    }
  }
  return os;
}

int &Mistral::Solution::operator[](Variable x) { return values[x.id()]; }

Mistral::SolverParameters::SolverParameters() { initialise(); }
Mistral::SolverParameters::~SolverParameters() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete solver parameters" << std::endl;
#endif
}
void Mistral::SolverParameters::initialise() {
  find_all = 0;
  node_limit = 0;
  backtrack_limit = 0;
  propagation_limit = 0;
  fail_limit = 0;
  // restart_limit = 0;
  limit = 0;
  seed = 11041979;

  verbosity = 0;
  time_limit = -1;
  seed = 11041979;
  restart_policy = GEOMETRIC;
  restart_base = 256;
  restart_limit = 256;
  restart_factor = 1.33;
  activity_increment = 1e-2;
  normalize_activity = 0;
  init_activity = 1;
  forgetfulness = .75;
  randomization = 1; // 2;
  shuffle = false;   // true;
  activity_decay = 1.0;
  checked = 1;
  backjump = 0;
  value_selection = 2;
  dynamic_value = 0; // 1;

  prefix_comment = "c";
  prefix_statistics = "d";
  prefix_objective = "o";
  prefix_solution = "v";
  prefix_outcome = "s";
}
Mistral::SolverParameters::SolverParameters(const SolverParameters& sp) {
  copy(sp);
}
void Mistral::SolverParameters::copy(const SolverParameters& sp) {
  propagation_limit = sp.propagation_limit;
  fail_limit = sp.fail_limit;
  seed = sp.seed;

  time_limit = sp.time_limit;
  restart_policy = sp.restart_policy;
  restart_base = sp.restart_base;
  restart_limit = sp.restart_limit;
  restart_factor = sp.restart_factor;
  activity_increment = sp.activity_increment;
  normalize_activity = sp.normalize_activity;
  init_activity = sp.init_activity;
  forgetfulness = sp.forgetfulness;
  randomization = sp.randomization;
  shuffle = sp.shuffle;
  activity_decay = sp.activity_decay;
  checked = sp.checked;
  backjump = sp.backjump;
  value_selection = sp.value_selection;
  dynamic_value = sp.dynamic_value;

  verbosity = sp.verbosity;
  find_all = sp.find_all;
  node_limit = sp.node_limit;
  backtrack_limit = sp.backtrack_limit;
  fail_limit = sp.fail_limit;
  restart_limit = sp.restart_limit;
  limit = sp.limit;
  time_limit = sp.time_limit;
}

Mistral::SolverStatistics::SolverStatistics(Solver *s) {

  // VARNAME = {"virtual", "constant", "boolean", "range", "bitset", "list"};
  // METHOD_NAME = {
  //   "get_size"         ,
  //   "get_degree"       ,
  //   "get_min"          ,
  //   "get_max"          ,
  //   "get_initial_min"  ,
  //   "get_initial_max"  ,
  //   "get_min_pos"      ,
  //   "get_max_neg"      ,
  //   "next"             ,
  //   "prev"             ,
  //   "is_range"         ,
  //   "is_ground"        ,
  //   "equal"            ,
  //   "contain"          ,
  //   "intersect_range"  ,
  //   "included_range"   ,
  //   "includes_range"   ,
  //   "intersect_set"    ,
  //   "included_set"     ,
  //   "includes_set"     ,
  //   "intersect_to"     ,
  //   "union_to"         ,
  //   "remove"           ,
  //   "set_domain"       ,
  //   "set_min"          ,
  //   "set_max"          ,
  //   "set_domain_bitset",
  //   "remove_set"       ,
  //   "remove_interval"  ,
  //   "restore"
  // };

  initialise(s);

#ifdef _PROFILING
  init_prof();
#endif
}
Mistral::SolverStatistics::~SolverStatistics() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete solver statistics" << std::endl;
#endif
}
void Mistral::SolverStatistics::initialise(Solver *s) {
  solver = s;

  objective_value = 0;
  max_depth = 0;
  num_variables = 0;
  num_values = 0;
  num_constraints = 0;
  num_clauses = 0;
  num_learned = 0;
  size_learned = 0;
  num_nodes = 0;
  num_decisions = 0;
  num_restarts = 0;
  num_backtracks = 0;
  num_failures = 0;
  num_propagations = 0;
  num_solutions = 0;
  num_filterings = 0;
  // start_time = 0.0;
  creation_time = get_run_time();
  end_time = -1.0;
  start_time = end_time;
  best_time = 0.0;

  outcome = UNKNOWN;

  avg_learned_size = 0;
  // base_avg_size = 0;
  // learnt_avg_size = 0;
  // literals = 0;
  // small = 0;

  avg_amsc_expl_size = 0;
  ///
  num_amsc_explanations = 0;

  num_branch_on_large_domains = 0;

  // negative_weight = false;
  max_arity = 0;
}

#ifdef _PROFILING

void Mistral::SolverStatistics::init_prof() {
  total_propag_time = 0;
  total_branching_time = 0;
  total_restore_time = 0;

#ifdef _PROFILING_PRIMITIVE

  for (int i = 0; i < NUM_METHODS; ++i) {
    for (int j = 0; j < NUM_VARTYPES; ++j) {
      prof_num[i][j] = 0;
      prof_time[i][j] = 0;
    }
  }

#endif
}


int **comparison_array;
int compar_prof(const void *a, const void *b)
{
  int xa = (*(int *)a) / NUM_VARTYPES;
  int ya = (*(int *)a) % NUM_VARTYPES;

  int xb = (*(int *)b) / NUM_VARTYPES;
  int yb = (*(int *)b) % NUM_VARTYPES;

  if (comparison_array[xa][ya] == comparison_array[xb][yb])
    return 0;
  else if (comparison_array[xa][ya] < comparison_array[xb][yb])
    return 1;
  else
    return -1;
}

#ifdef _PROFILING_PRIMITIVE

std::ostream& Mistral::SolverStatistics::print_profile(std::ostream& os) const {
  int N = NUM_METHODS * NUM_VARTYPES;
  int order[N];
  unsigned long long int total_calls = 0;
  double total_time = 0.0;

  int method;
  int vartype;

  // os << "\n" << N << std::endl;

  for (int i = 0; i < N; ++i) {

    // method = i/NUM_VARTYPES;
    // vartype = i%NUM_VARTYPES;

    // os << method << " " << vartype << std::endl;

    // os << VAR_NAME[vartype] << "." << METHOD_NAME[method] << "(): "
    //    << prof_num[method][vartype] << " calls ("
    //    << (double)(prof_num[method][vartype]*100) / (double)total_calls <<
    //    "%), "
    //    << prof_time[method][vartype] << "s ("
    //    << prof_time[method][vartype]*100 / (double)total_time << "%)" <<
    //    std::endl;

    total_calls += prof_num[i / NUM_VARTYPES][i % NUM_VARTYPES];
    total_time += prof_time[i / NUM_VARTYPES][i % NUM_VARTYPES];
    order[i] = i;
  }

  comparison_array = new int *[NUM_METHODS];
  for (int i = 0; i < NUM_METHODS; ++i) {
    comparison_array[i] = new int[NUM_VARTYPES];
    for (int j = 0; j < NUM_VARTYPES; ++j)
      comparison_array[i][j] = prof_num[i][j];
  }

  // comparison_array = &((&(prof_time[0]))[0])

  // comparison_array = &((&prof_time[0])[0]);
  qsort(order, N, sizeof(int), compar_prof);

  for (int i = 0; i < N; ++i) {
    method = order[i] / NUM_VARTYPES;
    vartype = order[i] % NUM_VARTYPES;

    if (prof_num[method][vartype] == 0)
      break;

    std::string vname = std::string(VAR_NAME[vartype]);
    std::string mname = std::string(METHOD_NAME[method]);

    //+"."+METHOD_NAME[method]+"()";

    os //<< VAR_NAME[vartype] << "." << METHOD_NAME[method] << "(): "
        << std::left << std::setw(24) << vname + "." + mname + "()"
        << std::right << std::setw(10) << prof_num[method][vartype]
        << " calls (" << std::setprecision(6) << std::right << std::setw(7)
        << (double)(prof_num[method][vartype] * 100) / (double)total_calls
        << "%), " << std::right << std::setw(6) << prof_time[method][vartype]
        << "s (" << std::setprecision(6) << std::right << std::setw(7)
        << prof_time[method][vartype] * 100 / (double)total_time << "%)"
        << std::endl;
  }

  for (int i = 0; i < NUM_METHODS; ++i) {
    delete[] comparison_array[i];
  }
  delete[] comparison_array;

  return os;
}

#endif

#endif

std::ostream& Mistral::SolverStatistics::print_full(std::ostream& os) const {
  os << " " << solver->parameters.prefix_comment << " +" << std::setw(89)
     << std::setfill('=')
     //"=============================================================================
     << "+" << std::endl
     << std::setfill(' ') << std::left << " "
     << solver->parameters.prefix_outcome
     << std::setw(46 - solver->parameters.prefix_outcome.size()) << "  ";

  switch (outcome) {
  case SAT:
    os << std::right << std::setw(45) << "SATISFIABLE";
    break;
  case OPT:
    os << std::right << std::setw(45) << "OPTIMAL";
    break;
  case UNSAT:
    os << std::right << std::setw(45) << "UNSATISFIABLE";
    break;
  case UNKNOWN:
    os << std::right << std::setw(45) << "UNKNOWN";
    break;
  case LIMITOUT:
    if (num_solutions > 0)
      os << std::right << std::setw(45) << "SUBOPTIMAL";
    else
      os << std::right << std::setw(45) << "LIMITOUT";
    // break;
  }

  // Solution sol(solver->variables);

  std::string ps = solver->parameters.prefix_statistics;
  int lps = ps.size();

  os << std::endl
     //<< std::left << std::setw(46) << " " <<
     // solver->parameters.prefix_solution << " " << sol << std::endl
     //<< std::left << std::setw(46) << " " <<
     // solver->parameters.prefix_statistics << "  OBJECTIVE"
     << std::left << " " << solver->parameters.prefix_solution << std::right
     << std::setw(90 - lps) << num_solutions << std::endl;

  if (solver->objective) {
    if (solver->objective->is_optimization()) {
      os << std::left << " " << solver->parameters.prefix_statistics
         << std::setw(44 - lps) << "  OBJECTIVE" << std::right << std::setw(46)
         << objective_value << std::endl;

      os << std::left << " " << solver->parameters.prefix_statistics
         << std::setw(44 - lps) << "  BESTTIME" << std::right << std::setw(46)
         << best_time << std::endl;

      if (outcome == OPT)
        os << std::left << " " << solver->parameters.prefix_statistics
           << std::setw(44 - lps) << "  OPT" << std::right << std::setw(46) << 1
           << std::endl;
      else
        os << std::left << " " << solver->parameters.prefix_statistics
           << std::setw(44 - lps) << "  OPT" << std::right << std::setw(46) << 0
           << std::endl;

    } else if (solver->objective->is_satisfaction()) {
      os << std::left << " " << solver->parameters.prefix_statistics
         << std::setw(44 - lps) << "  MAXDEPTH" << std::right << std::setw(46)
         << max_depth << std::endl;
    } else if (solver->objective->is_enumeration()) {
      os << std::left << " " << solver->parameters.prefix_statistics
         << std::setw(44 - lps) << "  NUMSOL" << std::right << std::setw(46)
         << num_solutions << std::endl;
    }
  }

  os << std::left << " " << solver->parameters.prefix_statistics
     << std::setw(44 - lps) << "  SUCCESS" << std::right << std::setw(46)
     << (outcome != LIMITOUT && outcome != UNKNOWN) << std::endl
     << std::left << " " << solver->parameters.prefix_statistics
     << std::setw(44 - lps) << "  RUNTIME" << std::right << std::setw(46)
     << (end_time - start_time) << std::endl
     << std::left << " " << solver->parameters.prefix_statistics
     << std::setw(44 - lps) << "  PREPROCTIME" << std::right << std::setw(46)
     << (start_time - creation_time) << std::endl
     << std::left << " " << solver->parameters.prefix_statistics
     << std::setw(44 - lps) << "  MEMORY" << std::right << std::setw(46)
     << (mem_used() / 1048576.0) << std::endl
     << std::left << " " << solver->parameters.prefix_statistics
     << std::setw(44 - lps) << "  NODES" << std::right << std::setw(46)
     << num_nodes << std::endl
     << std::left << " " << solver->parameters.prefix_statistics
     << std::setw(44 - lps) << "  RESTARTS" << std::right << std::setw(46)
     << num_restarts << std::endl
     << std::left << " " << solver->parameters.prefix_statistics
     << std::setw(44 - lps) << "  FAILURES" << std::right << std::setw(46)
     << num_failures << std::endl
     << std::left << " " << solver->parameters.prefix_statistics
     << std::setw(44 - lps) << "  BACKTRACKS" << std::right << std::setw(46)
     << num_backtracks << std::endl
     << std::left << " " << solver->parameters.prefix_statistics
     << std::setw(44 - lps) << "  PROPAGATIONS" << std::right << std::setw(46)
     << num_propagations << std::endl
     << std::left << " " << solver->parameters.prefix_statistics
     << std::setw(44 - lps) << "  VARIABLES" << std::right << std::setw(46)
     << num_variables << std::endl
     << std::left << " " << solver->parameters.prefix_statistics
     << std::setw(44 - lps) << "  CONSTRAINTS" << std::right << std::setw(46)
     << num_constraints + (solver->base ? solver->base->clauses.size : 0)
     << std::endl
     << std::left << " " << solver->parameters.prefix_statistics
     << std::setw(44 - lps) << "  ARITY" << std::right << std::setw(46)
     << max_arity << std::endl
     << std::left << " " << solver->parameters.prefix_statistics
     << std::setw(44 - lps) << "  WEAKDEC" << std::right << std::setw(46)
     << num_branch_on_large_domains
     << std::endl
     //<< std::left << " " << solver->parameters.prefix_statistics <<
     // std::setw(44-lps) << "  NOGOODSIZE"
     //<< std::right << std::setw(46) << avg_learned_size << std::endl
     //<< std::left << " " << solver->parameters.prefix_statistics <<
     // std::setw(44-lps) << "  LEARNEDSIZE"
     //<< std::right << std::setw(46) << (num_learned ? size_learned/num_learned
     //: 0) << std::endl
     //<< std::left << " " << solver->parameters.prefix_statistics <<
     // std::setw(44-lps) << "  AMSCEXPLSIZE"
     //<< std::right << std::setw(46) << (num_amsc_explanations ?
     // avg_amsc_expl_size/num_amsc_explanations : 0) << std::endl
     //<< std::left << " " << solver->parameters.prefix_statistics <<
     // std::setw(44-lps) << "  NEGWEIGHT"
     //<< std::right << std::setw(46) << negative_weight << std::endl
     << " " << solver->parameters.prefix_comment << " +" << std::setw(89)
     << std::setfill('=') << "+" << std::endl
     << std::setfill(' ');
  //<< " " << parameters.prefix_comment << "
  //+=============================================================================+"
  //<< std::endl;
  return os;
}
std::ostream& Mistral::SolverStatistics::print_short(std::ostream& os) const {
  os << " " << solver->parameters.prefix_comment << " |";

  if (solver->parameters.verbosity > 2) {
    os << std::right << std::setw(7) << num_variables << " ";
    os << std::right << std::setw(5) << num_values / num_variables << " ";
    os << std::right << std::setw(8) << num_constraints + num_clauses << " |";
    os << std::right << std::setw(9) << num_learned << " ";
    os << std::right << std::setw(4) << (int)avg_learned_size << " ";
    if (num_learned)
      os << std::right << std::setw(4)
         << (num_learned ? size_learned / num_learned : 0) << " |";
    else
      os << " n/a |";
  }
  os << std::right << std::setw(9) << num_nodes << " ";
  os << std::right << std::setw(11) << num_propagations << " ";
  os << std::right << std::setw(8) << std::setprecision(5);
  os << (get_run_time() - start_time) << " |";

  os << std::right;
  if (solver->objective) {
    if (solver->objective->is_optimization())
      if (num_solutions)
        os << std::setw(10) << objective_value;
      else
        os << std::setw(10) << "nill";
    else if (solver->objective->is_satisfaction())
      os << std::setw(10) << max_depth;
    else if (solver->objective->is_enumeration())
      os << std::setw(10) << num_solutions;
  }

  os << " | ";
  return os;
}
std::ostream& Mistral::SolverStatistics::display(std::ostream& os) const {
  if (end_time >= 0.0) {
    if (solver->parameters.verbosity > 2)
      print_full(os);
    else {
      std::cout << " c +===========================================+\n";
      print_short(os);
      std::cout << "\n c +===========================================+\n";
    }
  } else
    print_short(os);
  return os;
}
Mistral::SolverStatistics::SolverStatistics(const SolverStatistics& sp) {
  copy(sp);
}
void Mistral::SolverStatistics::copy(const SolverStatistics& sp) {
  objective_value = sp.objective_value;
  best_time = sp.best_time;
  num_variables = sp.num_variables;
  num_constraints = sp.num_constraints;
  num_values = sp.num_values;
  num_nodes = sp.num_nodes;
  num_decisions = sp.num_decisions;
  num_restarts = sp.num_restarts;
  num_backtracks = sp.num_backtracks;
  num_failures = sp.num_failures;
  num_propagations = sp.num_propagations;
  num_solutions = sp.num_solutions;
  num_filterings = sp.num_filterings;
  start_time = sp.start_time;
  end_time = sp.end_time;
}
void Mistral::SolverStatistics::update(const SolverStatistics& sp) {
  objective_value = sp.objective_value;
  best_time = sp.best_time;

  num_nodes += sp.num_nodes;
  num_decisions += sp.num_decisions;
  num_restarts += sp.num_restarts;
  num_backtracks += sp.num_backtracks;
  num_failures += sp.num_failures;
  num_propagations += sp.num_propagations;
  num_solutions += sp.num_solutions;
  num_filterings += sp.num_filterings;
  if (end_time < sp.end_time)
    end_time = sp.end_time;
}

Mistral::ConstraintQueue::ConstraintQueue()
{
  min_priority = 0;
  cardinality = -1;
  higher_priority = -1;
  triggers = NULL;
  // taboo_constraint = NULL;
  _set_.initialise();
}



void Mistral::ConstraintQueue::reset_higher_priority() {
  while (--higher_priority >= min_priority && triggers[higher_priority].empty())
    ;
}

Constraint Mistral::ConstraintQueue::select(Mistral::Vector<Constraint>& constraints) {
  int cons_id = triggers[higher_priority].pop();
  Constraint cons = constraints[cons_id];
  _set_.fast_remove(cons_id);
  if (triggers[higher_priority].empty())
    reset_higher_priority();
  // taboo_constraint = cons.freeze();
  return cons;
}
// inline void select(Constraint cons) {
//   //_set_.remove(cons->id);
//   taboo_constraint = cons.freeze();
// }

void Mistral::ConstraintQueue::clear() {
  while (higher_priority >= min_priority)
    triggers[higher_priority--].clear();
  _set_.clear();
  // taboo_constraint = NULL;
}

void Mistral::ConstraintQueue::declare(Constraint c, Solver *s) {

  // std::cout << "declare " << " << parameters.prefix_comment << " << "(" <<
  // c.priority() << ") to the constraint queue" << std::endl; std::cout <<
  // "was: [" << min_priority << ","
  //  	    << min_priority+cardinality-1 << "]" << std::endl;

  int cons_idx = c.id();
  int cons_priority = c.priority();

  int new_min_p = min_priority;
  int new_max_p = min_priority + cardinality - 1;

  if (cons_idx == 0)
    solver = s;

  if (cons_priority < new_min_p || cons_priority > new_max_p) {
    if (cardinality > 0) {
      if (cons_priority < new_min_p)
        new_min_p = cons_priority;
      if (cons_priority > new_max_p)
        new_max_p = cons_priority;
    } else {
      new_min_p = cons_priority;
      new_max_p = cons_priority;
    }
    //
    // std::cout << "now: [" << new_min_p << ","
    // 	      << new_max_p << "]" << std::endl;
    //

    int new_cardinality = (new_max_p - new_min_p + 1);
    Queue *aux_triggers = triggers;
    triggers = new Queue[new_cardinality];
    triggers -= new_min_p;
    for (int i = min_priority; i < min_priority + cardinality; ++i) {
      triggers[i] = aux_triggers[i];
      aux_triggers[i].cancel();
    }

    triggers[cons_priority].initialise(cons_idx, cons_idx + 7);

    aux_triggers += min_priority;
    delete[] aux_triggers;

    if (higher_priority < min_priority)
      higher_priority = new_min_p - 1;
    min_priority = new_min_p;
    cardinality = new_cardinality;

    //
    // std::cout << " ==> init " << triggers[cons_priority] << std::endl;

  } else {

    // std::cout << "no need to create a new trigger list" << std::endl;
    // std::cout << "extend " << triggers[cons_priority] << " with " << cons_idx
    // << std::endl;
    //
    if (!triggers[cons_priority].is_initialised()) {
      triggers[cons_priority].initialise(cons_idx, cons_idx + 7);
    } else {
      triggers[cons_priority].extend(cons_idx);
    }
  }

  if (_set_.table)
    _set_.extend(cons_idx);
  else
    _set_.initialise(cons_idx, cons_idx, BitSet::empt);
}

void Mistral::ConstraintQueue::initialise(Solver *s)
{
  solver = s;
  _set_.initialise(0, 2 * s->constraints.capacity, BitSet::empt);

  int min_p = INFTY;
  int max_p = -INFTY;
  for (unsigned int i = 0; i < s->constraints.size; ++i) {
    if (s->constraints[i].priority() < min_p)
      min_p = s->constraints[i].priority();
    if (s->constraints[i].priority() > max_p)
      max_p = s->constraints[i].priority();
  }
  initialise(min_p, max_p, s->constraints.size);
}

void Mistral::ConstraintQueue::initialise(const int min_p, const int max_p,
                                          const int size) {
  cardinality = max_p - min_p + 1;
  min_priority = min_p;
  triggers = new Queue[cardinality];
  for (int i = 0; i < cardinality; ++i) {
    triggers[i].initialise(size);
  }
  triggers -= min_priority;
}

Mistral::ConstraintQueue::~ConstraintQueue() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete constraint queue" << std::endl;
#endif
  triggers += min_priority;
  delete[] triggers;
}

void Mistral::ConstraintQueue::trigger(GlobalConstraint *cons)//;
{
  // #ifdef _DEBUG_AC
  //  if(_DEBUG_AC) {
  // std::cout << " initial trigger for " << cons << "(" << (cons->id) << ") "<<
  // empty() << std::endl;
  //  }
  // #endif

  // int priority = cons->priority, cons_id = cons->id, triggered=false;
  Event evt;
  Variable x;

  for (unsigned int i = 0; i < cons->scope.size; ++i) {
    x = cons->_scope[i];

    if (!x.is_void()) {

      evt = ( // cons->scope[i].domain_type != BOOL_VAR &&
          x.is_ground() ? VALUE_EVENT : (LB_EVENT | UB_EVENT));

      if (x.is_ground() || cons->is_triggered_on(i, EVENT_TYPE(evt))) {
        // std::cout << "trigger " << cons << " because " << event2str(evt) << "
        // on "<< x << std::endl;
        trigger(cons, i, evt);
      }
    }
  }
}




void Mistral::ConstraintQueue::trigger(BinaryConstraint *cons)//;
{
  add(cons);
  // int cons_id = cons->id;
  // if(!_set_.fast_contain(cons_id)) {
  //   _set_.fast_add(cons_id);
  //   triggers[2].add(cons_id);
  //   if(2 > higher_priority) higher_priority = 2;
  // }
}

void Mistral::ConstraintQueue::trigger(TernaryConstraint *cons)//;
{
  add(cons);
  // int cons_id = cons->id;
  // if(!_set_.fast_contain(cons_id)) {
  //   _set_.fast_add(cons_id);
  //   triggers[2].add(cons_id);
  //   if(2 > higher_priority) higher_priority = 2;
  // }
}

void Mistral::ConstraintQueue::add(ConstraintImplementation *cons)//;
{
  int cons_id = cons->id;
  if (!_set_.fast_contain(cons_id)) {
    _set_.fast_add(cons_id);
    triggers[2].add(cons_id);
    if (2 > higher_priority)
      higher_priority = 2;
  }
}

void Mistral::ConstraintQueue::add(Constraint cons)//;
{
  add(cons.propagator);
}


void Mistral::ConstraintQueue::trigger(GlobalConstraint *cons, const int var, const Event evt)//;
{

  // int var = c.index();
  // GlobalConstraint *cons = (GlobalConstraint*)c.propagator;

  // #ifdef _DEBUG_AC
  //  if(_DEBUG_AC) {
  //    std::cout << "  triggers " << cons << " after a "
  // 	     << event2str(evt)
  // 	     << " event on "
  // 	     << cons->scope[var] << " in " << cons->scope[var].get_domain()
  // 	     << std::endl;
  //  }
  // #endif

  if (cons != solver->taboo_constraint) {
    int priority = cons->priority, cons_id = cons->id;
    if (_set_.fast_contain(cons_id)) {

      // std::cout << "NOTIFY OTHER EVENT " << event2strc(evt) << std::endl;

      cons->notify_other_event(var, evt);

      // std::cout << " ==> " << event2strc(cons->event_type[var]) << std::endl;
      //  if(cons->events.contain(var)) {
      //  	cons->event_type[var] |= evt;
      //  } else {
      //  	cons->events.add(var);
      //  	cons->event_type[var] = evt;
      //  }
    } else {

      // std::cout << "\nadd " << cons_id << " to " << _set_  << " " <<
      // _set_.table << std::endl;

      _set_.fast_add(cons_id);

      // std::cout << "here " << _set_ << " " << priority << std::endl;

      if (priority > higher_priority)
        higher_priority = priority;
      triggers[priority].add(cons_id);

      // std::cout << "NOTIFY OTHER EVENT " << event2strc(evt)  << std::endl;

      cons->notify_first_event(var, evt);

      // std::cout << " ==> " << event2strc(cons->event_type[var]) << std::endl;
      //  cons->events.set_to(var);
      //  cons->event_type[var] = evt;
    }
  }

#ifdef _DEBUG_QUEUE

  else
    std::cout << cons << "(" << (cons->id) << ") is currently being processed "
              << std::endl;

#endif


#ifdef _DEBUG_QUEUE

  std::cout << "=> " << cons->events << std::endl;

#endif

}


std::ostream& Mistral::ConstraintQueue::display(std::ostream& os) {
  int elt = INFTY;
  int end;

  os
      //<< cardinality << " (" << min_priority << ".." << higher_priority << ")"
      << "[";
  for (int i = cardinality - 1; i >= 0; --i) {
    if (!triggers[i + min_priority].empty()) {
      if (elt != INFTY)
        os << " | ";

      elt = triggers[i + min_priority].first();
      end = triggers[i + min_priority]._head;

      os << solver->constraints[elt];

      while (triggers[i + min_priority].next[elt] != end) {
        elt = triggers[i + min_priority].next[elt];
        os << " " << solver->constraints[elt];
      }
    }
  }
  os << "]";

  // for(int i=0; i<cardinality; ++i) {
  //   if(!triggers[i+min_priority].empty()) {
  //     elt = triggers[i+min_priority].first();
  //     end = triggers[i+min_priority]._head;
  //     os << "P" << (i+min_priority) << ": "
  // 	//<< solver->constraints[elt];
  // 	 << "["<< solver->constraints[elt].id() << "]";

  //     if(!_set_.contain(solver->constraints[elt].id())) {
  // 	std::cout << "inconsistent constraint queue" <<std::endl;
  // 	exit(1);
  //     }

  //     while(triggers[i+min_priority].next[elt] != end) {
  // 	elt = triggers[i+min_priority].next[elt];
  // 	os << ", " //<< solver->constraints[elt];
  // 	   << " [" << solver->constraints[elt].id() << "]";

  // 	if(!_set_.contain(solver->constraints[elt].id())) {
  // 	  std::cout << "inconsistent constraint queue" <<std::endl;
  // 	  exit(1);
  // 	}
  //     }
  //     os << std::endl;
  //   }
  // }
  return os;
}

Mistral::Solver::Solver()
#ifdef _MONITOR
    : monitor_list(this)
#endif
{
  // statistics
  statistics.initialise(this);

  consolidate_manager = NULL;
  search_started = false;

  // search stuf
  heuristic = NULL;
  policy = NULL;
  objective = NULL; // new Goal(Goal::SATISFACTION);
  // backjump_policy = NULL;

  // variables & constraints
  domain_types.initialise(0, 128);
  variables.initialise(0, 128);
  removed_variables.initialise(0, 32);
  // declared_variables.initialise(0,128);
  assignment_level.initialise(0, 128);
  assignment_order.initialise(0, 128);
  assignment_rank.initialise(this, 0);
  assigned.initialise(0, 128);
  visited.initialise(0, 1023);
  // reason.initialise(0,128);
  // EXPL
  reason_for.initialise(0, 128);
  // reason_index.initialise(0,128);
  constraints.initialise(0, 256);
  // constraint_graph.initialise(128);
  posted_constraints.initialise(0, 255, 256, false);
  sequence.initialise(this);
  sequence.initialise(128);
  initialised_vars = 0;
  initialised_cons = 0;
  num_search_variables = 0;
  base = NULL;

  search_root = -2;

  // trail stuff
  level = -1;
  // saved_objs.initialise(0,4096);
  saved_vars.initialise(0, 4096);
  trail_.initialise(0, 4096);
  decisions.initialise(0, 4096);
  // con_trail_.initialise(0,512);
  // con_trail_.add(0);
  // con_trail_.add(-1);

  // lit_activity.initialise(0,8192);
  //  var_activity.initialise(0,4096);
  // lit_activity = NULL;
  // var_activity = NULL;

  active_variables.initialise(4096);

  // params
  parameters.initialise();

  heuristic =
      NULL; // new GenericHeuristic< GenericDVO< MinDomain >, MinValue >(this);
  policy = NULL; // new Geometric();

  if (!random_generator_is_ready())
    usrand(parameters.seed);

  wiped_idx = CONSISTENT;
  prev_wiped_idx = CONSISTENT;

#ifdef _PARALLEL
  // used only with parallelization to indicate an exist because a solution is
  // found by another thread
  solution_found_elsewhere = NULL;
#endif

  save();
}

void Mistral::Solver::parse_pbo(const char *filename) {
  unsigned int LARGENUMBER = 131072;
  std::ifstream infile(filename);
  char c = ' ';
  std::string word;
  // int N, M, l=0;
  Literal lit;
  Vector<Vector<Literal>> clauses;
  Vector<Literal> new_clause;
  Vector<int> weight;
  Vector<Variable> scope;

  new_clause.initialise(0, 10);

  Variable Goal;
  int obj_dir = 0;

  // skip comments
  infile >> c;
  while (c == '*') {
    infile.ignore(LARGENUMBER, '\n');
    infile >> c;
  }

  infile.unget();

  int aux, parse_objective;

  bool should_continue = infile.good();

  while (should_continue) {

    //    std::cout << "DEBUT DE LIGNE! ";

    weight.clear();
    scope.clear();
    new_clause.clear();

    bool all_pones = true;
    bool all_mones = true;

    parse_objective = 0;
    do {
      infile >> c;

      should_continue = infile.good();

      if (should_continue) {

        if (c == 'm') {

          // std::cout << " obj" ;

          infile >> word;

          if (word == "in:") {
            parse_objective = 1;
          } else if (word == "ax:") {
            parse_objective = 2;
          }

          infile.ignore(LARGENUMBER, ' ');
          infile >> c;
        }

        if (c == '+' || c == '-') {

          // std::cout << " trm" ;

          infile >> aux;
          weight.add((c == '+' ? aux : -aux));

          all_pones &= (weight.back() == 1);
          all_mones &= (weight.back() == -1);

          // std::cout << " " << weight.back() << std::endl;

          infile.ignore(LARGENUMBER, ' ');
          infile >> c;

          assert(c == 'x');

          infile >> aux;
          while (aux > (int)(variables.size)) {
            Variable x(0, 1);
            add(x);
          }

          --aux;
          scope.add(variables[aux]);
        } else {

          // std::cout << " end" ;

          infile.unget();

          if (!parse_objective) {

            // std::cout << " CON" ;

            infile >> word;
            infile >> aux;
            infile.ignore(LARGENUMBER, '\n');

            int bounds[2] = {-INFTY, INFTY};

            if (word == ">=") {
              bounds[0] = aux;
            } else if (word == ">") {
              bounds[0] = aux + 1;
            } else if (word == "<=") {
              bounds[1] = aux;
            } else if (word == "<") {
              bounds[1] = aux - 1;
            } else if (word == "=") {
              bounds[0] = bounds[1] = aux;
            } else {
              std::cout << "unknown connector: " << word << " exiting"
                        << std::endl;
              exit(1);
            }

            // is it a clause?
            bool is_clause = false;
            int max_weight = 0;

            // compute the minimum
            int reachable[2] = {0, 0}; // = 0, reachable[1] = 0;
            for (unsigned int i = 0; i < weight.size; ++i) {
              if (weight[i] >= 0)
                reachable[1] += weight[i];
              else
                reachable[0] += weight[i];
              int abs_weight = std::abs(weight[i]);
              if (max_weight < abs_weight)
                max_weight = abs_weight;
              // mean_weight = ((mean_weight*(double)i) + (double)(weight[i])) /
              // (double)(i+1);
            }

            // check that each variable can make it true
            if (bounds[1] == INFTY) {
              is_clause = true;
              for (unsigned int i = 0; is_clause && i < weight.size; ++i) {
                if (weight[i] >= 0) {
                  is_clause = (reachable[0] + weight[i] >= bounds[0]);
                  lit = scope[i].id() * 2 + 1;
                  new_clause.add(lit);
                } else {
                  is_clause = (reachable[0] - weight[i] >= bounds[0]);
                  lit = scope[i].id() * 2;
                  new_clause.add(lit);
                }
              }
            } else if (bounds[0] == -INFTY) {
              is_clause = true;
              for (unsigned int i = 0; is_clause && i < weight.size; ++i) {
                if (weight[i] >= 0) {
                  is_clause = (reachable[1] - weight[i] <= bounds[1]);
                  lit = scope[i].id() * 2;
                  new_clause.add(lit);
                } else {
                  is_clause = (reachable[1] + weight[i] <= bounds[1]);
                  lit = scope[i].id() * 2 + 1;
                  new_clause.add(lit);
                }
              }
            }

            if (is_clause) {
              clauses.add(new_clause);
            } else {

              if (all_mones) {
                int aux = -bounds[0];
                bounds[0] = -bounds[1];
                bounds[1] = aux;
                all_pones = true;
              }

              if (all_pones) {

                // std::cout << "cardinality constraint" << std::endl;

                add(BoolSum(scope, // weight,
                            bounds[0], bounds[1]));

                // int real_max = std::min((int)scope.size, bounds[1]);
                // int real_min = std::max(0, bounds[0]);
                // int span = (real_max - real_min);

                // double center = (double)(real_min + real_max) / 2;
                // double skew = center/(double)(scope.size);

                // double activity_increment = parameters.activity_increment /
                // (scope.size * (1 << span));

                // if(activity_increment > 0.0) {
                //   int i=scope.size;
                //   while(i--) {
                //     lit_activity[2*scope[i].id()] += (1-skew) *
                //     activity_increment; lit_activity[2*scope[i].id()+1] +=
                //     skew * activity_increment; var_activity[scope[i].id()] +=
                //     activity_increment;
                //   }
                // }

              } else {

                // int real_max = std::min(reachable[1], bounds[1]);
                // int real_min = std::max(reachable[0], bounds[0]);
                // double span = (real_max - real_min);

                // double center = (double)(real_min + real_max) / 2;
                // double skew = center/(double)(reachable[1] - reachable[0] +
                // 1);

                // int exp = (int)(scope.size * (span)/(double)(reachable[1] -
                // reachable[0] + 1));

                // double activity_increment = parameters.activity_increment /
                // (1 << exp); double ac_i; // = activity_increment;

                // if(activity_increment > 0.0) {
                //   int i=scope.size;
                //   while(i--) {
                //     ac_i = activity_increment * (double)(std::abs(weight[i]))
                //     / (double)max_weight;
                //     lit_activity[2*scope[i].id()+(weight[i]<0)] += (1-skew) *
                //     ac_i; lit_activity[2*scope[i].id()+(weight[i]>=0)] +=
                //     skew * ac_i; var_activity[scope[i].id()] += ac_i;
                //   }
                // }

                // std::cout << "linear constraint" << std::endl;

                add(BoolSum(scope, weight, bounds[0], bounds[1]));
              }
            }
          } else {

            // std::cout << " OBJ" << std::endl;

            // minimize(BoolSum(scope, weight));

            Goal = BoolSum(scope, weight);
            obj_dir = parse_objective;
            // objective = new Goal(Goal::MINIMIZATION,  );

            infile.ignore(LARGENUMBER, '\n');
            // std::cout << "OBJECTIVE!" << std::endl;
            // exit(1);
          }

          // std::cout << " END LINE" << std::endl;

          break;
        }
      }
    } while (should_continue);
  }

  // std::cout << clauses << std::endl;

  if (parameters.backjump || clauses.size) {
    base = new ConstraintClauseBase(variables);
    add(base);
  }

  for (unsigned int i = 0; i < clauses.size; ++i) {

    add(clauses[i]);
  }

  if (!Goal.is_void()) {
    if (obj_dir == 1) {
      minimize(Goal);
    } else {
      maximize(Goal);
    }
  }
}

void Mistral::Solver::parse_dimacs(const char* filename) {
  unsigned int LARGENUMBER = 131072;
  std::ifstream infile(filename);
  char c = ' ';
  std::string word;
  int N, M, l = 0;
  Literal lit;
  Vector<Literal> new_clause;

  // skip comments
  infile >> c;
  while (c != 'p') {
    infile.ignore(LARGENUMBER, '\n');
    infile >> c;
  }

  infile >> word;
  assert(word == "cnf");

  // get number of atoms and clauses
  infile >> N;
  infile >> M;

  for (int i = 0; i < N; ++i) {
    Variable x(0, 1);
    add(x);
  }

  new_clause.initialise(0, N);
  // ConstraintClauseBase *cbase = new ConstraintClauseBase(variables);
  // base->reason = reason.stack_;
  // add(cbase);

  for (int i = 0; i < M; ++i) {
    new_clause.clear();
    do {
      infile >> l;
      if (l != 0) {
        if (l > 0)
          lit = (l - 1) * 2 + 1;
        else
          lit = (l + 1) * -2;
        new_clause.add(lit);
        // if(parameters.init_activity == 1)
        //   base->activity[lit] += parameters.activity_increment;
      }
    } while (l && infile.good());
    // cbase->add( new_clause );
    add(new_clause);

    // std::cout << new_clause << std::endl;

    // if(params.checked) add_original_clause( new_clause );
  }

  // init_watchers();

  // if(params.normalize_activity != 0)
  // normalize_activity(params.normalize_activity);

  //  std::cout << base << std::endl;

  // cbase->reason = reason.stack_;
}


void Mistral::Solver::add(Vector< Literal >& clause) {
  if (!base) {
    base = new ConstraintClauseBase(variables);
    add(base);
    // base->reason = reason.stack_;
    // EXPL
    // base->reason_for = reason_for.stack_;
  } // else {
  //   while(base->scope.size < variables.sizes) {
  //     base->add()
  //   }
  // }

  // double activity_increment = parameters.activity_increment / (1 <<
  // clause.size); if(activity_increment > 0.0) {
  //   int i=clause.size;
  //   while(i--) {
  //     lit_activity[NOT(clause[i])] += activity_increment;
  //     //lit_activity[clause[i]] += activity_increment;
  //     var_activity[UNSIGNED(clause[i])] += activity_increment;
  //   }
  // }

  base->add(clause,
            (parameters.init_activity ? parameters.activity_increment : 0.0));
}

void Mistral::Solver::set_parameters(SolverParameters &p) { parameters = p; }

void Mistral::Solver::add(VarArray& x) {

#ifdef _DEBUG_BUILD
  std::cout << "add " << x << std::endl;
#endif

  for (unsigned int i = 0; i < x.size; ++i)
    x[i].initialise(this);
}

void Mistral::Solver::add(Variable x) { 

#ifdef _DEBUG_BUILD
  std::cout << "add " << x << std::endl;
#endif

  x.initialise(this);
}

void Mistral::Solver::remove(Variable x) {
  int idx = x.id(), i, j;
  domain_types[idx] |= REMOVED_VAR;
  for (Event trig = 0; trig < 3; ++trig) {
    for (i = constraint_graph[idx].on[trig].size; i--;) {
      j = constraint_graph[idx].on[trig][i].id();
      if (posted_constraints.contain(j))
        posted_constraints.remove(j);
    }
  }
}

int Mistral::Solver::declare(Variable x) {

#ifdef _DEBUG_BUILD
  std::cout << "declare the variable "; //<< x << std::endl;
#endif
//
// std::cout << "declare the variable " << x << " " << variables.size << " "<< variables.capacity << std::endl;
// std::cout << "constraint graph " << constraint_graph.size << " " << constraint_graph.capacity << std::endl;

if (x.domain_type > DYN_VAR)
  booleans.add(&x);

// add the variables to the set of vars
active_variables.declare(variables.size);

x.variable->initialise(this);
// x.variable->id = variables.size;
// x.variable->solver = this;
visited.extend(variables.size);
variables.add(x);

// declared_variables.add(x);
assignment_level.add(INFTY);
assignment_order.add(INFTY);
// reason.add(NULL);
// EXPL
reason_for.add(NULL);
// reason_index.add(-1);
domain_types.add(DYN_VAR | (x.is_range() ? RANGE_VAR : 0));

last_solution_lb.add(-INFTY);
last_solution_ub.add(INFTY);

ConstraintTriggerArray array;
bool extend_struct = (constraint_graph.capacity == constraint_graph.size);
constraint_graph.add(array);
if (extend_struct) {
  // constraint_graph.add(array);
  for (int i = constraint_graph.size - 1; i >= 0; --i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = constraint_graph[i].on[j].size - 1; k >= 0; --k) {
        constraint_graph[i].on[j][k].re_link_to(&constraint_graph[i].on[j]);
      }
    }
  }
}
constraint_graph.back().initialise(4);

// while(lit_activity.capacity < 2*variables.size)
//   lit_activity.extendStack();
// while(var_activity.capacity < variable.size)
//   var_activity.extendStack();

// lit_activity.add(0);
// lit_activity.add(0);
// var_activity.add(0);

// is_relevant.resize(variables.size, variables.back().is_ground());

notify_add_variable();

// std::cout << x.get_domain() << std::endl;

#ifdef _DEBUG_BUILD
std::cout << x << " in " << x.get_domain() << std::endl;
#endif

return variables.size - 1;
}

void Mistral::Solver::checkcg(const char *msg) {

  for (size_t i = 0; i < constraints.size; ++i) {
    if (constraints[i].propagator->index) {

      int nbactive = 0;
      for (int j = 0; j < constraints[i].arity(); ++j) {
        nbactive += (constraints[i].propagator->on[j] != NULL);
      }

      if (nbactive > 1) {
        for (int j = 0; j < constraints[i].arity(); ++j) {
          Variable x = constraints[i].get_scope()[j];
          if (constraints[i].propagator->on[j]) {
            // if(x.is_ground()) {
            // 	std::cout <<
            // }

            if (constraints[i]
                    .propagator->on[j]
                    ->stack_[constraints[i].propagator->index[j]]
                    .propagator != constraints[i].propagator) {
              std::cout << msg << " " << constraints[i].id() << " " << j << " "
                        << x << " in " << x.get_domain() << std::endl;

              std::cout << constraints[i].propagator->on[j] << std::endl;
              exit(1);
            }

            assert(constraints[i]
                       .propagator->on[j]
                       ->stack_[constraints[i].propagator->index[j]]
                       .propagator == constraints[i].propagator);
          } else {
            assert(x.is_ground());
          }
        }
      }
    }
  }
}

void Mistral::Solver::add(Constraint c) { 

#ifdef _DEBUG_BUILD
  std::cout << "add " << c << std::endl;
#endif

  // checkcg("start add");

  if (c.id() < 0) {

    // std::cout << "INIT" << std::endl;

    c.initialise(this);

    // checkcg("after init");

    // if(constraints.size>=680) {
    // 	std::cout << std::endl << 22 << c << std::endl;
    // 	std::cout << c.propagator->on.size << std::endl;
    // 	std::cout << c.propagator->on[14] << std::endl;
    // 	// exit(1);
    // }

    // get a name for the constraint and add it to the list
    c.set_id(constraints.size);
    constraints.add(c);

    // checkcg("after add");

    active_constraints.declare(c, this);

    // checkcg("after declare");

    notify_add_constraint(c);

    // checkcg("after notify");

    c.post(this);

    // checkcg("after post");

  } else {

    c.awaken();

    // checkcg("after awaken?");
  }

  // std::cout << "================\n" << active_constraints <<
  // "\n================" << std::endl;

  c.trigger();

  // checkcg("after trigger");

  // std::cout << "================\n" << active_constraints <<
  // "\n================" << std::endl;

  if (level <= 0 && !posted_constraints.safe_contain(c.id())) {
    posted_constraints.init_add(c.id());
  }

  // checkcg("end add");
}

Mistral::Outcome Mistral::Solver::solve() {
  BranchingHeuristic *heu =
      new GenericHeuristic<GenericDVO<MinDomainOverWeight, 1,
                                      // PruningCountManager
                                      FailureCountManager>,
                           RandomMinMax>(this);
  RestartPolicy *pol = new Geometric();
  return depth_first_search(variables, heu, pol, NULL);
  // return (search_outcome == SAT || search_outcome == OPT);
}


void Mistral::Solver::minimize(Variable X) {
  X.initialise(this, 1);
  objective = new Goal(Goal::MINIMIZATION, X.get_var());

  if (!initialised_vars) {
    consolidate_manager = new ConsolidateListener(this);
  }
  consolidate_manager->id_obj = X.id();
}

void Mistral::Solver::maximize(Variable X) {
  X.initialise(this, 1);
  objective = new Goal(Goal::MAXIMIZATION, X.get_var());

  if (!initialised_vars) {
    consolidate_manager = new ConsolidateListener(this);
  }
  consolidate_manager->id_obj = X.id();
}

Mistral::Outcome Mistral::Solver::search_minimize(Variable X) {
  BranchingHeuristic *heu =
      new GenericHeuristic<GenericDVO<MinDomainOverWeight, 1,
                                      // PruningCountManager
                                      FailureCountManager>,
                           RandomMinMax>(this);
  RestartPolicy *pol = new Geometric();
  Goal *goal = new Goal(Goal::MINIMIZATION, X.get_var());
  return depth_first_search(variables, heu, pol, goal);
  // return (search_outcome == OPT);
}


Mistral::Outcome Mistral::Solver::search_maximize(Variable X) {
  BranchingHeuristic *heu =
      new GenericHeuristic<GenericDVO<MinDomainOverWeight, 1,
                                      // PruningCountManager
                                      FailureCountManager>,
                           RandomMinMax>(this);
  RestartPolicy *pol = new Geometric();
  Goal *goal = new Goal(Goal::MAXIMIZATION, X.get_var());
  return depth_first_search(variables, heu, pol, goal);
  // return (search_outcome == OPT);
}

Mistral::Outcome Mistral::Solver::depth_first_search(bool _restore_) {

  assert(sequence.size > 0);
  assert(policy != NULL);
  assert(objective != NULL);
  assert(heuristic != NULL);

  statistics.start_time = get_run_time();

  search_started = true;

  return restart_search(0, _restore_);
}

Mistral::Outcome Mistral::Solver::depth_first_search(BranchingHeuristic *heu,
                                                     RestartPolicy *pol,
                                                     Goal *goal,
                                                     bool _restore_) {
  return depth_first_search(variables, heu, pol, goal, _restore_);
}

Mistral::Outcome Mistral::Solver::depth_first_search(Vector<Variable> &seq,
                                                     BranchingHeuristic *heu,
                                                     RestartPolicy *pol,
                                                     Goal *goal,
                                                     bool _restore_) {
  // std::cout << "\nINIT LEVEL = " << level << std::endl;
  initialise_search(seq, heu, pol, goal, false);

  // std::cout << "SEARCH\n" ;
  return depth_first_search(_restore_);
}

Mistral::Outcome
Mistral::Solver::sequence_search(Vector<Vector<Variable>> &sequences,
                                 Vector<BranchingHeuristic *> &heuristics,
                                 Vector<RestartPolicy *> &policies,
                                 Vector<Goal *> &goals) {
  statistics.start_time = get_run_time();

#ifdef _DEBUG_SEARCH
  std::cout << " " << parameters.prefix_comment
            << " start new sequence search (in " << sequences.size << " phases)"
            << std::endl;
#endif

  unsigned int phase = 0;
  Outcome satisfiability = UNKNOWN, phase_satisfiability = UNKNOWN;
  Vector<int> phase_level;

  VarStack<Variable, ReversibleNum<int>> *copy_sequences =
      new VarStack<Variable, ReversibleNum<int>>[sequences.size];

  for (unsigned int i = 0; i < sequences.size; ++i) {
    copy_sequences[i].initialise(this);
    copy_sequences[i].initialise(sequences[i].size);
    for (unsigned int j = 0; j < variables.size; ++j) {
      Variable x = variables[j];
      copy_sequences[i].declare(x);
    }
  }

  for (unsigned int i = 0; i < sequences.size; ++i) {
    copy_sequences[i].clear();
    for (unsigned int j = sequences[i].size; j;) {
      Variable x = sequences[i][--j].get_var();
      if (!x.is_ground() && !copy_sequences[i].safe_contain(x) &&
          !(domain_types[x.id()] & REMOVED_VAR))
        copy_sequences[i].add(x);
    }
  }

  phase_level.add(level);
  // int max_phase = -1;

  // repeat until
  while (satisfiability == UNKNOWN) {
    // phase_level.add(level);

#ifdef _DEBUG_SEARCH
    std::cout << " c";
    for (unsigned int i = 0; i < phase; ++i)
      std::cout << "    ";
    std::cout << " ss init phase " << phase << " " << level << " "
              << phase_level << std::endl;
    std::cout << " c";
    for (unsigned int i = 0; i < phase; ++i)
      std::cout << "    ";
    std::cout << " ss search on " << sequences[phase] << std::endl;
#endif

    // initialise with the parameters of the current phase
    // if(phase > max_phase) {
    initialise_search(copy_sequences[phase], heuristics[phase], policies[phase],
                      goals[phase], false);

    std::cout << heuristics[phase] << std::endl;

#ifdef _DEBUG_SEARCH
    std::cout << " c";
    for (unsigned int i = 0; i < phase; ++i)
      std::cout << "    ";
    std::cout << " ss search phase " << phase << ": " << std::endl;
#endif

    // search the subset of variables
    if (objective->has_function()) {
      objective->set_type(phase < sequences.size - 1 ? Goal::SATISFACTION
                                                     : Goal::OPTIMIZATION);
    }
    phase_satisfiability = chronological_dfs(phase_level.back());

#ifdef _DEBUG_SEARCH
    std::cout << " c";
    for (unsigned int i = 0; i < phase; ++i)
      std::cout << "    ";
    std::cout << " ss " << outcome2str(phase_satisfiability) << std::endl;
#endif

    if (phase_satisfiability == UNSAT || phase_satisfiability == OPT) {

      int lvl = phase_level.pop();
      // the current phase is not satisfiable
      if (phase_level.empty()) {
#ifdef _DEBUG_SEARCH
        std::cout << " c";
        for (unsigned int i = 0; i < phase; ++i)
          std::cout << "    ";
        std::cout << " ss UNSAT! " << std::endl;
#endif

        // we have exhausted the complete search, returns
        satisfiability = phase_satisfiability;
      } else {
        // go back to the previous phase
        // phase_level.pop();
        --phase;

#ifdef _DEBUG_SEARCH
        std::cout << " c";
        for (unsigned int i = 0; i < phase; ++i)
          std::cout << "    ";
        std::cout << " ss backtrack to phase " << phase << " (level " << lvl - 1
                  << ")" << std::endl;
// exit(1);
#endif

        std::cout << decisions << std::endl;

        // int lvl = phase_level.pop();
        restore(lvl);
        decisions.size = lvl;

        // 	    Mistral::Decision deduction = decisions[lvl-1];
        // 	    deduction.invert();

        // #ifdef _DEBUG_SEARCH
        //   if(_DEBUG_SEARCH) {
        //     std::cout << parameters.prefix_comment;
        //     for(unsigned int k=0; k<=decisions.size; ++k) std::cout << " ";
        //     std::cout << deduction << std::endl;
        //   }
        // #endif

        // 	    deduction.make();

        branch_right();
        // phase_satisfiability = UNKNOWN;
      }
    } else if (phase_satisfiability == SAT) {

      phase_level.add(level);
      // the current phase has a solution
      if (++phase == sequences.size) {

#ifdef _DEBUG_SEARCH
        std::cout << " c";
        for (unsigned int i = 0; i < phase; ++i)
          std::cout << "    ";
        std::cout << " ss SAT! " << std::endl;
#endif

        // we have found a complete solution
        satisfiability = SAT;
      }

#ifdef _DEBUG_SEARCH
      else {
        std::cout << " c";
        for (unsigned int i = 0; i < phase; ++i)
          std::cout << "    ";
        std::cout << " ss advance to phase " << phase << std::endl;
      }
#endif
    }
  }

  statistics.outcome = satisfiability;

  statistics.end_time = get_run_time();

  if (parameters.verbosity) {
    std::cout << statistics;
  }

  delete[] copy_sequences;

  return satisfiability;
}

double Mistral::Solver::get_current_target() {
  double target = 0;
  if (objective) {
    if (objective->is_optimization()) {
      target = objective->upper_bound - objective->lower_bound;
    } else if (objective->is_satisfaction()) {
      target = statistics.num_variables - statistics.max_depth;
    } else if (objective->is_enumeration()) {
      if (statistics.num_solutions)
        target = 1.0 / (double)(statistics.num_solutions);
      else
        target = 2.0;
    }
  }
  return target;
}


Mistral::Outcome Mistral::Solver::restart_search(const int root, const bool _restore_) { //const bool _restore_, const bool _exit_on_solution_) {

  // int initial_level = level;

  if (parameters.verbosity > 2) {
    std::cout << " " << parameters.prefix_comment << " +" << std::setw(89)
              << std::setfill('=') << "+" << std::endl
              << std::setfill(' ') << " " << parameters.prefix_comment
              << " |      INSTANCE STATS   |      LEARNING      |         "
                 "SEARCH STATS          |";
    if (objective) {
      if (objective->is_optimization())
        std::cout << " OBJECTIVE |";
      else if (objective->is_satisfaction())
        std::cout << " MAX DEPTH |";
      else if (objective->is_enumeration())
        std::cout << " #SOLUTION |";
    }
    std::cout << std::endl
              << " " << parameters.prefix_comment
              << " |   vars  vals     cons |  #learnt size kept |    nodes     "
                 "propags     time |           |"
              << std::endl;
  } else if (parameters.verbosity > 1) {
    std::cout << " " << parameters.prefix_comment << " +" << std::setw(44)
              << std::setfill('=') << "+" << std::endl
              << std::setfill(' ') << " " << parameters.prefix_comment
              << " |         SEARCH STATS          |";
    if (objective) {
      if (objective->is_optimization())
        std::cout << " OBJECTIVE |";
      else if (objective->is_satisfaction())
        std::cout << " MAX DEPTH |";
      else if (objective->is_enumeration())
        std::cout << " #SOLUTION |";
    }
    std::cout << std::endl
              << " " << parameters.prefix_comment
              << " |    nodes     propags     time |           |" << std::endl;
  }

  Outcome satisfiability = UNKNOWN;

  statistics.objective_value = objective->value();

  statistics.num_variables = sequence.size;

  double last_target = get_current_target();
  double target_i;
  double progress_i = 0;

  // double total_progress = 0;
  // unsigned int tprog = 0;

  // std::cout << "[" << std::right << std::setw(33) << "]";
  // std::cout.flush();

  while (satisfiability == UNKNOWN) {

    statistics.num_constraints = posted_constraints.size;
    if (base)
      statistics.num_clauses = base->clauses.size;
    if (base)
      statistics.num_learned = base->learnt.size;
    statistics.num_variables = sequence.size;
    statistics.num_values = 0;
    for (int i = 0; i < sequence.size; ++i)
      statistics.num_values += sequence[i].get_size();

    if (parameters.verbosity > 1) {
      statistics.print_short(std::cout);
      std::cout << " " << (int)(100 * progress_i) << "%" << std::endl;
    }

    // if(progress_i > 0) {
    //   total_progress += (1.0-total_progress)*progress_i;
    //   tprog = (4294977296.0*total_progress);
    //   unsigned int i, j=31;
    //   //for(i=0; i<33; ++i) std::cout << "\b";
    //   for(i=(1<<j); j && i <= tprog; i+=(1<<(--j))) ;

    //   std::cout << (32-j) << " "<< tprog << std::endl;

    //   i = tprog = 32-j;
    //   // while(i--) std::cout << "=";
    //   // std::cout << std::right << std::setw(33-tprog) << "]";
    //   // std::cout.flush();
    // }

    ++statistics.num_restarts;

    // std::cout << " c notify restart" << std::endl;
    notify_restart(progress_i);

    // std::cout << "seq: [";
    // for(int i=sequence.size; i<variables.size; ++i) {
    //   std::cout << (sequence[i].get_value() ? " +" : " -") << sequence[i];
    // }
    // std::cout << "]\n";

    satisfiability = //(parameters.backjump ?
                     // conflict_directed_backjump() :
        chronological_dfs(root); //_exit_on_solution_); //);

    // if(_exit_on_solution_ && objective)
    //   satisfiability = objective->notify_solution(this);

    if (satisfiability == LIMITOUT) {

      policy->reset(parameters.restart_limit);

      if (!limits_expired()) {
        satisfiability = UNKNOWN;
      }
    }

    if (_restore_ || satisfiability == UNKNOWN)
      restore(root);

    forget();

    target_i = get_current_target();
    progress_i = (last_target - target_i) / last_target;
    last_target = target_i;

    is_relevant.reset();
    assert(is_relevant.size() >= variables.size);
    for (auto i{0}; i < variables.size; ++i) {
      if (variables[i].is_ground())
        is_relevant.set(i);
    }
  }

  if (satisfiability == LIMITOUT)
    statistics.outcome = interrupted();
  else
    statistics.outcome = satisfiability;

  statistics.end_time = get_run_time();

  if (parameters.verbosity) {
    std::cout << statistics;
  }

  return satisfiability;
}

Mistral::Outcome Mistral::Solver::get_next_solution()  
{
  Outcome satisfiability = UNSAT;

  if (search_started) {
    if (decisions.size)
      branch_right();
    else {
      return satisfiability;
    }
  }

  search_started = true;

  statistics.num_variables = sequence.size;
  statistics.num_values = 0;
  for (int i = 0; i < sequence.size; ++i)
    statistics.num_values += sequence[i].get_size();

  // display(std::cout);
  satisfiability = chronological_dfs();

  if (parameters.verbosity) {
    statistics.print_short(std::cout);
  }

  return satisfiability;
}

void Mistral::Solver::BooleanMemoryManager::add(Variable *x) {
  if (size.back() < 1024) {
    x->bool_domain = slots.back() + size.back();
    ++size.back();
  } else {
    int *nslot = new int[1024];
    std::fill(nslot, nslot + 1024, 3);
    size.add(1);
    slots.add(nslot);
    x->bool_domain = nslot;
  }

  // std::cout << "zzz " << *x << ": " << x->domain_type << std::endl;
}

void Mistral::Solver::initialise_search(Vector<Variable> &seq,
                                        BranchingHeuristic *heu,
                                        RestartPolicy *pol, Goal *goal,
                                        bool interrupted) {

#ifdef _DEBUG_BUILD
  std::cout << "INIT SEARCH!" << std::endl;
#endif

  is_relevant.resize(variables.size, 0);

  // cout << is_relevant.size() << " / " << variables.size << endl;

  if (level < 0)
    save();

  active_solver = this;

  if (interrupted)
    signal(SIGINT, Mistral_SIGINT_handler);

  sequence.clear();
  // decisions.clear();
  for (unsigned int i = seq.size; i;) {
    Variable x = seq[--i].get_var();
    if (!x.is_ground())
      if (!sequence.contain(x) && !(domain_types[x.id()] & REMOVED_VAR))
        sequence.add(x);
    // if(x.is_ground()) sequence.remove(x);
  }
  num_search_variables = sequence.size;

  if (heu) { // delete heuristic
    ;
    heuristic = heu;
  } else if (!heuristic)
    heuristic = new GenericHeuristic<Lexicographic, MinValue>(this);
  if (pol) { // delete policy;
    policy = pol;
  } else if (!policy)
    policy = new NoRestart();
  if (goal) { // delete objective;
    objective = goal;
  } else if (!objective)
    objective = new Goal(Goal::SATISFACTION);

  heuristic->initialise(sequence);

  parameters.restart_limit = policy->base;
  parameters.limit = (policy->base > 0);

  statistics.num_constraints = posted_constraints.size;

  if (base)
    statistics.num_clauses = base->clauses.size;

  unsigned int arity;
  for (unsigned int i = 0; i < posted_constraints.size; ++i) {
    arity = constraints[posted_constraints[i]].arity();
    if (arity > statistics.max_arity)
      statistics.max_arity = arity;
  }
}

void Mistral::Solver::initialise_search(
    VarStack<Variable, ReversibleNum<int>> &seq, BranchingHeuristic *heu,
    RestartPolicy *pol, Goal *goal, bool interrupted) {

  // consolidate();

  if (level < 0)
    save();

  active_solver = this;

  if (interrupted)
    signal(SIGINT, Mistral_SIGINT_handler);

  sequence.point_to(seq);

  // std::cout << "init seq: " << sequence << std::endl;

  num_search_variables = sequence.size;

  if (heu) { // delete heuristic
    ;
    heuristic = heu;
  } else if (!heuristic)
    heuristic = new GenericHeuristic<Lexicographic, MinValue>(this);
  if (pol) { // delete policy;
    policy = pol;
  } else if (!policy)
    policy = new NoRestart();
  if (goal) { // delete objective;
    objective = goal;
  } else if (!objective)
    objective = new Goal(Goal::SATISFACTION);

  // std::cout << (int*)heu << " " << (int*)heuristic << std::endl;
  // std::cout << heuristic << std::endl;
  // heuristic->display(std::cout);
  // std::cout << std::endl << sequence << std::endl;

  heuristic->initialise(sequence);

  parameters.restart_limit = policy->base;
  parameters.limit = (policy->base > 0);

  // statistics.num_constraints = constraints.size;
  statistics.num_constraints = posted_constraints.size;

  if (base)
    statistics.num_clauses = base->clauses.size;

  unsigned int arity;
  for (unsigned int i = 0; i < posted_constraints.size; ++i) {
    arity = constraints[posted_constraints[i]].arity();
    if (arity > statistics.max_arity)
      statistics.max_arity = arity;
  }

  // if(parameters.verbosity)  std::cout << " " << parameters.prefix_comment <<
  // " +" << std::setw(90) << std::setfill('=')
  // 				      << "+" << std::endl << std::setfill(' ')
  // 				      << " " << parameters.prefix_comment << " |
  // INSTANCE STATS       |                    SEARCH STATS                 |
  // OBJECTIVE
  // |" << std::endl
  // 				      << " " << parameters.prefix_comment << " |   vars
  // | vals |   cons |    nodes | filterings | propagations | cpu time | |" <<
  // std::endl;
}

Mistral::Solver::~Solver() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete solver" << std::endl;
#endif

  delete heuristic;
  delete policy;
  delete objective;

  // std::cout << "c delete consolidate manager" << std::endl;

  delete consolidate_manager;

  // std::cout << "delete constraints" << std::endl;
  for (unsigned int i = 0; i < constraints.size; ++i) {

    // std::cout << " delete cons " << constraints[i] << std::endl;

    delete constraints[i].propagator;

    // std::cout << " ok delete\n" ;
  }

  // std::cout << "delete expressions" << std::endl;
  for (unsigned int i = expression_store.size; i;) {

    // Variable x(expression_store[i]);
    // std::cout << " delete exp " << expression_store[i-1] << std::endl;

    delete expression_store[--i];

    // std::cout << " ok delete\n" ;
  }

  // std::cout << "delete variables" << std::endl;
  for (unsigned int i = 0; i < variables.size; ++i) {

    // // if(i==139) {
    // 	std::cout << "free var " << i << std::endl;
    // 	std::cout << variables[i] << " in " << variables[i].get_domain() <<
    // std::endl;
    //
    //
    // 	std::cout << variables[i].bitset_domain->domain.values.neg_words << " "
    // 						<< variables[i].bitset_domain->domain.values.pos_words <<
    // "
    // "
    // 						<< (int*)(variables[i].bitset_domain->domain.values.table)
    // << std::endl;
    // // }

    //
    // std::cout << removed_variables[0] << std::endl;
    variables[i].free_object();

    // std::cout << variables[i].bitset_domain->domain.values.neg_words << " "
    // 					<< variables[i].bitset_domain->domain.values.pos_words <<
    // "
    // "
    // 					<< (int*)(variables[i].bitset_domain->domain.values.table)
    // << std::endl;
    //
    //     //std::cout << removed_variables[0] << std::endl;
    // std::cout << "freed\n";

    // int domain_type = variables[i].domain_type;
    // if     (domain_type ==  BITSET_VAR) delete variables[i].bitset_domain;
    // else if(domain_type ==    LIST_VAR) delete variables[i].list_domain;
    // else if(domain_type ==   RANGE_VAR) delete variables[i].range_domain;
    // else if(domain_type == VIRTUAL_VAR) delete variables[i].virtual_domain;
    // else if(domain_type ==  EXPRESSION) delete variables[i].expression;
    // else if(domain_type !=   CONST_VAR) delete variables[i].variable;

    // #ifdef _CONSTRAINT_LIST
    //     delete constraint_graph[i].propagator;
    // #endif
  }

  for (unsigned int i = 0; i < removed_variables.size; ++i) {
    // std::cout << removed_variables[i] << std::endl;

    removed_variables[i].free_object();

    // int domain_type = variables[i].domain_type;
    // if     (domain_type ==  BITSET_VAR) delete variables[i].bitset_domain;
    // else if(domain_type ==    LIST_VAR) delete variables[i].list_domain;
    // else if(domain_type ==   RANGE_VAR) delete variables[i].range_domain;
    // else if(domain_type == VIRTUAL_VAR) delete variables[i].virtual_domain;
    // else if(domain_type ==  EXPRESSION) delete variables[i].expression;
    // else if(domain_type !=   CONST_VAR) delete variables[i].variable;

    // #ifdef _CONSTRAINT_LIST
    //     delete constraint_graph[i].propagator;
    // #endif
  }
}





// void Mistral::Solver::trigger_event(const int var, 
// 				    const Mistral::Event evt) {

//   active_variables.add(var, evt);

// }

// void Mistral::Solver::save(const int idx) {
//   saved_vars.add(idx);
// }

// void Mistral::Solver::() {

// }

void Mistral::Solver::restore() {

  unsigned int previous_level;
  // Constraint c;

  previous_level = trail_.pop();

#ifdef _DEBUG_RESTORE
  std::cout << "Restore to level " << previous_level << std::endl;
#endif

  while (saved_cons.size > previous_level) {

#ifdef _DEBUG_RESTORE
    std::cout << "  (c) " << saved_cons.back() << " -> ";
#endif

    saved_cons.pop().restore();

#ifdef _DEBUG_RESTORE
    std::cout << "  (c) " << saved_cons.back(0) << std::endl;
#endif
  }

  previous_level = trail_.pop();
  while (saved_ints.size > previous_level) {

#ifdef _DEBUG_RESTORE
    std::cout << "  (i) " << saved_ints.back() << " -> ";
#endif

    saved_ints.pop()->restore();

#ifdef _DEBUG_RESTORE
    std::cout << "  (i) " << saved_ints.back(0) << std::endl;
#endif
  }

  previous_level = trail_.pop();
  while (saved_lists.size > previous_level) {

#ifdef _DEBUG_RESTORE
    std::cout << "  (l) " << *(saved_lists.back()) << " -> ";
#endif

    saved_lists.pop()->restore();

#ifdef _DEBUG_RESTORE
    std::cout << "  (l) " << *(saved_lists.back(0)) << std::endl;
#endif
  }

  previous_level = trail_.pop();
  while (saved_bools.size > previous_level) {

#ifdef _DEBUG_RESTORE
    std::cout << "  (b) " << saved_bools.back() << " -> ";
#endif

    *(saved_bools.pop()) = 3;

#ifdef _DEBUG_RESTORE
    std::cout << "  (b) " << saved_bools.back(0) << std::endl;
#endif
  }

  previous_level = trail_.pop();
  while (saved_vars.size > previous_level) {

#ifdef _DEBUG_RESTORE
    std::cout << "  (v) " << variables[saved_vars.back()] << " in "
              << variables[saved_vars.back()].get_domain() << " -> ";
#endif

    variables[saved_vars.pop()].restore();

#ifdef _DEBUG_RESTORE
    std::cout << "  (v) " << variables[saved_vars.back(0)] << " in "
              << variables[saved_vars.back(0)].get_domain() << std::endl;
#endif
  }

  --level;
  ++statistics.num_backtracks;

  // unsigned int previous_level;

  // previous_level = trail_.pop();
  // //std::cout << saved_relax.size << " -> " <<  previous_level << std::endl;

  // while( saved_relax.size > previous_level )
  //   saved_relax.pop()->post_on_array();

  // previous_level = trail_.pop();
  // //std::cout << saved_post.size << " -> " <<  previous_level << std::endl;

  // while( saved_post.size > previous_level )
  //   saved_post.pop()->relax_from_array();

  // previous_level = trail_.pop();
  // while( saved_vars.size > previous_level ) {
  //   // if(saved_vars.back() == 4) {
  //   //   std::cout << level << " RESTORE X4: " << variables[4] << " in " <<
  //   variables[4].get_domain()
  //   // 		<< " " << (variables[4].domain_type == RANGE_VAR ?
  //   // ((VariableRange*)(variables[4].variable))->trail_
  //   :
  //   // ((VariableBitmap*)(variables[4].variable))->trail_)
  //   // 		<< std::endl;
  //   // }

  //   variables[saved_vars.pop()].restore();

  //   // if(saved_vars.back(0) == 4) {
  //   //   std::cout << level << " ======> X4: " << variables[4] << " in " <<
  //   variables[4].get_domain()
  //   // 		<< " " << (variables[4].domain_type == RANGE_VAR ?
  //   // ((VariableRange*)(variables[4].variable))->trail_
  //   :
  //   // ((VariableBitmap*)(variables[4].variable))->trail_)
  //   // 		<< std::endl << std::endl;
  //   // }
  // }

  // previous_level = trail_.pop();
  // while( saved_objs.size > previous_level )
  //   saved_objs.pop()->restore();

  // if(con_trail_.back() == level) {
  //   con_trail_.pop();
  //   posted_constraints.size = con_trail_.pop();
  // }

  // //  decisions.pop();
  // sequence.size = trail_.pop();

  // ++statistics.num_backtracks;
  // --level;
}

void Mistral::Solver::restore(const int lvl) {
  if (lvl >= search_root && search_root >= -1) {
    // decisions.size = lvl;
    decisions.size = lvl - search_root;
    while (lvl < level)
      restore();
  }
}


void Mistral::Solver::add(Mistral::SolutionListener* l) {
  l->mid = solution_triggers.size;
  solution_triggers.add(l);
}
void Mistral::Solver::add(Mistral::RestartListener* l) {
  l->rid = restart_triggers.size;
  restart_triggers.add(l);
}
void Mistral::Solver::add(Mistral::SuccessListener* l) {
  l->sid = success_triggers.size;
  success_triggers.add(l);
}
void Mistral::Solver::add(Mistral::BacktrackListener* l) {
  l->fid = backtrack_triggers.size;
  backtrack_triggers.add(l);
}
// void Mistral::Solver::add(Mistral::FailureListener* l) {
//   l->fid = failure_triggers.size;
//   failure_triggers.add(l);
// }
void Mistral::Solver::add(Mistral::DecisionListener* l) {
  l->did = decision_triggers.size;
  decision_triggers.add(l);
}
void Mistral::Solver::add(Mistral::VariableListener* l) {
  l->vid = variable_triggers.size;
  variable_triggers.add(l);
}
void Mistral::Solver::add(Mistral::ConstraintListener* l) {
  l->cid = constraint_triggers.size;
  constraint_triggers.add(l);
}


void Mistral::Solver::remove(Mistral::SolutionListener* l) {
  unsigned int idx = l->mid;
  solution_triggers.remove(idx);
  if (solution_triggers.size > idx)
    solution_triggers[idx]->mid = idx;
}
void Mistral::Solver::remove(Mistral::RestartListener* l) {
  unsigned int idx = l->rid;
  restart_triggers.remove(idx);
  if (restart_triggers.size > idx)
    restart_triggers[idx]->rid = idx;
}
void Mistral::Solver::remove(Mistral::SuccessListener* l) {
  unsigned int idx = l->sid;
  success_triggers.remove(idx);
  if (success_triggers.size > idx)
    success_triggers[idx]->sid = idx;
}
void Mistral::Solver::remove(Mistral::BacktrackListener* l) {
  unsigned int idx = l->fid;
  backtrack_triggers.remove(idx);
  if (backtrack_triggers.size > idx)
    backtrack_triggers[idx]->fid = idx;
}
// void Mistral::Solver::remove(Mistral::FailureListener* l) {
//   unsigned int idx = l->fid;
//   failure_triggers.remove(idx);
//   if(failure_triggers.size>idx)
//     failure_triggers[idx]->fid = idx;
// }
void Mistral::Solver::remove(Mistral::DecisionListener* l) {
  unsigned int idx = l->did;
  decision_triggers.remove(idx);
  if (decision_triggers.size > idx)
    decision_triggers[idx]->did = idx;
}
void Mistral::Solver::remove(Mistral::VariableListener* l) {
  unsigned int idx = l->vid;
  variable_triggers.remove(idx);
  if (variable_triggers.size > idx)
    variable_triggers[idx]->vid = idx;
}
void Mistral::Solver::remove(Mistral::ConstraintListener* l) {
  unsigned int idx = l->cid;
  constraint_triggers.remove(idx);
  if (constraint_triggers.size > idx)
    constraint_triggers[idx]->cid = idx;
}


void Mistral::Solver::notify_backtrack() { //Constraint *con, const int idx) {
  for (unsigned int i = 0; i < backtrack_triggers.size; ++i) {
    backtrack_triggers[i]->notify_backtrack();
  }
} 

void Mistral::Solver::notify_success() { //Variable* changes, const int n) {
  for (unsigned int i = 0; i < success_triggers.size; ++i) {
    success_triggers[i]->notify_success();
  }
  if ((int)(statistics.max_depth) < assignment_rank) {
    // for(unsigned int i=0; i<variables.size; ++i)
    //   {
    // 	if(last_solution_lb[i] != -INFTY)
    // 	  std::cout << last_solution_lb[i];
    // 	else
    // 	  std::cout << ".";
    //   }
    // std::cout << std::endl;
    statistics.max_depth = assignment_rank;
  }
}
void Mistral::Solver::notify_decision() { //Decision d) {
  for (unsigned int i = 0; i < decision_triggers.size; ++i) {
    decision_triggers[i]->notify_decision();
  }
}

void Mistral::Solver::notify_restart(const double prog) {
  for (unsigned int i = 0; i < restart_triggers.size; ++i) {
    // std::cout << " c notify restart (2)" << std::endl;
    restart_triggers[i]->notify_restart(prog);
  }
}

void Mistral::Solver::notify_relax(Constraint c) {
  for (unsigned int i = 0; i < constraint_triggers.size; ++i) {
    constraint_triggers[i]->notify_relax(c);
  }
}

void Mistral::Solver::notify_add_constraint(Constraint c) {
  for (unsigned int i = 0; i < constraint_triggers.size; ++i) {
    constraint_triggers[i]->notify_add_con(c);
  }
}

void Mistral::Solver::notify_post(Constraint c) {
  for (unsigned int i = 0; i < constraint_triggers.size; ++i) {
    constraint_triggers[i]->notify_post(c);
  }
}

void Mistral::Solver::notify_change_variable(const int idx) {

  // std::cout << "notify change on " << variables[idx] << " in " <<
  // variables[idx].get_domain()
  // 	    << " to " <<  variable_triggers << std::endl;

  for (unsigned int i = 0; i < variable_triggers.size; ++i) {
    variable_triggers[i]->notify_change(idx);
  }
}

void Mistral::Solver::notify_add_variable() {
  for (unsigned int i = 0; i < variable_triggers.size; ++i) {
    variable_triggers[i]->notify_add_var();
  }
}

void Mistral::Solver::consolidate() 
{
  if (!initialised_vars) {
    consolidate_manager = new ConsolidateListener(this);
  }

  for (auto con : constraints) {
    con.propagator->mark_domain();
  }

  for (; initialised_vars < variables.size; ++initialised_vars) {
    variables[initialised_vars] = variables[initialised_vars].get_var();

    // std::cout << variables[initialised_vars].get_var() << " " <<
    // domain_types[initialised_vars] << " / " <<
    // variables[initialised_vars].domain_type << " (" << RANGE_VAR << ")\n";

    if (!(domain_types[initialised_vars] & RANGE_VAR) &&
        variables[initialised_vars].domain_type == RANGE_VAR &&
        !variables[initialised_vars].is_ground() &&
        variables[initialised_vars].get_degree() > 0) {

      // std::cout << "change " << domain_types[initialised_vars] << std::endl;

      int minval = variables[initialised_vars].get_min();
      int maxval = variables[initialised_vars].get_max();
      Variable X(minval, maxval, domain_types[initialised_vars]);

      X.variable->solver = this;
      X.variable->id = initialised_vars;
      variables[initialised_vars] = X;

      if (X.domain_type > DYN_VAR) {
        booleans.add(variables.stack_ + initialised_vars);
      }
    }
  }

  while (initialised_cons < constraints.size)
    constraints[initialised_cons++].consolidate();
}


void Mistral::Solver::make_non_convex(const int idx) 
{

  // std::cout << variables << std::endl;

  // if(idx==139)
  // std::cout << "\nmake " << variables[idx] << " non convex\n";

  if (variables[idx].domain_type == RANGE_VAR) {
    Variable r = variables[idx];

    Variable X(r, true);

    // if(variables[idx].is_expression()) {

    //   std::cout << "(exp)" << std::endl;

    //   Variable x = variables[idx];
    //   while (x.expression->self.is_expression()) {
    // 	x = x.expression->self;
    //   }
    //   x.expression->self = X;

    // if(idx==17) {
    // 	int ids = sequence.index(idx);
    // 	std::cout << "change " << variables[idx] << " <- " << X << " (" <<
    // sequence.list_[ids] << ")" << std::endl;
    // }

    // } else {
    variables[idx] = X; //.get_var();
    //
    // std::cout << X.bitset_domain->domain.values.neg_words << " "
    // 					<< X.bitset_domain->domain.values.pos_words <<
    // "
    // "
    // 					<< (int*)(X.bitset_domain->domain.values.table)
    // << std::endl;
    //

    // if(ids>=0) sequence.list_[ids] = X;

    // for(int i=0; i<3; ++i) {
    //   for(int j=constraint_graph[idx].on[i].size-1; j>=0; --j)
    // 	constraint_graph[idx].on[i][j].consolidate();
    // }

    // std::cout << "MADE " << variables[idx] << " IN " <<
    // variables[idx].get_domain() << " NON CONVEX: NOTIFY CHANGE" << std::endl;

    notify_change_variable(idx);
    // }

    // if(idx==17) {
    // 	int ids = sequence.index(idx);
    // 	std::cout << "res: " << variables[idx] << " (" << sequence.list_[ids] <<
    // ")" << std::endl;
    // }

    // std::cout << parameters.prefix_comment << std::endl;
    // unsigned int k=0;
    // // for(; k<monitored_variables.size; ++k) {
    // //   monitored_variables[k] = monitored_variables[k].get_var();
    // // }
    // // k = 0;
    // for(unsigned int q=0; q<monitored_index.size; ++q) {
    //   std::cout << "c "; for(int lvl=0; lvl<level; ++lvl) std::cout << " ";
    //   for(; k<monitored_index[q]; ++k) {
    // 	std::cout << "|" << variables[monitored[k]] ;
    // 	std::cout.flush();
    // 	std::cout << variables[monitored[k]].get_domain();
    // 	std::cout.flush();
    // 	std::cout << variables[monitored[k]].get_history() << " ";
    //   }
    //   std::cout << std::endl;
    // }

    // r.free_object();

    removed_variables.add(r);
  }

  // std::cout << "done\n";
}

void Mistral::Solver::parity_processing(const int Knst) {
  std::vector<BitSet *> scopes;
  std::vector<int> sizes;
  std::vector<int> parities;
  // Get all parity scopes

  Variable *scope;
  int arity, min_idx, max_idx;
  for (unsigned int i = 0; i < constraints.size; ++i) {
    if (constraints[i].symbol() == "parity") {
      scope = constraints[i].get_scope();
      arity = constraints[i].arity();
      min_idx = max_idx = scope[0].id();
      for (int j = arity; --j;) {
        if (scope[j].id() < min_idx)
          min_idx = scope[j].id();
        else if (scope[j].id() > max_idx)
          max_idx = scope[j].id();
      }
      BitSet *e = new BitSet(min_idx, max_idx, BitSet::empt);
      for (int j = arity; j--;) {
        e->add(scope[j].id());
      }
      // std::cout << e << std::endl;
      scopes.push_back(e);
      sizes.push_back(arity);
      parities.push_back(
          ((ConstraintParity *)(constraints[i].propagator))->target_parity);
    }
  }

  BitSet f1(0, variables.size, BitSet::empt);
  BitSet f2(0, variables.size, BitSet::empt);
  int a1, a2;

  int original_size = scopes.size();

  // for(std::vector< BitSet* >::iterator it=scopes.begin(); it<scopes.end();
  // ++it) {
  //   std::cout << (*it) << std::endl;
  // }

  // // int i=0;
  // // for(std::vector< BitSet >::iterator it=scopes.begin();
  // it<(scopes.end()-1); ++it, ++i) {
  // //   //std::cout << (*it) << std::endl;
  // //   int j=i+1;
  // //   for(std::vector< BitSet >::iterator jt=(it+1); jt!=scopes.end(); ++jt,
  // ++j) {

  for (unsigned int i = 0; i < scopes.size() - 1; ++i) {
    for (unsigned int j = i + 1; j < scopes.size(); ++j) {

      // std::cout << i << " x " << j << ": " << scopes[i] << " x " << scopes[j]
      // << std::endl;

      scopes[i]->union_to(f1);
      scopes[j]->setminus_to(f1);
      scopes[j]->union_to(f2);
      scopes[i]->setminus_to(f2);

      // f1 = ei \ ej;
      // f2 = ej \ ei;

      a1 = f1.size();
      a2 = f2.size();

      int cooling = (scopes.size() / (Knst * original_size)) - 1;

      if (!a1) {

        // std::cout << a2 << std::endl;

        if (!a2) {
          if (parities[i] != parities[j])
            fail();
        } else {

          BitSet *ne = new BitSet(f2);

          std::cout << scopes[j] << " \\ " << scopes[i] << " => " << ne
                    << std::endl;

          // std::cout << scopes.size() << " + "<< ne << " - " << scopes[j] <<
          // std::endl;

          // ei is included in ej, hence we can remove ej and replace scopes[i]
          // by f2
          scopes.push_back(ne);
          sizes.push_back(a2);
          parities.push_back(parities[i] != parities[j]);
          // std::cout << "reduction!\n";
        }
        // exit(1);

      } else if (!a2) {

        // std::cout << a1 << std::endl;

        if (!a1) {
          if (parities[i] != parities[j])
            fail();
        } else {

          BitSet *ne = new BitSet(f1);

          std::cout << scopes[i] << " \\ " << scopes[j] << " => " << ne
                    << std::endl;

          // std::cout << scopes.size() << "  + "<< ne << " - " << scopes[i] <<
          // std::endl;

          // ej is included in ei, hence we can remove ei and replace it by f1
          scopes.push_back(ne);
          sizes.push_back(a1);
          parities.push_back(parities[i] != parities[j]);
          // std::cout << "reduction!\n";
          // exit(1);
        }

      } else if ((a1 + a2 < sizes[i] - cooling) &&
                 (a1 + a2 < sizes[j] - cooling)) {
        // std::cout << scopes.size() << "  + "<< f1 << " u " << f2 <<
        // std::endl;

        BitSet *ne = new BitSet(0, variables.size, BitSet::empt);
        ne->union_with(f1);
        ne->union_with(f2);

        std::cout << scopes[i] << " <> " << scopes[j] << " => " << ne
                  << std::endl;

        bool already_in = false;
        for (unsigned int k = 0; !already_in && k < scopes.size(); ++k) {
          already_in |= (*(scopes[i]) == *ne);
          // std::cout << "check " << scopes[i] << " == " << ne << " " <<
          // (already_in ? "yes" : "no") << std::endl;
        }

        if (!already_in) {
          scopes.push_back(ne);
          sizes.push_back(a1 + a2);
          parities.push_back(parities[i] != parities[j]);
        }
      } //  else {
      // 	std::cout << " nothing\n";
      // }

      f1.clear();
      f2.clear();
    }
  }

  // for(std::vector< BitSet* >::iterator it=scopes.begin(); it<scopes.end();
  // ++it) {
  //   std::cout << std::setw(3) << sizes[(int)(it-scopes.begin())] << (*it) <<
  //   std::endl;
  // }

  VarArray pscope;
  for (unsigned int i = original_size; i < scopes.size(); ++i) {
    int nxt = scopes[i]->min(), j;
    do {
      j = nxt;
      pscope.add(variables[j]);
      nxt = scopes[i]->next(j);
    } while (nxt > j);

    add(Parity(pscope, parities[i]));
    pscope.clear();
  }
}

bool Mistral::Solver::is_pseudo_boolean() const {

  bool is_pb = true;

  for (unsigned int i = 0; is_pb && i < variables.size; ++i) {
    is_pb &= ((domain_types[i] & REMOVED_VAR) || variables[i].is_boolean() ||
              ((objective != NULL) &&
               (objective->objective.operator_equal(variables[i]))));

    if (!is_pb) {
      std::cout << " " << parameters.prefix_comment
                << " not pseudo boolean because of " << variables[i] << " in "
                << variables[i].get_domain() << std::endl;
    }
  }

  for (unsigned int i = 0; is_pb && i < posted_constraints.size; ++i) {
    is_pb &= constraints[posted_constraints[i]].explained();

    if (!is_pb) {
      std::cout << " " << parameters.prefix_comment
                << " not pseudo boolean because "
                << constraints[posted_constraints[i]] << " is not explained"
                << std::endl;
    }
  }

  if (is_pb) {
    std::cout << " " << parameters.prefix_comment
              << " the problem is pseudo boolean!" << std::endl;
  }

  return is_pb;
}

bool Mistral::Solver::simple_rewrite()
{

  int i = 0;
  unsigned int con_i = 0;
  IntStack to_rewrite(0, 2 * constraints.size, false);
  Constraint con;
  RewritingOutcome rewritten;

#ifdef _DEBUG_REWRITE
  std::cout << "\n===================START REWRITING=====================\n";
#endif

  do {

#ifdef _DEBUG_REWRITE
    std::cout << " propagate  ";
#endif

    if (!propagate())
      break;

#ifdef _DEBUG_REWRITE
    std::cout << (*this) << std::endl;
// std::cout << "Collect rewritable constraints: " << std::endl;
#endif

    while (con_i < constraints.size) {
      // add the unseen constraints that might be rewritten
      for (; con_i < constraints.size; ++con_i) {

        // #ifdef _DEBUG_REWRITE
        //       std::cout << "   " << constraints[con_i] ;
        // #endif

        if (constraints[con_i].simple_rewritable()) {
          to_rewrite.safe_add(con_i);

          // #ifdef _DEBUG_REWRITE
          // 	std::cout << " in" << std::endl;
          // #endif
        }

        // #ifdef _DEBUG_REWRITE
        //       else {
        // 	std::cout << " out" << std::endl;
        //       }
        // #endif
      }

      // // add all rewritable constraints to the stack
      // for(i=0; i<to_rewrite.size; ++i)
      //   active_constraints.add(constraints[to_rewrite[i]]);

#ifdef _DEBUG_REWRITE
      std::cout << " rewrite: " << to_rewrite << std::endl;
#endif

      // // rewriting
      // while(!active_constraints.empty()) {

      //   con = active_constraints.select(constraints);
      for (i = to_rewrite.size; i--;) {
        con = constraints[to_rewrite[i]];

#ifdef _DEBUG_REWRITE
        std::cout << "   [rewrite c" << con.id() << " " << con << std::endl;
#endif

        rewritten = con.rewrite();

        switch (rewritten) {
        case NO_EVENT: {

#ifdef _DEBUG_REWRITE
          std::cout << "    -> no event] " << std::endl;
#endif

        } break;

        case SUPPRESSED: {

#ifdef _DEBUG_REWRITE
          std::cout << "    -> suppressed] " << std::endl;
#endif

          to_rewrite.remove(con.id());
          if (posted_constraints.contain(con.id()))
            posted_constraints.remove(con.id());

          // #ifdef _DEBUG_REWRITE
          // 	 std::cout << to_rewrite << std::endl;
          // 	 std::cout << posted_constraints << std::endl;
          // #endif

        } break;
        case FAIL_EVENT: {

#ifdef _DEBUG_REWRITE
          std::cout << "    -> fail event] " << std::endl;
#endif
          return FAILURE(0);
        } break;
        default: {

#ifdef _DEBUG_REWRITE
          std::cout << "   -> replaced] " << std::endl;
// exit(1);
#endif

          to_rewrite.remove(con.id());
          if (posted_constraints.contain(con.id()))
            posted_constraints.remove(con.id());

          // if(constraints[rewritten].rewritable()) {
          //   to_rewrite.add(rewritten);
          //   active_constraints.add(constraints[rewritten]);
          // }
        }
        }
      }
    }
#ifdef _DEBUG_REWRITE

    std::cout << "==>\n" << (*this) << std::endl;

#endif

#ifdef _DEBUG_REWRITE
    std::cout << " propagate: " << active_variables << std::endl;
    std::cout << "       and: " << active_constraints << std::endl;
#endif

    // fix_point = active_variables.empty();
  } while (!active_variables.empty() || !active_constraints.empty());

#ifdef _DEBUG_REWRITE
  std::cout << "\n====================END REWRITING======================\n";
  std::cout << active_variables.empty() << std::endl
            << active_constraints.empty() << std::endl;
#endif

  bool consistent = IS_OK(wiped_idx);
  wiped_idx = CONSISTENT;
  return consistent;

  // return IS_OK(wiped_idx);
}


bool Mistral::Solver::rewrite() 
{

  int i = 0;
  unsigned int con_i = 0;
  IntStack to_rewrite(0, 2 * constraints.size, false);
  Constraint con;
  RewritingOutcome rewritten;
  // bool fix_point;
  // Vector<Constraint> to_rewrite;

  // fix_point =

  // enforce AC
  // while( !fix_point ) {

#ifdef _DEBUG_REWRITE
  std::cout << "\n===================START REWRITING=====================\n";
#endif

  do {

#ifdef _DEBUG_REWRITE
    std::cout << " propagate  ";
#endif

    if (!propagate())
      break;

#ifdef _DEBUG_REWRITE
    std::cout << (*this) << std::endl;
// std::cout << "Collect rewritable constraints: " << std::endl;
#endif

    while (con_i < constraints.size) {
      // add the unseen constraints that might be rewritten
      for (; con_i < constraints.size; ++con_i) {

        // #ifdef _DEBUG_REWRITE
        //       std::cout << "   " << constraints[con_i] ;
        // #endif

        if (constraints[con_i].rewritable()) {
          to_rewrite.safe_add(con_i);

          // #ifdef _DEBUG_REWRITE
          // 	std::cout << " in" << std::endl;
          // #endif
        }

        // #ifdef _DEBUG_REWRITE
        //       else {
        // 	std::cout << " out" << std::endl;
        //       }
        // #endif
      }

      // // add all rewritable constraints to the stack
      // for(i=0; i<to_rewrite.size; ++i)
      //   active_constraints.add(constraints[to_rewrite[i]]);

#ifdef _DEBUG_REWRITE
      std::cout << " rewrite: " << to_rewrite << std::endl;
#endif

      // // rewriting
      // while(!active_constraints.empty()) {

      //   con = active_constraints.select(constraints);
      for (i = to_rewrite.size; i--;) {
        con = constraints[to_rewrite[i]];

#ifdef _DEBUG_REWRITE
        std::cout << "   [rewrite c" << con.id() << " " << con << std::endl;
#endif

        rewritten = con.rewrite();

        switch (rewritten) {
        case NO_EVENT: {

#ifdef _DEBUG_REWRITE
          std::cout << "    -> no event] " << std::endl;
#endif

        } break;

        case SUPPRESSED: {

#ifdef _DEBUG_REWRITE
          std::cout << "    -> suppressed] " << std::endl;
#endif

          to_rewrite.remove(con.id());
          if (posted_constraints.contain(con.id()))
            posted_constraints.remove(con.id());

          // #ifdef _DEBUG_REWRITE
          // 	 std::cout << to_rewrite << std::endl;
          // 	 std::cout << posted_constraints << std::endl;
          // #endif

        } break;
        case FAIL_EVENT: {

#ifdef _DEBUG_REWRITE
          std::cout << "    -> fail event] " << std::endl;
#endif
          return FAILURE(0);
        } break;
        default: {

#ifdef _DEBUG_REWRITE
          std::cout << "   -> replaced] " << std::endl;
// exit(1);
#endif

          to_rewrite.remove(con.id());
          if (posted_constraints.contain(con.id()))
            posted_constraints.remove(con.id());

          // if(constraints[rewritten].rewritable()) {
          //   to_rewrite.add(rewritten);
          //   active_constraints.add(constraints[rewritten]);
          // }
        }
        }
      }
    }
#ifdef _DEBUG_REWRITE

    std::cout << "==>\n" << (*this) << std::endl;

#endif

#ifdef _DEBUG_REWRITE
    std::cout << " propagate: " << active_variables << std::endl;
    std::cout << "       and: " << active_constraints << std::endl;
#endif

    // fix_point = active_variables.empty();
  } while (!active_variables.empty() || !active_constraints.empty());

#ifdef _DEBUG_REWRITE
  std::cout << "\n====================END REWRITING======================\n";
#endif

  bool consistent = IS_OK(wiped_idx);
  wiped_idx = CONSISTENT;
  return consistent;

  // return IS_OK(wiped_idx);
}

Mistral::PropagationOutcome Mistral::Solver::propagate(Constraint c, 
const bool force_trigger,
const bool trigger_self) {

  // std::cout << "solver.cpp: specific propagate(" << c << ")" << std::endl;

  // std::cout << active_variables << std::endl
  // 	    << active_constraints << std::endl;

  int trig, cons;
  bool fix_point;
  Triplet<int, Event, ConstraintImplementation *> var_evt;

  // wiped_idx = CONSISTENT;
  if (force_trigger)
    c.trigger();

  fix_point = (active_variables.empty() && active_constraints.empty());

  // std::cout << IS_OK(wiped_idx) << " & " << fix_point << std::endl;

  // std::cout << active_variables << std::endl;

  // std::cout << active_constraints << std::endl;

  while (IS_OK(wiped_idx) && !fix_point) {

    // empty the var stack first
    while (IS_OK(wiped_idx) && !active_variables.empty()) {

      var_evt = active_variables.pop_front();

      if (ASSIGNED(var_evt.second) &&
          sequence.contain(variables[var_evt.first])) {

        sequence.remove(variables[var_evt.first]);
        // last_solution_lb[var_evt.first] = last_solution_ub[var_evt.first] =
        // variables[var_evt.first].get_value();
        assignment_level[var_evt.first] = level;
        assignment_order[var_evt.first] = assignment_rank;
        ++assignment_rank;
      }

      // std::cout << var_evt << " "
      // 		<< variables[var_evt.first] << " \n0 "
      // 		<< constraint_graph[var_evt.first].on[0] << "\n1 "
      // 		<< constraint_graph[var_evt.first].on[1] << "\n2 "
      // 		<< constraint_graph[var_evt.first].on[2] << "\n"
      // 		<< c.propagator << std::endl;

      if (trigger_self || var_evt.third != c.propagator) {

        // for each triggered constraint
        for (trig = EVENT_TYPE(var_evt.second); IS_OK(wiped_idx) && trig < 3;
             ++trig) {

          for (cons = constraint_graph[var_evt.first].on[trig].size;
               IS_OK(wiped_idx) && --cons >= 0;) {

            culprit = constraint_graph[var_evt.first].on[trig][cons];

            if (culprit == c) {

              // if the constraints asked to be pushed on the constraint stack
              // we do that
              if (culprit.pushed()) {
                active_constraints.trigger(
                    (GlobalConstraint *)culprit.propagator, culprit.index(),
                    var_evt.second);
              }

              // if the constraint is not postponed, we propagate it
              if (!culprit.postponed()) {
                ++statistics.num_propagations;
                taboo_constraint = culprit.freeze();
                wiped_idx = culprit.propagate(var_evt.second);
                taboo_constraint = culprit.defrost();
              }
            }
          }
        }
      }
    }

    if (IS_OK(wiped_idx) && !active_constraints.empty()) {
      // propagate postponed constraint
      ++statistics.num_propagations;
      culprit = active_constraints.select(constraints);
      taboo_constraint = culprit.freeze();
      wiped_idx = culprit.propagate();
      taboo_constraint = culprit.defrost();
    } else if (active_variables.empty())
      fix_point = true;
  }

  taboo_constraint = NULL;
  active_constraints.clear();
  active_variables.clear();

  // std::cout << wiped_idx << std::endl;
  PropagationOutcome consistent = wiped_idx;
  wiped_idx = CONSISTENT;
  return consistent;
}

Mistral::PropagationOutcome Mistral::Solver::checker_propagate(Constraint c, 
const bool force_trigger,
const bool trigger_self) {
  int trig, cons;
  bool fix_point;
  Triplet<int, Event, ConstraintImplementation *> var_evt;

#ifdef _DEBUG_AC
  if (_DEBUG_AC) {
    std::cout << "c start (checker) propagation" << std::endl;
  }
#endif

  // wiped_idx = CONSISTENT;
  if (force_trigger)
    c.trigger();

  fix_point = (active_variables.empty() && active_constraints.empty());

#ifdef _DEBUG_AC
  int iteration = 0;
  if (_DEBUG_AC) {
    std::cout << std::endl;
  }
#endif

  while (IS_OK(wiped_idx) && !fix_point) {

#ifdef _DEBUG_AC
    if (_DEBUG_AC) {
      std::cout << "c ";
      std::cout << "var stack: " << active_variables << std::endl;
      ++iteration;
    }
#endif

    // empty the var stack first
    while (IS_OK(wiped_idx) && !active_variables.empty()) {

      var_evt = active_variables.pop_front();

#ifdef _DEBUG_AC
      if (_DEBUG_AC) {
        std::cout
            << "c "; // for(int lvl=0; lvl<iteration; ++lvl) std::cout << " ";
        std::cout << "react to " << event2str(var_evt.second) << " on "
                  << variables[var_evt.first] << " in "
                  << variables[var_evt.first].get_domain();
        if (var_evt.third)
          std::cout << " because of " << var_evt.third;
        std::cout << ". var stack: " << active_variables << std::endl;
      }
#endif

      if (ASSIGNED(var_evt.second) &&
          sequence.contain(variables[var_evt.first])) {
        sequence.remove(variables[var_evt.first]);
        // last_solution_lb[var_evt.first] = last_solution_ub[var_evt.first] =
        // variables[var_evt.first].get_value();
        assignment_level[var_evt.first] = level;
        assignment_order[var_evt.first] = assignment_rank;
        ++assignment_rank;
      }

      // std::cout << var_evt.third << " <-> " << c.propagator << std::endl;

      if (trigger_self || var_evt.third != c.propagator) {

        // std::cout << "cons list of " << variables[var_evt.first] << " = "
        // 	  << constraint_graph[var_evt.first].on << std::endl;

        // for each triggered constraint
        for (trig = EVENT_TYPE(var_evt.second); IS_OK(wiped_idx) && trig < 3;
             ++trig) {

          // std::cout << " trig[" << trig << "] = " <<
          // constraint_graph[var_evt.first].on[trig] << std::endl;

          for (cons = constraint_graph[var_evt.first].on[trig].size;
               IS_OK(wiped_idx) && --cons >= 0;) {

            culprit = constraint_graph[var_evt.first].on[trig][cons];

            // if(ASSIGNED(var_evt.second)) {
            //   culprit.notify_assignment();
            // }

#ifdef _DEBUG_AC
            if (_DEBUG_AC) {
              std::cout << "c "; // for(int lvl=0; lvl<iteration; ++lvl)
                                 // std::cout << " ";
              std::cout << "  -awake " << culprit << ": ";
              std::cout.flush();
            }
#endif

            if (culprit == c) {

              // if the constraint asks to be pushed on the constraint stack we
              // do that
              if (culprit.pushed()) {
#ifdef _DEBUG_AC
                if (_DEBUG_AC) {
                  std::cout << "pushed on the stack";
                }
#endif
                active_constraints.trigger(
                    (GlobalConstraint *)culprit.propagator, culprit.index(),
                    var_evt.second);
              }

              // if the constraint is not postponed, we propagate it
              if (!culprit.postponed()) {
#ifdef _DEBUG_AC
                if (_DEBUG_AC) {
                  if (IS_OK(wiped_idx)) {
                    Variable *scp = culprit.get_scope();
                    int arity = culprit.arity();
                    for (int i = 0; i < arity; ++i)
                      std::cout << scp[i].get_domain() << " ";
                    std::cout << "ok";
                  } else
                    std::cout << "fail";
                }
#endif
                ++statistics.num_propagations;
                taboo_constraint = culprit.freeze();

                // std::cout << "taboo: " << taboo_constraint << std::endl;

                wiped_idx = culprit.checker_propagate(var_evt.second);
                taboo_constraint = culprit.defrost();
              }
            }
#ifdef _DEBUG_AC
            else {
                                if(_DEBUG_AC) {
					std::cout << "c "; //for(int lvl=0; lvl<iteration; ++lvl) std::cout << " ";
                                        std::cout << "  -does not propagate "
                                                  << culprit
                                                  << " (we propagate " << c
                                                  << ")";
                                }
            }
            if (_DEBUG_AC) {
              std::cout << std::endl;
            }
#endif
          }
        }
      }
#ifdef _DEBUG_AC
      else {
        if (_DEBUG_AC) {
          std::cout
              << "c "; // for(int lvl=0; lvl<iteration; ++lvl) std::cout << " ";
          std::cout << "  -does not awake " << culprit << " (idempotent)"
                    << std::endl;
        }
      }
#endif
    }

    if (IS_OK(wiped_idx) && !active_constraints.empty()) {
      // propagate postponed constraint
      culprit = active_constraints.select(constraints);
      taboo_constraint = culprit.freeze();

      // if(taboo_constraint)
      // 	std::cout << "[constraint]" << std::endl;
      // else
      // 	std::cout << "[nill]" << std::endl;

      wiped_idx = culprit.checker_propagate();
      taboo_constraint = culprit.defrost();
    } else if (active_variables.empty())
      fix_point = true;
  }

  taboo_constraint = NULL;
  active_constraints.clear();
  active_variables.clear();

#ifdef _DEBUG_AC
  if (_DEBUG_AC) {
    std::cout << "c end (checker) propagation\n" << std::endl;
  }
#endif

  // return wiped_idx;
  PropagationOutcome consistent = wiped_idx;
  wiped_idx = CONSISTENT;
  return consistent;
}

Mistral::PropagationOutcome Mistral::Solver::bound_checker_propagate(Constraint c, 
const bool force_trigger,
const bool trigger_self) {

  // std::cout << "solver.cpp: bound checker propagate(" << c << ")" <<
  // std::endl;

  // std::cout << active_variables << std::endl
  // 	    << active_constraints << std::endl;

  int trig, cons;
  bool fix_point;
  Triplet<int, Event, ConstraintImplementation *> var_evt;
  // VarEvent var_evt;

  // wiped_idx = CONSISTENT;
  if (force_trigger)
    c.trigger();

  fix_point = (active_variables.empty() && active_constraints.empty());

  // std::cout << active_variables << std::endl;

  // std::cout << active_constraints << std::endl;

  // std::cout << (IS_OK(wiped_idx)) << " " << fix_point <<std::endl;

  while (IS_OK(wiped_idx) && !fix_point) {

    // empty the var stack first
    while (IS_OK(wiped_idx) && !active_variables.empty()) {

      var_evt = active_variables.pop_front();

      if (ASSIGNED(var_evt.second) &&
          sequence.contain(variables[var_evt.first])) {
        sequence.remove(variables[var_evt.first]);
        // last_solution_lb[var_evt.first] = last_solution_ub[var_evt.first] =
        // variables[var_evt.first].get_value();
        assignment_level[var_evt.first] = level;
        assignment_order[var_evt.first] = assignment_rank;
        ++assignment_rank;
      }

      // std::cout << var_evt << " "
      // 		<< variables[var_evt.first] << " \n0 "
      // 		<< constraint_graph[var_evt.first].on[0] << "\n1 "
      // 		<< constraint_graph[var_evt.first].on[1] << "\n2 "
      // 		<< constraint_graph[var_evt.first].on[2] << "\n"
      // 		<< c.propagator << std::endl;

      if (trigger_self || var_evt.third != c.propagator) {

        // std::cout << " " <<
        // constraint_graph[var_evt.first].on[EVENT_TYPE(var_evt.second)] <<
        // std::endl;

        // for each triggered constraint
        for (trig = EVENT_TYPE(var_evt.second); IS_OK(wiped_idx) && trig < 3;
             ++trig) {

          for (cons = constraint_graph[var_evt.first].on[trig].size;
               IS_OK(wiped_idx) && --cons >= 0;) {

            culprit = constraint_graph[var_evt.first].on[trig][cons];

            // std::cout << culprit << std::endl;

            if (culprit == c) {

              // if the constraints asked to be pushed on the constraint stack
              // we do that
              if (culprit.pushed()) {
                active_constraints.trigger(
                    (GlobalConstraint *)culprit.propagator, culprit.index(),
                    var_evt.second);
              }

              // if the constraint is not postponed, we propagate it
              if (!culprit.postponed()) {
                ++statistics.num_propagations;
                taboo_constraint = culprit.freeze();

                // if(taboo_constraint)
                // 	std::cout << "[constraint]" << std::endl;
                // else
                // 	std::cout << "[nill]" << std::endl;

                wiped_idx = culprit.bound_checker_propagate(var_evt.second);
                taboo_constraint = culprit.defrost();
              }
            }
          }
        }
      }
    }

    // std::cout << active_constraints << std::endl;
    // std::cout << active_constraints.empty() << std::endl;

    if (IS_OK(wiped_idx) && !active_constraints.empty()) {
      // propagate postponed constraint
      culprit = active_constraints.select(constraints);
      taboo_constraint = culprit.freeze();

      // if(taboo_constraint)
      // 	std::cout << "[constraint]" << std::endl;
      // else
      // 	std::cout << "[nill]" << std::endl;

      // std::cout << "solver.cpp: bound checker prop" << std::endl;

      wiped_idx = culprit.bound_checker_propagate();
      taboo_constraint = culprit.defrost();
    } else if (active_variables.empty())
      fix_point = true;
  }

  taboo_constraint = NULL;
  active_constraints.clear();
  active_variables.clear();

  // std::cout << wiped_idx << std::endl;

  // return wiped_idx;
  PropagationOutcome consistent = wiped_idx;
  wiped_idx = CONSISTENT;
  return consistent;
}


// Mistral::PropagationOutcome Mistral::Solver::propagate(Constraint *cons) {
//   culprit = cons;
//   //active_constraints.taboo_constraint = culprit->freeze();
//   active_constraints.select(culprit);
//   wiped_idx = culprit->propagate();
//   culprit->defrost();
//   if(IS_OK(wiped_idx)) culprit = NULL;
//   return wiped_idx;
// }

//#define _DEBUG_AC (statistics.num_nodes >= 489)

bool Mistral::Solver::propagate() 
{


#ifdef _DEBUG_AC
  if (_DEBUG_AC) {
    std::cout << "c start propagation" << std::endl;
  }
#endif

  bool fix_point;
  int trig, cons, vidx;
  Triplet<int, Event, ConstraintImplementation *> var_evt;

  // wiped_idx = CONSISTENT;
  culprit.clear();

  ++statistics.num_filterings;

  // TODO, we shouldn't have to do that
  if (IS_OK(wiped_idx) && objective && objective->enforce()) {
    wiped_idx = objective->objective.id();
  }

  fix_point = (active_variables.empty() && active_constraints.empty());

#ifdef _DEBUG_AC
  int iteration = 0;
  if (_DEBUG_AC) {
    std::cout << std::endl;
  }
#endif

  // if(!IS_OK(wiped_idx)) {

  //   std::cout << "FAILED WHEN ENFORCING UPPER BOUND" << std::endl;

  // }

  while (IS_OK(wiped_idx) && !fix_point) {

#ifdef _DEBUG_AC
    if (_DEBUG_AC) {

      std::cout
          << "c "; // for(int lvl=0; lvl<iteration; ++lvl) std::cout << " ";
      std::cout << "var stack: " << active_variables << std::endl;
      // std::cout << "c "; for(int lvl=0; lvl<iteration; ++lvl) std::cout << "
      // "; std::cout << "con stack: " << active_constraints << std::endl
      //   ;

      ++iteration;
    }
#endif

    assigned.clear();

    // empty the var stack first
    while (IS_OK(wiped_idx) && !active_variables.empty()) {

      // get the variable event
      var_evt = active_variables.pop_front();
      vidx = var_evt.first;

      assert(vidx >= 0);
      assert(vidx < is_relevant.size());
      is_relevant.set(vidx);

#ifdef _DEBUG_AC
      if (_DEBUG_AC) {
        std::cout
            << "c "; // for(int lvl=0; lvl<iteration; ++lvl) std::cout << " ";
        std::cout << "react to " << event2str(var_evt.second) << " on "
                  << variables[vidx] << " in " << variables[vidx].get_domain();
        if (var_evt.third)
          std::cout << " because of " << var_evt.third;
        std::cout << ". var stack: " << active_variables << std::endl;
      }
#endif

      if (ASSIGNED(var_evt.second)) {

        assigned.add(vidx);

        // std::cout << " -remove " << variables[var_evt.first] << " in " <<
        // variables[var_evt.first].get_domain() << std::endl;

        if (sequence.contain(variables[vidx]))
          sequence.remove(variables[vidx]);

        // last_solution_lb[var_evt.first] = last_solution_ub[var_evt.first] =
        // variables[var_evt.first].get_value();
        assignment_level[vidx] = level;
        assignment_order[vidx] = assignment_rank;
        ++assignment_rank;

        reason_for[vidx] = var_evt.third;
      }

      // for each triggered constraint
      for (trig = EVENT_TYPE(var_evt.second); IS_OK(wiped_idx) && trig < 3;
           ++trig) {

#ifdef _DEBUG_AC
        if (_DEBUG_AC) {
          std::cout << "trigger type: " << trig << ": "
                    << (constraint_graph[vidx].on[trig].size) << " constraints"
                    << std::endl;
        }
#endif

        for (cons = constraint_graph[vidx].on[trig].size;
             IS_OK(wiped_idx) && --cons >= 0;) {

          culprit = constraint_graph[vidx].on[trig][cons];

#ifdef _DEBUG_AC
          if (_DEBUG_AC) {
            std::cout << " -- " << culprit << std::endl;
          }
#endif

          // idempotency, if the event was triggered by itself
          if (var_evt.third != culprit.propagator) {
#ifdef _DEBUG_AC
            if (_DEBUG_AC) {
              // std::cout << variables[86] << " in " <<
              // variables[86].get_domain() << "  " << variables[17] << " in "
              // << variables[17].get_domain() << " ";
              std::cout << "c [" << culprit.id()
                        << "]"; // for(int lvl=0; lvl<iteration; ++lvl)
                                // std::cout << " ";
              std::cout << "  -awake " << culprit << ": ";
              std::cout.flush();
            }
#endif

            // if the constraints asked to be pushed on the constraint stack we
            // do that
            if (culprit.pushed()) {

#ifdef _DEBUG_AC
              if (_DEBUG_AC) {
                std::cout << "pushed on the stack";
              }
#endif
              active_constraints.trigger((GlobalConstraint *)culprit.propagator,
                                         culprit.index(), var_evt.second);
            }

            // if the constraint is not postponed, we propagate it
            if (!culprit.postponed()) {

              if (ASSIGNED(var_evt.second)) {
                culprit.notify_assignment();
              }

#ifdef _DEBUG_AC
              if (_DEBUG_AC) {
                if (culprit.pushed())
                  std::cout << ", ";
                std::cout << "propagated: ";
                Variable *scp = culprit.get_scope();
                int arity = culprit.arity();
                for (int i = 0; i < arity; ++i)
                  std::cout << scp[i].get_domain() << " ";
                std::cout << "-> ";
                std::cout.flush();
              }
#endif
              ++statistics.num_propagations;
              taboo_constraint = culprit.freeze();
              wiped_idx = culprit.propagate(var_evt.second);
              taboo_constraint = culprit.defrost();

#ifdef _DEBUG_AC
              if (_DEBUG_AC) {
                if (IS_OK(wiped_idx)) {
                  Variable *scp = culprit.get_scope();
                  int arity = culprit.arity();
                  for (int i = 0; i < arity; ++i)
                    std::cout << scp[i].get_domain() << " ";
                  std::cout << "ok";
                } else
                  std::cout << "fail";
              }
#endif

#ifdef _REF_SOL_
              check_reference();
              if (!IS_OK(wiped_idx)) {
                std::cout << "wrong fail\n";
                exit(1);
              }
#endif
            }
#ifdef _DEBUG_AC
            if (_DEBUG_AC) {
              std::cout << std::endl;
            }
#endif
          }
#ifdef _DEBUG_AC
          else {
            if (_DEBUG_AC) {
              std::cout << "c [" << culprit.id()
                        << "]"; // for(int lvl=0; lvl<iteration; ++lvl)
                                // std::cout << " ";
              std::cout << "  -does not awake " << culprit << " (idempotent)"
                        << std::endl;
            }
          }
#endif
        }
      }
    }

    if (IS_OK(wiped_idx)) {
      while (!assigned.empty()) {
        vidx = assigned.pop();
        for (trig = 0; trig < 3; ++trig) {
          for (cons = constraint_graph[vidx].on[trig].size; --cons >= 0;) {
            culprit = constraint_graph[vidx].on[trig][cons];
            if (culprit.postponed()) {
              culprit.notify_assignment();
            }
          }
        }
      }
    }

    if (IS_OK(wiped_idx) && !active_constraints.empty()) {

#ifdef _DEBUG_AC
      if (_DEBUG_AC) {
        if (level > 0) {
          std::cout << "\nc propagate postponed constraints: "
                    << active_constraints << std::endl;
        } else {
          std::cout << "\nc propagate postponed constraints: "
                    << active_constraints._set_ << std::endl;
        }
      }
#endif

      // propagate postponed constraint
      culprit = active_constraints.select(constraints);

#ifdef _DEBUG_AC
      if (_DEBUG_AC) {
        std::cout
            << "c [" << culprit.id()
            << "]"; // for(int lvl=0; lvl<iteration; ++lvl) std::cout << " ";
        std::cout << "  -propagate " << culprit << " (";

        if (culprit.global()) {
          GlobalConstraint *gc = (GlobalConstraint *)(culprit.propagator);

          std::cout << gc->events << " ";
          std::cout.flush();

          std::cout << event2str(gc->event_type[gc->events[0]]) << " on "
                    << gc->scope[gc->events[0]];
          for (unsigned int i = 1; i < gc->events.size; ++i) {
            std::cout << ", " << event2str(gc->event_type[gc->events[i]])
                      << " on " << gc->scope[gc->events[i]];
          }
        }
        std::cout << ") ";
        Variable *scp = culprit.get_scope();
        int arity = culprit.arity();
        for (int i = 0; i < arity; ++i)
          std::cout << scp[i].get_domain() << " ";
        std::cout << "-> ";
        std::cout.flush();
      }
#endif

      ++statistics.num_propagations;
      taboo_constraint = culprit.freeze();
      wiped_idx = culprit.propagate();
      taboo_constraint = culprit.defrost();

#ifdef _DEBUG_AC
      if (_DEBUG_AC) {
        if (IS_OK(wiped_idx)) {
          Variable *scp = culprit.get_scope();
          int arity = culprit.arity();
          for (int i = 0; i < arity; ++i)
            std::cout << scp[i].get_domain() << " ";
          std::cout << "ok";
        } else
          std::cout << "fail";
        std::cout << std::endl;
      }
#endif

#ifdef _REF_SOL_
      check_reference();
      if (!IS_OK(wiped_idx)) {
        std::cout << "wrong fail\n";
        exit(1);
      }
#endif

    } else if (active_variables.empty())
      fix_point = true;
  }

  taboo_constraint = NULL;
  active_constraints.clear();
  if (!parameters.backjump) {
    active_variables.clear();
  }

#ifdef _DEBUG_AC
  if (_DEBUG_AC) {
    if (!IS_OK(wiped_idx)) {
      std::cout << "inconsistency found!" << std::endl;
    } else {
      std::cout << "done" << std::endl;
    }
  }
#endif 


#ifdef _DEBUG_AC
  if (_DEBUG_AC) {
    std::cout << "c end propagation" << std::endl;
  }
#endif

  if (IS_OK(wiped_idx)) {
#ifdef _DEBUG_SEARCH
    if (_DEBUG_SEARCH) {
      std::cout << parameters.prefix_comment;
      for (unsigned int k = 0; k <= decisions.size; ++k)
        std::cout << " ";
      std::cout << "success! "
                //<< sequence
                << std::endl;
    }
#endif
    culprit.clear();
    notify_success();
    return true;
  } else {
    ++statistics.num_failures;

#ifdef _DEBUG_SEARCH
    if (_DEBUG_SEARCH) {
      std::cout << parameters.prefix_comment;
      for (unsigned int k = 0; k <= decisions.size; ++k)
        std::cout << " ";
      std::cout << "failure!" << std::endl;
    }
#endif

#ifdef _DEBUG_FAIL
    std::cout << "failure of c" << culprit.id() << " " << culprit << std::endl;
#endif

    // std::cout << "solver: notify failure" << std::endl;

    if (parameters.backjump)
      close_propagation();

    prev_wiped_idx = wiped_idx;
    // notify_failure();
    wiped_idx = CONSISTENT;
    return false;
  }
}


void Mistral::Solver::fail() {

#ifdef _DEBUG_SEARCH
  if (_DEBUG_SEARCH) {
    std::cout << parameters.prefix_comment;
    for (unsigned int k = 0; k <= decisions.size; ++k)
      std::cout << " ";
    std::cout << "failurex!" << std::endl;
  }
#endif

  wiped_idx = FAILURE(0);
  prev_wiped_idx = wiped_idx;
}

// bool Mistral::Solver::rewrite() 
// {

//   bool fix_point;
//   int trig, cons;
//   Triplet < int, Event, ConstraintImplementation* > var_evt;

//   wiped_idx = CONSISTENT;
//   culprit.clear();

//   ++statistics.num_filterings;  

//   // TODO, we shouldn't have to do that
//   if(objective && objective->enforce())
//     wiped_idx = objective->objective.id();

//   fix_point =  (active_variables.empty() && active_constraints.empty());


// #ifdef _DEBUG_AC
//   int iteration = 0;
//   std::cout << std::endl;
// #endif

//   while(IS_OK(wiped_idx) && !fix_point) {

// #ifdef _DEBUG_AC

//     std::cout << "c "; for(int lvl=0; lvl<iteration; ++lvl) std::cout << " ";
//     std::cout << "var stack: " << active_variables << std::endl;
//     // std::cout << "c "; for(int lvl=0; lvl<iteration; ++lvl) std::cout << " ";
//     // std::cout << "con stack: " << active_constraints << std::endl
//     //   ;

//     ++iteration;
// #endif

//     // empty the var stack first
//     while( IS_OK(wiped_idx) && 
// 	   !active_variables.empty() ) {

//       // get the variable event
//       var_evt = active_variables.pop_front();

// #ifdef _DEBUG_AC
//       std::cout << "c "; for(int lvl=0; lvl<iteration; ++lvl) std::cout << " ";
//       std::cout << "react to " << event2str(var_evt.second) << " on " << variables[vidx] 
// 		<< " in " << variables[vidx] ;
//       if(var_evt.third)
// 	std::cout << " because of " << var_evt.third ;
//       std::cout << ". var stack: " << active_variables << std::endl;
// #endif

//       if(ASSIGNED(var_evt.second) && sequence.contain(variables[vidx])) {
// 	sequence.remove(variables[vidx]);
// 	assignment_level[vidx] = level;
//       }

//       // for each triggered constraint
//       for(trig = EVENT_TYPE(var_evt.second); IS_OK(wiped_idx) &&
// 	    trig<3; ++trig) {

// 	// std::cout << variables[vidx] << ": " << constraint_graph[vidx].on[trig] << std::endl
// 	// 	  << constraint_graph[vidx].on[trig].size << std::endl;

// 	for(cons = constraint_graph[vidx].on[trig].size; IS_OK(wiped_idx) &&
// 	      --cons>=0;) {

// 	  culprit = constraint_graph[vidx].on[trig][cons];

// 	  // idempotency, if the event was triggered by itself
// 	  if(var_evt.third != culprit.propagator) {
// #ifdef _DEBUG_AC
// 	    std::cout << "c "; for(int lvl=0; lvl<iteration; ++lvl) std::cout << " ";
// 	    std::cout << "  -awake " << culprit << ": "; 
// #endif
// 	    // if the constraints asked to be pushed on the constraint stack we do that
// 	    if(culprit.pushed()) {
// #ifdef _DEBUG_AC
// 	      std::cout << "pushed on the stack" ;
// #endif    
// 	      active_constraints.trigger((GlobalConstraint*)culprit.propagator, 
// 					 culprit.index(), var_evt.second);
// 	    }
// 	    // if the constraint is not postponed, we propagate it
// 	    if(!culprit.postponed()) {
// #ifdef _DEBUG_AC
// 	      if(culprit.pushed()) std::cout << ", ";
// 	      std::cout << "propagated: ";
// 	      Variable *scp = culprit.get_scope();
// 	      int arity = culprit.arity();
// 	      for(int i=0; i<arity; ++i)
// 		std::cout << scp[i].get_domain() << " ";
// 	      std::cout << "-> ";
// #endif
// 	      ++statistics.num_propagations;  
// 	      taboo_constraint = culprit.freeze();
// 	      wiped_idx = culprit.propagate(var_evt.second); 
// 	      taboo_constraint = culprit.defrost();

// 	      if(IS_OK(wiped_idx)) {
// 		culprit.rewrite();
// 	      }

// #ifdef _DEBUG_AC
// 	      if(IS_OK(wiped_idx)) {
// 		Variable *scp = culprit.get_scope();
// 		int arity = culprit.arity();
// 		for(int i=0; i<arity; ++i)
// 		  std::cout << scp[i].get_domain() << " ";
// 		std::cout << "ok";
// 	      } else std::cout << "fail";
// #endif
// 	    }
// #ifdef _DEBUG_AC
// 	    std::cout << std::endl; 
// #endif    
// 	  } 
// #ifdef _DEBUG_AC
// 	  else {
// 	    std::cout << "c "; for(int lvl=0; lvl<iteration; ++lvl) std::cout << " ";
// 	    std::cout << "  -does not awake " << culprit << " (idempotent)" << std::endl; 
// 	  }
// #endif
// 	}
//       }
//     }


// #ifdef _DEBUG_AC
//     //std::cout << std::endl;
//     std::cout << "c "; for(int lvl=0; lvl<iteration; ++lvl) std::cout << " ";
//     //std::cout << " var stack: " << active_variables << std::endl;
//     //std::cout << "c2 "; for(int lvl=0; lvl<level; ++lvl) std::cout << " ";
//     //std::cout << "con stack: " << active_constraints << std::endl
//       ;
// #endif

//     if(IS_OK(wiped_idx) && !active_constraints.empty()) {

// // #ifdef _DEBUG_AC
// //       std::cout << "\npropagate postponed constraints: " 
// // 		<< active_constraints << std::endl;
// // #endif

//       // propagate postponed constraint
//       culprit = active_constraints.select(constraints);

// #ifdef _DEBUG_AC
//       std::cout << "c "; for(int lvl=0; lvl<iteration; ++lvl) std::cout << " ";
//       std::cout << "  -propagate " << culprit << " (" ;

//       if(culprit.global()) {
// 	GlobalConstraint *gc = (GlobalConstraint*)(culprit.propagator);

// 	std::cout << gc->events << " ";
// 	std::cout.flush();

// 	std::cout << event2str(gc->event_type[gc->events[0]]) << " on " << gc->scope[gc->events[0]] ;
// 	for(unsigned int i=1; i<gc->events.size; ++i) {
// 	  std::cout << ", " << event2str(gc->event_type[gc->events[i]]) << " on " << gc->scope[gc->events[i]] ;
// 	}
//       }
//       std::cout << ") ";
//       Variable *scp = culprit.get_scope();
//       int arity = culprit.arity();
//       for(int i=0; i<arity; ++i)
// 	std::cout << scp[i].get_domain() << " ";
//       std::cout << "-> ";

// #endif

//       taboo_constraint = culprit.freeze();
//       wiped_idx = culprit.propagate(); 
//       taboo_constraint = culprit.defrost();

// #ifdef _DEBUG_AC
//       if(IS_OK(wiped_idx)) {
// 	Variable *scp = culprit.get_scope();
// 	int arity = culprit.arity();
// 	for(int i=0; i<arity; ++i)
// 	  std::cout << scp[i].get_domain() << " ";
// 	std::cout << "ok";
//       } else std::cout << "fail";
//       std::cout << std::endl;
// #endif

//     } else if(active_variables.empty()) fix_point = true;
//   }

//   taboo_constraint = NULL;
//   active_constraints.clear();
//   if(!parameters.backjump) {
//     active_variables.clear();
//   }

// #ifdef _DEBUG_AC
//   if(!IS_OK(wiped_idx)) {
//     std::cout << "inconsistency found!" << std::endl;
//   } else {
//     std::cout << "done" << std::endl;
//   }
// #endif

//   if(IS_OK(wiped_idx)) {
//     notify_success();   
//     return true;
//   } else {
//     ++statistics.num_failures;

//     //std::cout << "solver: notify failure" << std::endl;

//     notify_backtrack();
//     return false;
//   }
// }



// bool Mistral::Solver::propagate() 
// {

//   wiped_idx = -1;
//   culprit = NULL;

//   // if there is an objective function, it is propagated first, 
//   // since upon backtrack, it could fail even if n-1 variables are set
//   if(objective && objective->enforce())
//     wiped_idx = objective->objective.id();

//   ++statistics.num_filterings;

//   int *v_index;
//   Constraint **c_iter, **end, *c;
//   bool finished = active_constraints.empty() && active_variables.empty();
//   while( IS_OK(wiped_idx) && (!finished) ) {

//     while(!active_variables.empty()) {      
//       var = active_variables.pop( evt );

//       c_iter = constraint_graph_array[var]->first(EVENT_TYPE(evt));
//       v_index = constraint_graph_array[var]->get_index(c_iter);
//       end = constraint_graph_array[var]->first();

//       if(ASSIGNED(evt)) {
// 	while(c_iter < end) {
// 	  active_constraints.trigger(*c_iter, *v_index, evt);
//        	  c = (*c_iter)->notify_assignment(*v_index, level);
//        	  if(c) saved_cons.add(c);
// 	  ++c_iter;
// 	  ++v_index;
//       	}
//       	if(sequence.contain(variables[var])) {
//       	  sequence.remove(variables[var]);
//       	  assignment_level[var] = level;
//       	}
//       } else {
// 	//if(!(ASSIGNED(evt)))
// 	while(c_iter < end) {
//   	  active_constraints.trigger(*c_iter, *v_index, evt);
// 	  ++c_iter;
// 	  ++v_index;
// 	}
//       }
//     }

//     if(IS_OK(wiped_idx) && !active_constraints.empty()) {

//       culprit = active_constraints.select(constraints);

//       ++statistics.num_propagations;

//       wiped_idx = culprit->propagate();

//       culprit->defrost();

//      }

//   }

//   //  taboo_constraint = NULL;
//   active_constraints.clear();
//   active_variables.clear();

//   if(IS_OK(wiped_idx)) {
//     notify_success();
//   } else {
//     ++statistics.num_failures;
//     notify_failure();
//   }
//   return IS_OK(wiped_idx); //!wiped_out;
// }


#ifdef _MONITOR
void Mistral::Solver::monitor(Vector< Variable >& X) {
  for (unsigned int i = 0; i < X.size; ++i)
    monitored.add(X[i].id());
  monitored_index.add(monitored.size);
}
void Mistral::Solver::monitor(Variable X) {
  monitored.add(X.id());
  monitored_index.add(monitored.size);
}
#endif


void Mistral::Solver::full_print() {
  for (int j = 0; j < level; ++j)
    std::cout << "  ";
  std::cout << sequence << std::endl;
  for (unsigned int i = 0; i < variables.size; ++i) {
    for (int j = 0; j < level; ++j)
      std::cout << "  ";
    // variables[i].full_print();
    std::cout << std::endl;
  }
}

void Mistral::Solver::debug_print() {

  //   std::cout << std::endl << this << std::endl;

  // //   std::cout << "    currently changed: " << std::endl;
  // //   int k = changed_objs.size;
  // //   while(k --> 0) {
  // //     changed_objs[k]->debug_print();
  // //     //std::cout << std::endl;
  // //   }
  // //   std::cout << std::endl;
  //   int i, k;
  //   for(i=var_trail_size.size-1; i>0; --i) {
  //     std::cout << "    changed at level " << (i-1) << ": " << std::endl;
  //     k = var_trail_size[i];
  //     while(k --> var_trail_size[i-1])
  //       saved_vars[k]->debug_print();
  //     //std::cout << std::endl;
  //   }

  //   std::cout << "    AC Queue: " << std::endl;
  // //   std::cout << active_constraints.active << std::endl;
  // //   std::cout << active_constraints.triggers[0] << std::endl;
  // //   std::cout << active_constraints.triggers[1] << std::endl;
  // //   std::cout << active_constraints.triggers[2] << std::endl;

  //   for(unsigned int i=0; i<3; ++i) {
  //     if(!active_constraints.triggers[i].empty()) {
  //       std::cout << "priority " << i << std::endl;
  //       int elt = active_constraints.triggers[i].first();
  //       while(elt != active_constraints.triggers[i]._head) {
  // 	Constraint *cons = constraints[elt];
  // 	std::cout << "\t" << cons << ": ";
  // 	for(unsigned int j=0; j<cons->changes.size; ++j) {
  // 	  int var = cons->changes[j];
  // 	  int evt = cons->event_type[var];
  // 	  std::cout << cons->scope[var] << "/"
  // 		    << (is_value(evt) ? "v" : (is_range(evt) ? "r" : "d") )
  // 		    << (is_upper_bound(evt) ? "u" : "")
  // 		    << (is_lower_bound(evt) ? "l" : "")
  // 		    << " ";
  // 	}
  // 	std::cout << std::endl;
  // 	elt = active_constraints.triggers[i].next[elt];
  //       }
  //     }
  //   }

  // //   for(unsigned int i=0; i<active_constraints.size; ++i) {
  // //     Constraint *cons = constraints[active_constraints[i]];
  // //     std::cout << cons << ": ";
  // //     for(unsigned int j=0; j<cons->changes.size; ++j) {
  // //       int var = cons->changes[j];
  // //       int evt = cons->event_type[var];
  // //       std::cout << cons->scope[var] << "/"
  // // 		<< (is_value(evt) ? "v" : (is_range(evt) ? "r" : "d") )
  // // 		<< (is_upper_bound(evt) ? "u" : "")
  // // 		<< (is_lower_bound(evt) ? "l" : "")
  // // 		<< " ";
  // //     }
  // //     std::cout << std::endl;
  // //   }
}

// std::string Mistral::Solver::getString() const {
//   std::string return_str = "Variables:\n";
//   for(unsigned int i=0; i<variables.size; ++i)
//     return_str += ("\t"+toString(variables[i])+" in "+toString(variables[i]->domain)+"\n");

//   return_str += "\nConstraints:\n";
//   for(unsigned int i=0; i<constraints.size; ++i)
//     return_str += ("\t"+toString(constraints[i])+"\n");

//   return_str += ("\nSearch on "+toString(search.sequence)+"\n");

//   return return_str;
// }

void Mistral::Solver::print_clist(int k) {
  for (Event trig = 0; trig < 3; ++trig)
    for (int cons = constraint_graph[k].on[trig].size; --cons >= 0;)
      std::cout << "[" << constraint_graph[k].on[trig][cons].id() << "]";
  std::cout << "\n";
}

std::ostream& Mistral::Solver::display(std::ostream& os, const int current) {

  if (objective)
    os << objective << std::endl;

  os << "Variables:\n";

  Variable *scope;
  int arity;

  Vector<Variable> rem_vars;
  for (unsigned int i = 0; i < variables.size; ++i) {
    if (!(domain_types[i] & REMOVED_VAR) &&
        (current != 1 || sequence.contain(i)) &&
        (current != 2 || !(variables[i].is_ground()))) {

      os << "  " << variables[i] << " in "
         << variables[i].get_domain(); //<< "\n";

      os << ": ";
      for (Event trig = 0; trig < 3; ++trig)
        for (int cons = constraint_graph[i].on[trig].size; --cons >= 0;) {
          if (current) {
            os << "[" << constraint_graph[i].on[trig][cons].id()
               << constraint_graph[i].on[trig][cons].symbol();

            scope = constraint_graph[i].on[trig][cons].get_scope();
            arity = constraint_graph[i].on[trig][cons].arity();

            int k = 0, j = 0;
            for (; j < arity && k < 1; ++j) {
              if (!scope[j].is_ground()) {
                ++k;
                std::cout << scope[j].id() << " ";
              }
            }
            for (; j < arity && k < 2; ++j) {
              if (!scope[j].is_ground()) {
                ++k;
                std::cout << scope[j].id(); //<< " " ;
              }
            }

            os
                //<< constraint_graph[i].on[trig][cons].symbol()
                // constraint_graph[i].on[trig][cons].id()
                << "]";
          } else {
            os << "[" <<
                // constraint_graph[i].on[trig][cons].symbol()
                constraint_graph[i].on[trig][cons].id() << "]";
          }
        }
      // if(lit_activity)
      // 	os << lit_activity[2*i] << "/" << lit_activity[2*i+1] << ": " <<
      // var_activity[i] ;
      os << "\n";

    } else {

      rem_vars.add(variables[i]);
    }
  }

  if (!current) {
    if (!rem_vars.empty()) {
      os << "  (suppressed:";
      for (unsigned int i = 0; i < rem_vars.size; ++i)
        os << " " << rem_vars[i];
      os << ")";
    }

    os << "\nConstraints:\n" << posted_constraints << std::endl;
    for (unsigned int i = 0; i < posted_constraints.size; ++i) {

      // os << posted_constraints[i] << std::endl;

      os << "  [" << constraints[posted_constraints[i]].id()
         << "]: " << constraints[posted_constraints[i]] << "\n";
    }
  }

  // os << lit_activity << " " << lit_activity[0] << " " << lit_activity[1] <<
  // std::endl;

  return os;
}

//#define _DEBUG_EXPLANATION

// Mistral::Decision
void Mistral::Solver::learn_nogood() {

  visited_literals.clear();

  // unsigned int vidx;

  // // first, resolve the unresolved events, so that 'reason_for' and
  // 'assignment_level' are up to date Triplet < int, Event,
  // ConstraintImplementation* > var_evt; while(!active_variables.empty()) {
  //   var_evt = active_variables.pop_front();
  //   vidx = var_evt.first;
  //   if(ASSIGNED(var_evt.second) && sequence.contain(variables[vidx])) {
  //     sequence.remove(variables[vidx]);
  //     assignment_level[vidx] = level;
  //     assignment_order[vidx] = assignment_rank;
  //     ++assignment_rank;
  //     reason_for[vidx] = var_evt.third; //->explain();
  //   }
  // }

  int pathC = 0, index = sequence.size - 1;
  Literal p = 0, q;
  Atom a = NULL_ATOM;
  Variable x;
  int lvl;

  // double *lit_activity = base->lit_activity.stack_;
  // double *var_activity = base->var_activity.stack_;

  // We start from the constraint that failed
  Explanation *current_explanation = culprit.propagator;

  // Variable *scope = culprit.get_scope();
  // int arity = culprit.arity();
  // for(int i=0; i<arity; ++i) {
  //   var_activity[scope[i].id()] += 10 * parameters.activity_increment;
  // }

#ifdef _DEBUG_NOGOOD
  int vidx;
  if (_DEBUG_NOGOOD) {
    // int depth = 0;
    std::cout << "\nExplain failure (sequence = [\n";
    for (vidx = sequence.size; vidx < variables.size; ++vidx) {
      int var_idx = sequence[vidx].id();
      std::cout << std::setw(3) << assignment_order[var_idx] << " "
                << std::setw(3) << assignment_level[var_idx] << " "
                << sequence[vidx] << " == " << sequence[vidx].get_value()
                << " <- ";

      std::cout.flush();

      Explanation *expl = reason_for[sequence[vidx].id()];

      if (!expl) {
        if (assignment_level[var_idx])
          std::cout << "decision";
        else
          std::cout << "deduction";
      } else if (expl->is_clause()) {
        std::cout << "l: ";
        print_clause(std::cout, (Clause *)expl);
      } else if (expl != base) {
        std::cout << parameters.prefix_comment
                  << ((ConstraintImplementation *)(expl))->id << ": (";

        std::cout.flush();

        Explanation::iterator stop;
        Explanation::iterator lit =
            expl->get_reason_for(var_idx, assignment_level[var_idx], stop);

        if (lit != stop) {

          print_literal(std::cout, *lit);
          ++lit;

          while (lit < stop) {
            std::cout << " v ";
            print_literal(std::cout, *lit);
            ++lit;
          }
        }

        std::cout << ")";

      } else {
        std::cout << "b: ";
        print_clause(std::cout, base->reason_for[var_idx]);
      }

      std::cout << std::endl;
    }
    std::cout << " ])\n";
    //<< sequence << ")\n";
  }
#endif

  backtrack_level = 0;

  // the resulting nogood is stored in the vector 'learnt_clause'
  learnt_clause.clear();
  learnt_clause.add(p);

  do {
    // add the parents of the conflict to the current set of visited atoms

#ifdef _DEBUG_NOGOOD
    if (_DEBUG_NOGOOD) {
      std::cout << " reason: ";

      std::cout.flush();

      if (!current_explanation) {
        if (assignment_level[a])
          std::cout << "decision";
        else
          std::cout << "deduction";
      } else if (current_explanation->is_clause()) {
        std::cout << "l: ";
        print_clause(std::cout, (Clause *)current_explanation);
      } else if (current_explanation != base) {
        std::cout << parameters.prefix_comment
                  << ((ConstraintImplementation *)(current_explanation))->id
                  << ": (";

        // std::cout << "\nhere\n";

        // std::cout << current_explanation << std::endl;

        Explanation::iterator stop;
        Explanation::iterator lit = current_explanation->get_reason_for(
            a, ((a != NULL_ATOM) ? assignment_level[a] : level), stop);

        //      std::cout << "\nthere\n";

        std::cout.flush();

        if (lit != stop) {

          print_literal(std::cout, *lit);
          ++lit;

          while (lit < stop) {
            std::cout << " v ";
            print_literal(std::cout, *lit);
            ++lit;
          }
        }

        std::cout << ")";

      } else {

        std::cout << "b: ";
        if (a == NULL_ATOM)
          print_clause(std::cout, base->conflict);
        else
          print_clause(std::cout, base->reason_for[a]);
      }

      std::cout << std::endl;

      // if(current_explanation != base) {
      //   std::cout << "*" << //((ConstraintImplementation*)
      // 	current_explanation
      // 	//)->id
      // 		<< std::endl;
      // } else if(a == NULL_ATOM) {
      //   print_clause(std::cout, base->conflict);
      //   std::cout <<  std::endl;
      // } else {
      //   print_clause(std::cout, base->reason_for[a]);
      //   std::cout << std::endl;
      // }
      //++depth;
    }
#endif

    // Explanation::iterator lit = current_explanation->begin(a);
    // Explanation::iterator stop = current_explanation->end(a);

    if (a == NULL_ATOM || assignment_level[a]) {

      Explanation::iterator stop;

      if (current_explanation == NULL) {
        std::cout << "NULL POINTER!!! " << (statistics.num_filterings)
                  << std::endl;
        exit(1);
      }

      // std::cerr << a << std::endl;

      // std::cerr << assignment_level[a] << std::endl;

#ifdef _CHECK_NOGOOD

      // std::cout << (int*)current_explanation << " " ;
      // std::cout.flush();
      // std::cout << current_explanation << std::endl;

      store_reason(current_explanation, a);
#endif

      Explanation::iterator lit = current_explanation->get_reason_for(
          a, (a != NULL_ATOM ? assignment_level[a] : level), stop);

      while (lit < stop) {
        q = *lit;
        ++lit;

        a = UNSIGNED(q);
        x = variables[a];
        lvl = assignment_level[a];

        // if(a == 5247) {
        // std::cout << a << "  cureent level=" << level << " / asgnmnt lvl=" <<
        // lvl << std::endl;
        // }

#ifdef _DEBUG_NOGOOD
        if (_DEBUG_NOGOOD) {
          // for(int i=0; i<depth; ++i)
          std::cout << "   " << lvl << " ";
          print_literal(std::cout, q);
          std::cout << ": ";
        }
#endif

        if (!visited.fast_contain(a)) {
          // if(lit_activity) {

          // #ifdef _DEBUG_ACTIVITY
          // 	    std::cout << "x" << a << " (" << lit_activity[NEG(a)] << "/"
          // << lit_activity[POS(a)] << ")/" << var_activity[a]
          // 		   << " -> " ;
          // #endif

          // //lit_activity[q] += 0.5 * parameters.activity_increment;
          // lit_activity[NOT(q)] += // 0.5 *
          //   parameters.activity_increment;
          // var_activity[a] += parameters.activity_increment;
          visited_literals.add(NOT(q));

          // #ifdef _DEBUG_ACTIVITY
          // 	    std::cout << " (" << lit_activity[NEG(a)] << "/" <<
          // lit_activity[POS(a)] << ")/" << var_activity[a]
          // 		   << std::endl ;
          // #endif

          //}
          visited.fast_add(a);

          if (lvl >= level) {
// we'll need to replace 'a' by its parents since its level is too high
#ifdef _DEBUG_NOGOOD
            if (_DEBUG_NOGOOD) {
              std::cout << "to expend later" << std::endl;
            }
#endif
            ++pathC;
          } else {
            // q's level is below the current level, we are not expending it
            // further
            learnt_clause.add(q);
#ifdef _DEBUG_NOGOOD
            if (_DEBUG_NOGOOD) {
              std::cout << "add to the clause";
              // std::cout << ": {" ;
              // for(unsigned int k=1; k<learnt_clause.size; ++k) {
              //   std::cout << " ";
              //   print_literal(std::cout, learnt_clause[k]);
              // }
              // std::cout << " }" ;
              std::cout << std::endl;
            }
#endif
            if (lvl > backtrack_level)
              backtrack_level = lvl;
          }
        }
#ifdef _DEBUG_NOGOOD
        else {
          if (_DEBUG_NOGOOD) {
            std::cout << "visited " // << visited
                      << std::endl;
          }
        }
#endif
      }
    }

#ifdef _DEBUG_NOGOOD
    if (_DEBUG_NOGOOD) {
      std::cout << "current literals to expend:";
      for (vidx = index + 1; vidx < variables.size; ++vidx) {
        if (visited.fast_contain(sequence[vidx].id())) {
          std::cout << " " << sequence[vidx];
        }
      }
      std::cout << std::endl;
    }
#endif

    // jump to the next visited atom that need be further expended
    while (!visited.fast_contain(sequence[++index].id())) {

#ifdef _DEBUG_NOGOOD
      if (_DEBUG_NOGOOD) {
        if (index >= variables.size - 1) {
          std::cout << "reached the end of the stack!!" << std::endl;
        }
      }
#endif
    };

    x = sequence[index];
    a = x.id();
    p = ((2 * a) | (x.get_min()));
    lvl = assignment_level[a];

    if (pathC > 1) {
      // there are still atoms to expend, we start with 'a'

      // EXPL
      current_explanation = reason_for[a];
      visited.fast_add(a);

#ifdef _DEBUG_NOGOOD
      if (_DEBUG_NOGOOD) {
        // for(int i=0; i<depth; ++i) std::cout << "  ";
        std::cout << pathC << " - replace ";
        print_literal(std::cout, p);
        std::cout << " by its ";
        std::cout.flush();
      }
#endif

    }
#ifdef _DEBUG_NOGOOD
    else {
      if (_DEBUG_NOGOOD) {
        std::cout << std::endl;
      }
    }
#endif

  } while (--pathC);
  // p is the last decision, since all atoms above it in the
  // assumption stack have been skipped or expended.
  learnt_clause[0] = NOT(p);

#ifdef _DEBUG_SEARCH
  if (_DEBUG_SEARCH) {
    for (int i = 0; i < level; ++i)
      std::cout << " ";
    std::cout << "learn " << learnt_clause.size << " (";
    print_literal(std::cout, learnt_clause[0]);
    for (unsigned int i = 1; i < learnt_clause.size; ++i) {
      std::cout << " v "; //<< learnt_clause[i];
      print_literal(std::cout, learnt_clause[i]);
    }
    std::cout << " ) " << (backtrack_level < level - 1 ? "-backjump" : "")
              << std::endl;
  }
#endif

  // exit(1);

#ifdef _CHECK_NOGOOD
  store_nogood(learnt_clause);
#endif

  statistics.size_learned += learnt_clause.size;
  statistics.avg_learned_size =
      ((statistics.avg_learned_size * (double)(statistics.num_failures)) +
       (double)(learnt_clause.size)) /
      (double)(statistics.num_failures + 1);

  ++statistics.num_failures;

  if (learnt_clause.size != 1) {

    // if(lit_activity) {
    //   int i=learnt_clause.size;
    //   while(i--) {
    // 	var_activity[UNSIGNED(learnt_clause[i])] +=
    // parameters.activity_increment;
    //   }
    // }

    base->learn(
        learnt_clause,
        (parameters.init_activity ? parameters.activity_increment : 0.0));
    // add_clause( learnt, learnt_clause, stats.learnt_avg_size );
    // reason[UNSIGNED(p)] = base->learnt.back();

    // EXPL
    // base->reason_for[UNSIGNED(p)] = base->learnt.back();

    // base->reason_for[UNSIGNED(p)] = base->learnt.back();
    // reason_for[UNSIGNED(p)] = base;
    // taboo_constraint = base;

    taboo_constraint = (ConstraintImplementation *)(base->learnt.back());
    // reason_for[UNSIGNED(p)].store_reason_for_change(VALUE_EVENT,
    // base->learnt.back());
  } else {
    taboo_constraint = NULL;
  }
  visited.clear();

  // std::cout << visited_literals << std::endl;

  // int real_size = 0;
  // for(int i=0; i<base->learnt.size; ++i) {
  //   real_size += base->learnt[i]->size;
  // }
  // if(real_size != statistics.size_learned) {
  //   std::cout << "discrepancy after learning!!\n" ;
  //   exit(1);
  // }

  // backjump_decision = decision(variables[UNSIGNED(p)], Decision::REMOVAL,
  // SIGN(p));

#ifdef _DEBUG_NOGOOD
  if (_DEBUG_NOGOOD) {
    // for(int i=0; i<level; ++i) std::cout << " ";
    std::cout << "backtrackLevel = " << backtrack_level << "/"
              << (decisions.size) << std::endl;
  }
#endif

  //   while(level>backtrack_level) {
  //     restore();
  //     decisions.pop();
  //   }

  // return decision;

  // exit(1);
}

void Mistral::Solver::forget() {

  // std::cout << lit_activity << " "  << lit_activity[0] << " "  <<
  // lit_activity[1] << std::endl;

  if (base)
    statistics.size_learned -=
        base->forget(parameters.forgetfulness, NULL, NULL);
  //(var_activity, lit_activity);

  // exit(1);

  // int real_size = 0;
  //  for(int i=0; i<base->learnt.size; ++i) {
  //    real_size += base->learnt[i]->size;
  //  }
  //  if(real_size != statistics.size_learned) {
  //    std::cout << "discrepancy after forget!!\n" ;
  //    exit(1);
  //  }
}


void Mistral::Solver::close_propagation() {
  unsigned int vidx;

  // first, resolve the unresolved events, so that 'reason_for' and
  // 'assignment_level' are up to date
  Triplet<int, Event, ConstraintImplementation *> var_evt;
  while (!active_variables.empty()) {
    var_evt = active_variables.pop_front();
    vidx = var_evt.first;
    if (ASSIGNED(var_evt.second) && sequence.contain(variables[vidx])) {
      sequence.remove(variables[vidx]);
      assignment_level[vidx] = level;
      assignment_order[vidx] = assignment_rank;
      ++assignment_rank;

      // std::cout << "REASON FOR x" << vidx << ":";
      // if(var_evt.third)
      // 	std::cout << " c" << var_evt.third->id << std::endl;
      // else
      // 	std::cout << " decision" << std::endl;

      reason_for[vidx] = var_evt.third; //->explain();
    }
  }
}


Mistral::Outcome Mistral::Solver::branch_right() {
  // std::cout << "level :: " << level << std::endl;

  int status = UNKNOWN;

  if (level == search_root)
    status = exhausted(); // objective);

#ifdef _OLD_

#else
  else if (limits_expired()) {

#ifdef _DEBUG_SEARCH
    if (_DEBUG_SEARCH) {
      std::cout << parameters.prefix_comment;
      for (unsigned int k = 0; k <= decisions.size; ++k)
        std::cout << " ";
      std::cout << "limit out!" << std::endl;
    }
#endif

    status = LIMITOUT;
  }
#endif
  else {

    // #ifdef _DEBUG_SEARCH
    //     if(_DEBUG_SEARCH) {
    //       std::cout << parameters.prefix_comment;
    //       for(unsigned int k=0; k<=decisions.size; ++k) std::cout << " ";
    //       //std::cout << "limit fine: " << statistics.num_failures << " < "
    //       << parameters.restart_limit << std::endl;
    //     }
    // #endif

    Mistral::Decision deduction;

    if (parameters.backjump && !culprit.empty()) {

#ifdef _OLD_

#else

      learn_nogood();

#endif

      Literal p = learnt_clause[0];
      deduction =
          Decision(variables[UNSIGNED(p)], Decision::REMOVAL, NOT(SIGN(p)));
    } else {
      backtrack_level = level - 1;
      deduction = decisions.back();
      deduction.invert();
    }

    notify_backtrack();
    restore(backtrack_level);

#ifdef _DEBUG_SEARCH
    if (_DEBUG_SEARCH) {
      std::cout << parameters.prefix_comment;
      for (unsigned int k = 0; k <= decisions.size; ++k)
        std::cout << " ";
      std::cout << "backtrack to lvl " << level << " and deduce " << deduction
                << " (" << statistics.num_filterings << ")" << std::endl;
    }
#endif

#ifdef _MONITOR
    monitor_list.display(std::cout);
    std::cout << std::endl;
// display(std::cout, 2);
// check_constraint_graph_integrity();
#endif


#ifdef _OUTPUT_TIKZ
    std::cout << "node (Refutation" << statistics.num_nodes
              << ") [refutation] {";
    deduction.display_latex(std::cout);
    std::cout << "}\n";
#endif

    // decisions.back(-1).make();
    // decision.make();
    deduction.make();
    // taboo_constraint = NULL;
    // }
  }

  return status;
}


void Mistral::Solver::backjump() {
  int backtrack_level = culprit.get_backtrack_level();
  decisions.size -= (level - backtrack_level);
  restore(backtrack_level);
  Decision decision = culprit.get_decision();

#ifdef _DEBUG_SEARCH
  if (_DEBUG_SEARCH) {
    std::cout << parameters.prefix_comment;
    for (unsigned int k = 0; k <= decisions.size; ++k)
      std::cout << " ";
    std::cout << decision << std::endl;
  }
#endif

  decision.make();
}

void Mistral::Solver::branch_left() {

  save();

  Mistral::Decision decision = heuristic->branch();

  if (decision.var.get_size() >= 20) {
    ++statistics.num_branch_on_large_domains;
  }

  reason_for[decision.var.id()] = NULL;

  decisions.add(decision);

#ifdef _DEBUG_SEARCH
  if (_DEBUG_SEARCH) {
    std::cout << parameters.prefix_comment;
    for (unsigned int k = 0; k <= decisions.size; ++k)
      std::cout << " ";
    std::cout << decision << std::endl;
  }
#endif

#ifdef _OUTPUT_TIKZ
  std::cout << "node (Decision" << statistics.num_nodes << ") [decision] {";
  decision.display_latex(std::cout);
  std::cout << "}\n";
#endif



#ifdef _SAFE
  if (decision.var.is_ground()) {

    std::cerr << "The variable " << decision.var
              << " is ground but still in the sequence!! (abort)" << std::endl;

    exit(1);
  }
#endif

  decision.make();

  ++statistics.num_decisions;
  notify_decision();
}


Mistral::Outcome Mistral::Solver::satisfied() {    
#ifdef _DEBUG_SEARCH
  if (_DEBUG_SEARCH) {
    std::cout << parameters.prefix_comment;
    for (unsigned int k = 0; k <= decisions.size; ++k)
      std::cout << " ";
    std::cout << " SAT! " << std::endl;
  }
#endif

  unsigned int i, j, k;

  if (parameters.checked) {

    /// check the current solution
    Vector<int> tmp_sol;
    Variable *scope;
    Constraint C;
    bool all_assigned;
    int real_arity;

    for (i = 0; i < posted_constraints.size; ++i) {
      all_assigned = true;
      C = constraints[posted_constraints[i]];
      // C.consolidate();

      real_arity = 0;
      k = C.arity();
      scope = C.get_scope();
      for (j = 0; j < k; ++j) {
        if (scope[j].is_ground())
          tmp_sol.add(scope[j].get_value());
        else {
          ++real_arity;
          tmp_sol.add(scope[j].get_min());
          all_assigned = false;
          // break;
        }
      }

      bool consistent = true;

      if (!all_assigned) {

        /// !!! This checks that all values are AC, this is too strong
        // for(j=0; j<k && consistent; ++j) {
        //   if(!scope[j].is_ground()) {
        //     int vali, vnext = scope[j].get_min();
        //     do {
        //       vali = vnext;
        //       if(!C.find_support(j, vali)) consistent = false;
        //       vnext = scope[j].next(vali);
        //     } while( consistent && vali<vnext );
        //   }
        // }

        // #ifdef _DEG_SEARCH
        // 	if(_DEBUG_SEARCH) {
        // 	  std::cout << "c check incomplete assignment of " << C << " ("
        // << real_arity << ")" << std::endl;
        // 	}
        // #endif

        if (real_arity < 5) {
          /// This checks that all bounds are BC
          for (j = 0; j < k && consistent; ++j) {
            if (!scope[j].is_ground()) {
              if (!C.find_bound_support(j, scope[j].get_min()))
                consistent = false;
              else if (!C.find_bound_support(j, scope[j].get_max()))
                consistent = false;
            }
          }
        }

      } else {
        consistent = !C.check(tmp_sol.stack_);
      }

      if (!consistent) {

        std::cerr << "\n@ p" << statistics.num_propagations;
        if (tmp_sol.size < k) {
          std::cerr << ", Error: solution does not satisfy c" << C.id() << ": "
                    << C << tmp_sol << " (backtracking)" << std::endl;
          exit(0);
        } else {
          std::cerr << ", Error: solution does not satisfy c" << C.id() << ": "
                    << C;
          for (j = 0; j < k; ++j) {
            std::cerr << " " << scope[j].get_domain();
          }
          std::cerr << " (backtracking)" << std::endl;
          exit(0);
        }
        if (decisions.empty())
          return UNSAT;
        else if (limits_expired())
          return LIMITOUT;
        else {
          branch_right();
          return UNKNOWN;
        }
      }
      tmp_sol.clear();
    }
  }

  /// store the solution
  for (i = 0; i < variables.size; ++i) {
    last_solution_lb[i] = variables[i].get_min();
    last_solution_ub[i] = variables[i].get_max();

    // std::cout << variables[i] << " := " << last_solution_lb[i] << " " <<
    // variables[i].get_domain() << std::endl;
  }
  // std::cout << std::endl;
  ++statistics.num_solutions;

#ifdef _DEBUG_SEARCH
  if (_DEBUG_SEARCH) {
    std::cout << parameters.prefix_comment;
    for (unsigned int k = 0; k <= decisions.size; ++k)
      std::cout << " ";
    if (objective && !objective->objective.is_void())
      std::cout << "notify solution to goal " << objective << " "
                << objective->objective << " "
                << objective->objective.get_domain() << std::endl;
  }
#endif

  /// notify the objective and return the outcome
  Outcome
      result = //(_exit_on_solution_ ? SAT : objective->notify_solution(this));
      objective->notify_solution(this);

  statistics.objective_value = objective->value();
  // std::cout << statistics.objective_value << std::endl;

  statistics.best_time = get_run_time() - statistics.start_time;
  // std::cout << " c best_time  " <<  statistics.best_time << std::endl;

  for (i = 0; i < solution_triggers.size; ++i) {
    solution_triggers[i]->notify_solution();
  }

  // if(statistics.objective_value == 933)
  //   exit(1);

#ifdef _DEBUG_SEARCH
  if (_DEBUG_SEARCH) {
    std::cout << parameters.prefix_comment;
    for (unsigned int k = 0; k <= decisions.size; ++k)
      std::cout << " ";
    std::cout << "=> " << outcome2str(result) << std::endl;
  }
#endif

  return result;

  // return SAT;
}


Mistral::Outcome Mistral::Solver::exhausted() {    
#ifdef _DEBUG_SEARCH
  if (_DEBUG_SEARCH) {
    std::cout << "c UNSAT!" << std::endl;
  }
#endif

  Outcome value = UNSAT;
  if (statistics.num_solutions)
    value = objective->notify_exhausted();
  return value;
}

bool Mistral::Solver::limits_expired() {

#ifdef _DEBUG_SEARCH
  if (_DEBUG_SEARCH) {
    if (parameters.limit &&
        ((parameters.time_limit > 0.0 &&
          (get_run_time() - statistics.start_time) > parameters.time_limit) ||
         (parameters.node_limit > 0 &&
          (statistics.num_nodes > parameters.node_limit)) ||
         (parameters.fail_limit > 0 &&
          (statistics.num_failures > parameters.fail_limit)) ||
         (parameters.restart_limit > 0 &&
          (statistics.num_failures > parameters.restart_limit)) ||
         (parameters.propagation_limit > 0 &&
          (statistics.num_propagations > parameters.propagation_limit)) // ||
         // (parameters.backtrack_limit > 0 && (statistics.num_backtracks >
         // parameters.backtrack_limit))
         ))
      std::cout << "c LIMIT REACHED, RESTART!" << std::endl;
  }
#endif

#ifdef _PARALLEL
if (*solution_found_elsewhere){
  // std::cout << "c LIMIT REACHED, BECAUSE solution_found_elsewhere !" <<
  // std::endl;
  return true;
}
#endif

// std::cout << parameters.time_limit << " " << (get_run_time() -
// statistics.start_time) << " " << parameters.time_limit << std::endl;

auto res =
    (parameters.limit &&
     ((parameters.time_limit > 0.0 &&
       (get_run_time() - statistics.start_time) > parameters.time_limit) ||
      (parameters.node_limit > 0 &&
       (statistics.num_nodes > parameters.node_limit)) ||
      (parameters.fail_limit > 0 &&
       (statistics.num_failures > parameters.fail_limit)) ||
      (parameters.restart_limit > 0 &&
       (statistics.num_failures > parameters.restart_limit)) ||
      (parameters.propagation_limit > 0 &&
       (statistics.num_propagations > parameters.propagation_limit)) // ||
      // (parameters.backtrack_limit > 0 && (statistics.num_backtracks >
      // parameters.backtrack_limit))
      ));

// if(res) {
// 	if((parameters.time_limit > 0.0 && (get_run_time() -
// statistics.start_time) > parameters.time_limit)) 		std::cout <<
// "TIME!\n"; 	if(parameters.node_limit > 0 && (statistics.num_nodes >
// parameters.node_limit)) 		std::cout << "NODE!\n";
// if(parameters.fail_limit > 0
// && (statistics.num_failures > parameters.fail_limit)) 		std::cout <<
// "FAIL!\n"; 	if(parameters.restart_limit > 0 && (statistics.num_failures >
// parameters.restart_limit)) 		std::cout << "RESTART!\n";
// 	if(parameters.propagation_limit > 0 && (statistics.num_propagations >
// parameters.propagation_limit)) 		std::cout << "PROPAG!\n";
// }

return res;
}

// void Mistral::Search::init_search(Vector< Variable >& seq, VarOrdering *h, RestartPolicy *p) {
//   for(unsigned int i=0; i<seq.size; ++i) {
//     std::cout << "insert " << seq[i] << std::endl;
//     sequence.insert(seq[i]);    
//   }
//   if(heuristic) delete heuristic;
//   if(policy) delete policy;
//   heuristic = (h ? h : new NoOrder(sequence));
//   policy = p;
// }


// std::string Mistral::toString(const Mistral::Solver& x) {
//   return x.getString();
// }

/**
Searches for an assignment of the variables in sequences

Stops when the objecive is satisfied, or when the search tree is exhausted
- SATISFACTION OBJECTIVE: returns 'SAT' when a solution is found, 'UNSAT' if the
search tree is exhausted, and 'UNKNOWN' otherwise
- OPTIMIZATION OBJECTIVE: returns 'OPT' when the search tree is exhausted and
'UNKNOWN' otherwise
- ENUMERATION  OBJECTIVE: returns 'ALL' when the search tree is exhausted and
'UNKNOWN' otherwise [TODO]

The argument 'root' controls how deep we can backtrack. i.e., the root of the
search tree. It is possible to start the search with the decision stack and
trail non empty, and backtrack on them (if 'root' is set to something lower than
their size), or not.
*/
Mistral::Outcome Mistral::Solver::chronological_dfs(const int _root) 
{
  search_root = _root;

#ifdef _OUTPUT_TIKZ
  int tikz_level = 0;
  Vector<int> deficit;
  deficit.add(0);
  // int last = tikz_level;
  if (!_OUTPUT_TIKZ) {
    std::cout << "\\node (Root) [refutation] {$\\emptyset$}\n";
  }
#endif

  int status = UNKNOWN;
  while (status == UNKNOWN) {

    if (propagate()) {

#ifdef _OUTPUT_TIKZ
      if (_OUTPUT_TIKZ) {
        if (tikz_level > 0)
          std::cout << "child{\n";
        else
          std::cout << "\\";
        std::cout << "node (state" << statistics.num_filterings
                  << ") [state] {\n";
        std::cout << "\\begin{tabular}{ll}\n";
        for (unsigned int i = 0; i < variables.size; ++i) {
          std::cout << "${\\cal D}(x_" << variables[i].id() << ")$ & $"
                    << variables[i].get_domain(true) << "$ \\\\\n";
        }
        std::cout << "\\end{tabular}\n}\n";
      }
#endif

#ifdef _MONITOR
      monitor_list.display(std::cout);
      std::cout << std::endl;
// display(std::cout, 2);
// check_constraint_graph_integrity();
#endif

#ifdef _OUTPUT_TIKZ
      deficit.add(0);
      ++tikz_level;
      // last = tikz_level;
      std::cout << "child{ [sibling distance="
                << (tikz_level > 5 ? 4 : (1 << (7 - tikz_level))) << "mm]\n";
#endif

      ++statistics.num_nodes;
      if (sequence.empty())
        status = satisfied();
      else
        branch_left();

    } else {

#ifdef _OUTPUT_TIKZ
      std::cout << "child{ \nnode[fail] (Fail" << statistics.num_failures
                << ") {\\bf fail}\n}\n";

      // std::cout << "LEVEL=" << level << " deficit=" << deficit << std::endl;

      int open = deficit.pop();
      while (open--) {
        --tikz_level;

        if (_OUTPUT_TIKZ) {
          std::cout << "}\n";
        }

        std::cout << "}\n";
      }
      ++deficit[level - 1];

      // std::cout << "LEVEL=" << level << " deficit=" << deficit << std::endl;

      if (level != search_root) {
        // if(_OUTPUT_TIKZ) {
        //   std::cout << "}\n";
        // }
        std::cout << "}\nchild{ [sibling distance="
                  << (tikz_level > 5 ? 4 : (1 << (7 - tikz_level))) << "mm]\n";
      }
#endif

#ifdef _OLD_
      if (parameters.backjump)
        learn_nogood();

      if (limits_expired()) {
        status = // interrupted();
            LIMITOUT;
      } else
        status = branch_right();
#else
      status = branch_right();
#endif
    }
  }

  // std::cout << outcome2str(status) << std::endl;

  //  std::cout << deficit << std::endl;

#ifdef _OUTPUT_TIKZ
  if (deficit.size) {
    int open = deficit.pop();
    if (open > 0)
      while (--open) {
        if (_OUTPUT_TIKZ) {
          std::cout << "}\n";
        }
        std::cout << "}\n";
      }
  }
  std::cout << ";\n";
#endif

  return status;
}


Mistral::Outcome Mistral::Solver::interrupted() {
#ifdef _DEBUG_SEARCH
  if (_DEBUG_SEARCH) {
    std::cout << "c UNSAT!" << std::endl;
  }
#endif

  Outcome value = LIMITOUT;
  // if(objective && statistics.num_solutions) value = SAT;
  return value;
}


Mistral::Outcome Mistral::Solver::conflict_directed_backjump()
{
  int status = UNKNOWN;
  while (status == UNKNOWN) {
    if (propagate()) {
      if (sequence.empty())
        status = satisfied();
      else
        branch_left();
    } else {
      if (decisions.empty())
        status = UNSAT;
      else if (limits_expired())
        status = LIMITOUT;
      else
        backjump();
    }
  }

  return status;
}

// double *Mistral::Solver::get_literal_activity() {
//   return lit_activity.stack_;
// }

std::ostream& Mistral::operator<< (std::ostream& os, Mistral::Solution& x) {
  return x.display(os);
}
std::ostream& Mistral::operator<< (std::ostream& os, Mistral::Solution* x) {
  return x->display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, Mistral::Solver& x) {
  return x.display(os);
}
std::ostream& Mistral::operator<< (std::ostream& os, Mistral::Solver* x) {
  return x->display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, Mistral::ConstraintQueue& x) {
  return x.display(os);
}
std::ostream& Mistral::operator<< (std::ostream& os, Mistral::ConstraintQueue* x) {
  return x->display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::SolverStatistics& x) {
  return x.display(os);
}
std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::SolverStatistics* x) {
  return x->display(os);
}


void Mistral::SearchMonitor::add(Variable x) {
  sequence.add(0);
  sequence.add(x.id());

  sequence_type.add(0);
  sequence_var.add(x);
}
void Mistral::SearchMonitor::add(Constraint x) {
  sequence.add(1);
  sequence.add(x.id());

  sequence_type.add(1);
  sequence_con.add(x);
}
// void Mistral::SearchMonitor::add(std::string& x) {
//   sequence.add(2);
//   sequence.add(strs.size);
//   strs.push_back(x);
// }
void Mistral::SearchMonitor::add(const char* x) {
  sequence.add(2);
  sequence.add(strs.size());

  strs.push_back(x);
  sequence_type.add(2);
  // strs.push_back(x);
}

std::ostream& Mistral::SearchMonitor::display( std::ostream& os ) const {
  // for(unsigned int i=0; i<sequence.size; i+=2) {
  //   if(sequence[i] == 0) {
  //     os << solver->variables[sequence[i+1]].get_domain();
  //   } else if(sequence[i] == 1) {
  //     os << solver->constraints[sequence[i+1]];
  //   } else {
  //     os << strs[sequence[i+1]];
  //   }
  // }
  // return os;

  int v = 0;
  int c = 0;
  int s = 0;
  for (unsigned int i = 0; i < sequence_type.size; ++i) {
    if (sequence_type[i] == 0) {
      if (sequence_var[v].is_ground() || !(sequence_var[v].is_boolean()))
        os << sequence_var[v].get_domain();
      else
        os << ".";
      ++v;
    } else if (sequence_type[i] == 1) {
      os << sequence_con[c++];
    } else {
      os << strs[s++];
    }
  }
  return os;
}

Mistral::SearchMonitor& Mistral::operator<< (Mistral::SearchMonitor& os, VarArray& x) {
  os.add("(");
  os.add(x[0]);
  for (unsigned int i = 1; i < x.size; ++i) {
    // os.add(" ");
    os.add(x[i]);
  }
  os.add(")");
  return os;
}
Mistral::SearchMonitor& Mistral::operator<< (Mistral::SearchMonitor& os, Variable& x) {
  os.add(x);
  return os;
}
Mistral::SearchMonitor& Mistral::operator<< (Mistral::SearchMonitor& os, Constraint& x) {
  os.add(x);
  return os;
}
Mistral::SearchMonitor& Mistral::operator<< (Mistral::SearchMonitor& os, const char* x) {
  os.add(x);
  return os;
}
// Mistral::SearchMonitor& Mistral::operator<< (Mistral::SearchMonitor& os, const char* x) {
//   std::string s(x);
//   os.add(s);
//   return os;
// }
Mistral::SearchMonitor& Mistral::operator<< (Mistral::SearchMonitor& os, const int x) {
  // std::ostringstream o_aux;
  // o_aux << x;
  // const char *y = ;
  //  os.add((o_aux.str().c_str()));
  // o_aux.close();

  // std::string number = boost::lexical_cast(x);
  //  char buf[255];
  //  sprintf(buf,"%d",x);

  os << " "; // int2str(x).c_str() ;

  // os << number.c_str() ;
  return os;
}



Mistral::BranchingHeuristic *Mistral::Solver::heuristic_factory(std::string var_ordering, std::string branching, const int randomness) {

  BranchingHeuristic *heu = NULL;

  if (var_ordering.find("MBA") == 0) {
    var_ordering.erase(0, 3);

    if (var_ordering == "wdeg" || var_ordering == "WDEG" ||
        var_ordering == "Weighted Degree" || var_ordering == "WeightedDegree") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, FailureCountManager>,
              AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, FailureCountManager>,
              AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, FailureCountManager>,
              AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, FailureCountManager>,
              AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, FailureCountManager>,
              AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, FailureCountManager>,
              MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, FailureCountManager>,
              MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, FailureCountManager>,
              MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, FailureCountManager>,
              MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, FailureCountManager>,
              MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, FailureCountManager>,
              MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, FailureCountManager>,
              MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, FailureCountManager>,
              MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, FailureCountManager>,
              MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, FailureCountManager>,
              MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, FailureCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, FailureCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, FailureCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, FailureCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, FailureCountManager>,
              HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, FailureCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, FailureCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, FailureCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, FailureCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, FailureCountManager>,
              RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, FailureCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, FailureCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, FailureCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, FailureCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, FailureCountManager>,
              RandomMinMax>(this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, FailureCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, FailureCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, FailureCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, FailureCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, FailureCountManager>,
              MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, FailureCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, FailureCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, FailureCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, FailureCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, FailureCountManager>,
              Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, FailureCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, FailureCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, FailureCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, FailureCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, FailureCountManager>,
              Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, FailureCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, FailureCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, FailureCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, FailureCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, FailureCountManager>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, FailureCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, FailureCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, FailureCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, FailureCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, FailureCountManager>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "dom/wdeg" || var_ordering == "dwdeg" ||
        var_ordering == "DWDEG" ||
        var_ordering == "Domain Over Weighted Degree" ||
        var_ordering == "DomainOverWeightedDegree") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         FailureCountManager>,
              AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         FailureCountManager>,
              AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         FailureCountManager>,
              AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         FailureCountManager>,
              AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         FailureCountManager>,
              AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         FailureCountManager>,
              MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         FailureCountManager>,
              MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         FailureCountManager>,
              MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         FailureCountManager>,
              MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         FailureCountManager>,
              MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         FailureCountManager>,
              MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         FailureCountManager>,
              MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         FailureCountManager>,
              MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         FailureCountManager>,
              MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         FailureCountManager>,
              MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         FailureCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         FailureCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         FailureCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         FailureCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         FailureCountManager>,
              HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         FailureCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         FailureCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         FailureCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         FailureCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         FailureCountManager>,
              RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         FailureCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         FailureCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         FailureCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         FailureCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         FailureCountManager>,
              RandomMinMax>(this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         FailureCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         FailureCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         FailureCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         FailureCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         FailureCountManager>,
              MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         FailureCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         FailureCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         FailureCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         FailureCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         FailureCountManager>,
              Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         FailureCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         FailureCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         FailureCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         FailureCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         FailureCountManager>,
              Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         FailureCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         FailureCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         FailureCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         FailureCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         FailureCountManager>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         FailureCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         FailureCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         FailureCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         FailureCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         FailureCountManager>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "gwdeg" || var_ordering == "GWDEG" ||
        var_ordering == "Global Weighted Degree" ||
        var_ordering == "GlobalWeightedDegree") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, ConflictCountManager>,
              AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, ConflictCountManager>,
              MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, ConflictCountManager>,
              MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, ConflictCountManager>,
              HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, ConflictCountManager>,
              RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, ConflictCountManager>,
              RandomMinMax>(this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, ConflictCountManager>,
              MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, ConflictCountManager>,
              Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, ConflictCountManager>,
              Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 1, ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 2, ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 3, ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 4, ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MaxWeight>, 5, ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "dom/gwdeg" || var_ordering == "dgwdeg" ||
        var_ordering == "DGWDEG" ||
        var_ordering == "Domain Over Global Weighted Degree" ||
        var_ordering == "DomainOverGlobalWeightedDegree") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         ConflictCountManager>,
              AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         ConflictCountManager>,
              MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         ConflictCountManager>,
              MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         ConflictCountManager>,
              HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         ConflictCountManager>,
              RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         ConflictCountManager>,
              RandomMinMax>(this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         ConflictCountManager>,
              MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         ConflictCountManager>,
              Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         ConflictCountManager>,
              Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "Impact" || var_ordering == "impact") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, ImpactManager>,
              AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, ImpactManager>,
              AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, ImpactManager>,
              AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, ImpactManager>,
              AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, ImpactManager>,
              AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, ImpactManager>,
              MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, ImpactManager>,
              MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, ImpactManager>,
              MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, ImpactManager>,
              MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, ImpactManager>,
              MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, ImpactManager>,
              MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, ImpactManager>,
              MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, ImpactManager>,
              MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, ImpactManager>,
              MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, ImpactManager>,
              MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, ImpactManager>,
              HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, ImpactManager>,
              HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, ImpactManager>,
              HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, ImpactManager>,
              HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, ImpactManager>,
              HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, ImpactManager>,
              RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, ImpactManager>,
              RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, ImpactManager>,
              RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, ImpactManager>,
              RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, ImpactManager>,
              RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, ImpactManager>,
              RandomMinMax>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, ImpactManager>,
              RandomMinMax>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, ImpactManager>,
              RandomMinMax>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, ImpactManager>,
              RandomMinMax>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, ImpactManager>,
              RandomMinMax>(this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, ImpactManager>,
              MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, ImpactManager>,
              MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, ImpactManager>,
              MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, ImpactManager>,
              MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, ImpactManager>,
              MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, ImpactManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, ImpactManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, ImpactManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, ImpactManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, ImpactManager>,
              Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, ImpactManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, ImpactManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, ImpactManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, ImpactManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, ImpactManager>,
              Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, ImpactManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, ImpactManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, ImpactManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, ImpactManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, ImpactManager>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, ImpactManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, ImpactManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, ImpactManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, ImpactManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, ImpactManager>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, ImpactManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, ImpactManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, ImpactManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, ImpactManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, ImpactManager>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, ImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, ImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, ImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, ImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, ImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "RIBS" || var_ordering == "ribs" ||
        var_ordering == "RealImpact" || var_ordering == "realimpact") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, RealImpactManager>,
              AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, RealImpactManager>,
              AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, RealImpactManager>,
              AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, RealImpactManager>,
              AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, RealImpactManager>,
              AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, RealImpactManager>,
              MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, RealImpactManager>,
              MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, RealImpactManager>,
              MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, RealImpactManager>,
              MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, RealImpactManager>,
              MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, RealImpactManager>,
              MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, RealImpactManager>,
              MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, RealImpactManager>,
              MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, RealImpactManager>,
              MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, RealImpactManager>,
              MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, RealImpactManager>,
              HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, RealImpactManager>,
              HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, RealImpactManager>,
              HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, RealImpactManager>,
              HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, RealImpactManager>,
              HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, RealImpactManager>,
              RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, RealImpactManager>,
              RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, RealImpactManager>,
              RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, RealImpactManager>,
              RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, RealImpactManager>,
              RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, RealImpactManager>,
              RandomMinMax>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, RealImpactManager>,
              RandomMinMax>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, RealImpactManager>,
              RandomMinMax>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, RealImpactManager>,
              RandomMinMax>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, RealImpactManager>,
              RandomMinMax>(this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, RealImpactManager>,
              MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, RealImpactManager>,
              MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, RealImpactManager>,
              MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, RealImpactManager>,
              MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, RealImpactManager>,
              MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, RealImpactManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, RealImpactManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, RealImpactManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, RealImpactManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, RealImpactManager>,
              Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, RealImpactManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, RealImpactManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, RealImpactManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, RealImpactManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, RealImpactManager>,
              Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, RealImpactManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, RealImpactManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, RealImpactManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, RealImpactManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, RealImpactManager>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, RealImpactManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, RealImpactManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, RealImpactManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, RealImpactManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, RealImpactManager>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, RealImpactManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, RealImpactManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, RealImpactManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, RealImpactManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, RealImpactManager>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 1, RealImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 2, RealImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 3, RealImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 4, RealImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinWeight>, 5, RealImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    // if(var_ordering == "IBS" || var_ordering == "ibs" || var_ordering ==
    // "Impact" || var_ordering == "impact" ) {
    //   if( branching == "No" || branching == "no" || branching == "Any" ||
    //   branching == "any" ) {
    //     if( randomness <= 1 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 1, RealImpactManager >, AnyValue > (this);
    //     }
    //     else if( randomness <= 2 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 2, RealImpactManager >, AnyValue > (this);
    //     }
    //     else if( randomness <= 3 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 3, RealImpactManager >, AnyValue > (this);
    //     }
    //     else if( randomness <= 4 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 4, RealImpactManager >, AnyValue > (this);
    //     }
    //     else if( randomness <= 5 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 5, RealImpactManager >, AnyValue > (this);
    //     }
    //   }
    //   if( branching == "Lex" || branching == "lex" || branching == "minvalue"
    //   || branching == "min value" || branching == "MinValue" || branching ==
    //   "Min Value" || branching == "lexicographic" || branching ==
    //   "Lexicographic" ) {
    //     if( randomness <= 1 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 1, RealImpactManager >, MinValue > (this);
    //     }
    //     else if( randomness <= 2 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 2, RealImpactManager >, MinValue > (this);
    //     }
    //     else if( randomness <= 3 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 3, RealImpactManager >, MinValue > (this);
    //     }
    //     else if( randomness <= 4 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 4, RealImpactManager >, MinValue > (this);
    //     }
    //     else if( randomness <= 5 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 5, RealImpactManager >, MinValue > (this);
    //     }
    //   }
    //   if( branching == "AntiLex" || branching == "antilex" || branching ==
    //   "maxvalue" || branching == "max value" || branching == "MaxValue" ||
    //   branching == "Max Value" || branching == "antilexicographic" ||
    //   branching == "Antilexicographic" ) {
    //     if( randomness <= 1 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 1, RealImpactManager >, MaxValue > (this);
    //     }
    //     else if( randomness <= 2 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 2, RealImpactManager >, MaxValue > (this);
    //     }
    //     else if( randomness <= 3 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 3, RealImpactManager >, MaxValue > (this);
    //     }
    //     else if( randomness <= 4 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 4, RealImpactManager >, MaxValue > (this);
    //     }
    //     else if( randomness <= 5 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 5, RealImpactManager >, MaxValue > (this);
    //     }
    //   }
    //   if( branching == "HalfSplit" || branching == "halfsplit" ) {
    //     if( randomness <= 1 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 1, RealImpactManager >, HalfSplit > (this);
    //     }
    //     else if( randomness <= 2 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 2, RealImpactManager >, HalfSplit > (this);
    //     }
    //     else if( randomness <= 3 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 3, RealImpactManager >, HalfSplit > (this);
    //     }
    //     else if( randomness <= 4 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 4, RealImpactManager >, HalfSplit > (this);
    //     }
    //     else if( randomness <= 5 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 5, RealImpactManager >, HalfSplit > (this);
    //     }
    //   }
    //   if( branching == "RandomSplit" || branching == "RandSplit" || branching
    //   == "randomsplit" || branching == "randsplit" ) {
    //     if( randomness <= 1 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 1, RealImpactManager >, RandomSplit >
    //       (this);
    //     }
    //     else if( randomness <= 2 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 2, RealImpactManager >, RandomSplit >
    //       (this);
    //     }
    //     else if( randomness <= 3 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 3, RealImpactManager >, RandomSplit >
    //       (this);
    //     }
    //     else if( randomness <= 4 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 4, RealImpactManager >, RandomSplit >
    //       (this);
    //     }
    //     else if( randomness <= 5 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 5, RealImpactManager >, RandomSplit >
    //       (this);
    //     }
    //   }
    //   if( branching == "RandomMinMax" || branching == "randomminmax" ||
    //   branching == "RandMinMax" || branching == "randminmax" ) {
    //     if( randomness <= 1 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 1, RealImpactManager >, RandomMinMax >
    //       (this);
    //     }
    //     else if( randomness <= 2 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 2, RealImpactManager >, RandomMinMax >
    //       (this);
    //     }
    //     else if( randomness <= 3 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 3, RealImpactManager >, RandomMinMax >
    //       (this);
    //     }
    //     else if( randomness <= 4 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 4, RealImpactManager >, RandomMinMax >
    //       (this);
    //     }
    //     else if( randomness <= 5 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 5, RealImpactManager >, RandomMinMax >
    //       (this);
    //     }
    //   }
    //   if( branching == "minweight" || branching == "MinWeight" ) {
    //     if( randomness <= 1 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 1, RealImpactManager >, MinWeightValue >
    //       (this);
    //     }
    //     else if( randomness <= 2 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 2, RealImpactManager >, MinWeightValue >
    //       (this);
    //     }
    //     else if( randomness <= 3 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 3, RealImpactManager >, MinWeightValue >
    //       (this);
    //     }
    //     else if( randomness <= 4 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 4, RealImpactManager >, MinWeightValue >
    //       (this);
    //     }
    //     else if( randomness <= 5 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 5, RealImpactManager >, MinWeightValue >
    //       (this);
    //     }
    //   }
    //   if( branching == "MinVal+Guided" || branching == "minval+guided" ) {
    //     if( randomness <= 1 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 1, RealImpactManager >, Guided< MinValue >
    //       > (this);
    //     }
    //     else if( randomness <= 2 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 2, RealImpactManager >, Guided< MinValue >
    //       > (this);
    //     }
    //     else if( randomness <= 3 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 3, RealImpactManager >, Guided< MinValue >
    //       > (this);
    //     }
    //     else if( randomness <= 4 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 4, RealImpactManager >, Guided< MinValue >
    //       > (this);
    //     }
    //     else if( randomness <= 5 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 5, RealImpactManager >, Guided< MinValue >
    //       > (this);
    //     }
    //   }
    //   if( branching == "MaxVal+Guided" || branching == "maxval+guided" ) {
    //     if( randomness <= 1 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 1, RealImpactManager >, Guided< MaxValue >
    //       > (this);
    //     }
    //     else if( randomness <= 2 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 2, RealImpactManager >, Guided< MaxValue >
    //       > (this);
    //     }
    //     else if( randomness <= 3 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 3, RealImpactManager >, Guided< MaxValue >
    //       > (this);
    //     }
    //     else if( randomness <= 4 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 4, RealImpactManager >, Guided< MaxValue >
    //       > (this);
    //     }
    //     else if( randomness <= 5 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 5, RealImpactManager >, Guided< MaxValue >
    //       > (this);
    //     }
    //   }
    //   if( branching == "MinWeight+Guided" || branching == "minweight+guided"
    //   ) {
    //     if( randomness <= 1 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 1, RealImpactManager >, Guided<
    //       MinWeightValue > > (this);
    //     }
    //     else if( randomness <= 2 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 2, RealImpactManager >, Guided<
    //       MinWeightValue > > (this);
    //     }
    //     else if( randomness <= 3 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 3, RealImpactManager >, Guided<
    //       MinWeightValue > > (this);
    //     }
    //     else if( randomness <= 4 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 4, RealImpactManager >, Guided<
    //       MinWeightValue > > (this);
    //     }
    //     else if( randomness <= 5 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 5, RealImpactManager >, Guided<
    //       MinWeightValue > > (this);
    //     }
    //   }
    //   if( branching == "MaxWeightVal+Guided" || branching ==
    //   "maxweightval+guided" ) {
    //     if( randomness <= 1 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 1, RealImpactManager >, Guided<
    //       MaxWeightValue > > (this);
    //     }
    //     else if( randomness <= 2 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 2, RealImpactManager >, Guided<
    //       MaxWeightValue > > (this);
    //     }
    //     else if( randomness <= 3 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 3, RealImpactManager >, Guided<
    //       MaxWeightValue > > (this);
    //     }
    //     else if( randomness <= 4 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 4, RealImpactManager >, Guided<
    //       MaxWeightValue > > (this);
    //     }
    //     else if( randomness <= 5 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 5, RealImpactManager >, Guided<
    //       MaxWeightValue > > (this);
    //     }
    //   }
    //   if( branching == "Random+Guided" || branching == "random+guided" ) {
    //     if( randomness <= 1 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 1, RealImpactManager >, Guided<
    //       RandomMinMax > > (this);
    //     }
    //     else if( randomness <= 2 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 2, RealImpactManager >, Guided<
    //       RandomMinMax > > (this);
    //     }
    //     else if( randomness <= 3 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 3, RealImpactManager >, Guided<
    //       RandomMinMax > > (this);
    //     }
    //     else if( randomness <= 4 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 4, RealImpactManager >, Guided<
    //       RandomMinMax > > (this);
    //     }
    //     else if( randomness <= 5 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 5, RealImpactManager >, Guided<
    //       RandomMinMax > > (this);
    //     }
    //   }
    //   if( branching == "Adpated" || branching == "adapted" ) {
    //     if( randomness <= 1 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 1, RealImpactManager >, ConditionalOnSize<
    //       GuidedSplit< HalfSplit >, Guided< MinValue > > > (this);
    //     }
    //     else if( randomness <= 2 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 2, RealImpactManager >, ConditionalOnSize<
    //       GuidedSplit< HalfSplit >, Guided< MinValue > > > (this);
    //     }
    //     else if( randomness <= 3 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 3, RealImpactManager >, ConditionalOnSize<
    //       GuidedSplit< HalfSplit >, Guided< MinValue > > > (this);
    //     }
    //     else if( randomness <= 4 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 4, RealImpactManager >, ConditionalOnSize<
    //       GuidedSplit< HalfSplit >, Guided< MinValue > > > (this);
    //     }
    //     else if( randomness <= 5 ) {
    //       heu = new GenericHeuristic < GenericDVO < MultiArmedBandit<
    //       MinDomainTimesWeight >, 5, RealImpactManager >, ConditionalOnSize<
    //       GuidedSplit< HalfSplit >, Guided< MinValue > > > (this);
    //     }
    //   }
    // }
    if (var_ordering == "ABS" || var_ordering == "abs" ||
        var_ordering == "activity based search" ||
        var_ordering == "Activity Based Search" ||
        var_ordering == "dom/activity") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         PruningCountManager>,
              AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         PruningCountManager>,
              AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         PruningCountManager>,
              AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         PruningCountManager>,
              AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         PruningCountManager>,
              AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         PruningCountManager>,
              MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         PruningCountManager>,
              MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         PruningCountManager>,
              MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         PruningCountManager>,
              MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         PruningCountManager>,
              MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         PruningCountManager>,
              MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         PruningCountManager>,
              MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         PruningCountManager>,
              MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         PruningCountManager>,
              MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         PruningCountManager>,
              MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         PruningCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         PruningCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         PruningCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         PruningCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         PruningCountManager>,
              HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         PruningCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         PruningCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         PruningCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         PruningCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         PruningCountManager>,
              RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         PruningCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         PruningCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         PruningCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         PruningCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         PruningCountManager>,
              RandomMinMax>(this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         PruningCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         PruningCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         PruningCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         PruningCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         PruningCountManager>,
              MinWeightValue>(this);
        }
      }
      if (branching == "maxweight" || branching == "MaxWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         PruningCountManager>,
              MaxWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         PruningCountManager>,
              MaxWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         PruningCountManager>,
              MaxWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         PruningCountManager>,
              MaxWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         PruningCountManager>,
              MaxWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         PruningCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         PruningCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         PruningCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         PruningCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         PruningCountManager>,
              Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         PruningCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         PruningCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         PruningCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         PruningCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         PruningCountManager>,
              Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         PruningCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         PruningCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         PruningCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         PruningCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         PruningCountManager>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         PruningCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         PruningCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         PruningCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         PruningCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         PruningCountManager>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         PruningCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         PruningCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         PruningCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         PruningCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         PruningCountManager>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 1,
                         PruningCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 2,
                         PruningCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 3,
                         PruningCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 4,
                         PruningCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainOverWeight>, 5,
                         PruningCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "mindomain" || var_ordering == "MinDomain" ||
        var_ordering == "min domain" || var_ordering == "minimum domain") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, AnyValue>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, AnyValue>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, AnyValue>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, AnyValue>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, AnyValue>(
              this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MinValue>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MinValue>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MinValue>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MinValue>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MinValue>(
              this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MaxValue>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MaxValue>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MaxValue>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MaxValue>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MaxValue>(
              this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, HalfSplit>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, HalfSplit>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, HalfSplit>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, HalfSplit>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, HalfSplit>(
              this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomSplit>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomSplit>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomSplit>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomSplit>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomSplit>(
              this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomMinMax>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomMinMax>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomMinMax>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomMinMax>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomMinMax>(
              this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "Neighbor" || var_ordering == "Neighbor" ||
        var_ordering == "Neighbour" || var_ordering == "neighbour") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              RandomMinMax>(this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "mindom->maxdeg" ||
        var_ordering == "MinDomainMaxDegree") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     RandomMinMax>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     RandomMinMax>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     RandomMinMax>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     RandomMinMax>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     RandomMinMax>(this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MultiArmedBandit<
                                         LexCombination<MinDomain, MaxDegree>>>,
                                     Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<
                  MultiArmedBandit<LexCombination<MinDomain, MaxDegree>>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<
                  MultiArmedBandit<LexCombination<MinDomain, MaxDegree>>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<
                  MultiArmedBandit<LexCombination<MinDomain, MaxDegree>>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<
                  MultiArmedBandit<LexCombination<MinDomain, MaxDegree>>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<
                  MultiArmedBandit<LexCombination<MinDomain, MaxDegree>>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "mindomain" || var_ordering == "MinDomain" ||
        var_ordering == "min domain" || var_ordering == "minimum domain" ||
        var_ordering == "MinDomnMaxDeg") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, AnyValue>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, AnyValue>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, AnyValue>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, AnyValue>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, AnyValue>(
              this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MinValue>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MinValue>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MinValue>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MinValue>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MinValue>(
              this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MaxValue>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MaxValue>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MaxValue>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MaxValue>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, MaxValue>(
              this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, HalfSplit>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, HalfSplit>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, HalfSplit>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, HalfSplit>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, HalfSplit>(
              this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomSplit>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomSplit>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomSplit>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomSplit>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomSplit>(
              this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomMinMax>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomMinMax>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomMinMax>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomMinMax>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>, RandomMinMax>(
              this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MultiArmedBandit<MinDomainTimesWeight>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }

    // MBA ICI
  } else {

    if (var_ordering == "wdeg" || var_ordering == "WDEG" ||
        var_ordering == "Weighted Degree" || var_ordering == "WeightedDegree") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, FailureCountManager>, AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, FailureCountManager>, AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, FailureCountManager>, AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, FailureCountManager>, AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, FailureCountManager>, AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, FailureCountManager>, MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, FailureCountManager>, MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, FailureCountManager>, MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, FailureCountManager>, MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, FailureCountManager>, MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, FailureCountManager>, MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, FailureCountManager>, MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, FailureCountManager>, MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, FailureCountManager>, MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, FailureCountManager>, MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, FailureCountManager>, HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, FailureCountManager>, HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, FailureCountManager>, HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, FailureCountManager>, HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, FailureCountManager>, HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, FailureCountManager>, RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, FailureCountManager>, RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, FailureCountManager>, RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, FailureCountManager>, RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, FailureCountManager>, RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, FailureCountManager>, RandomMinMax>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, FailureCountManager>, RandomMinMax>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, FailureCountManager>, RandomMinMax>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, FailureCountManager>, RandomMinMax>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, FailureCountManager>, RandomMinMax>(
              this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, FailureCountManager>, MinWeightValue>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, FailureCountManager>, MinWeightValue>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, FailureCountManager>, MinWeightValue>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, FailureCountManager>, MinWeightValue>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, FailureCountManager>, MinWeightValue>(
              this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, FailureCountManager>, Guided<MinValue>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, FailureCountManager>, Guided<MinValue>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, FailureCountManager>, Guided<MinValue>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, FailureCountManager>, Guided<MinValue>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, FailureCountManager>, Guided<MinValue>>(
              this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, FailureCountManager>, Guided<MaxValue>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, FailureCountManager>, Guided<MaxValue>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, FailureCountManager>, Guided<MaxValue>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, FailureCountManager>, Guided<MaxValue>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, FailureCountManager>, Guided<MaxValue>>(
              this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, FailureCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, FailureCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, FailureCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, FailureCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, FailureCountManager>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, FailureCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, FailureCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, FailureCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, FailureCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, FailureCountManager>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "dom/wdeg" || var_ordering == "dwdeg" ||
        var_ordering == "DWDEG" ||
        var_ordering == "Domain Over Weighted Degree" ||
        var_ordering == "DomainOverWeightedDegree") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, FailureCountManager>,
              AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, FailureCountManager>,
              AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, FailureCountManager>,
              AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, FailureCountManager>,
              AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, FailureCountManager>,
              AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, FailureCountManager>,
              MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, FailureCountManager>,
              MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, FailureCountManager>,
              MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, FailureCountManager>,
              MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, FailureCountManager>,
              MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, FailureCountManager>,
              MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, FailureCountManager>,
              MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, FailureCountManager>,
              MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, FailureCountManager>,
              MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, FailureCountManager>,
              MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, FailureCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, FailureCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, FailureCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, FailureCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, FailureCountManager>,
              HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, FailureCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, FailureCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, FailureCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, FailureCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, FailureCountManager>,
              RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, FailureCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, FailureCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, FailureCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, FailureCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, FailureCountManager>,
              RandomMinMax>(this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, FailureCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, FailureCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, FailureCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, FailureCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, FailureCountManager>,
              MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, FailureCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, FailureCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, FailureCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, FailureCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, FailureCountManager>,
              Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, FailureCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, FailureCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, FailureCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, FailureCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, FailureCountManager>,
              Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, FailureCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, FailureCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, FailureCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, FailureCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, FailureCountManager>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, FailureCountManager>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, FailureCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, FailureCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, FailureCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, FailureCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, FailureCountManager>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, FailureCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "gwdeg" || var_ordering == "GWDEG" ||
        var_ordering == "Global Weighted Degree" ||
        var_ordering == "GlobalWeightedDegree") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, ConflictCountManager>, AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, ConflictCountManager>, AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, ConflictCountManager>, AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, ConflictCountManager>, AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, ConflictCountManager>, AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, ConflictCountManager>, MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, ConflictCountManager>, MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, ConflictCountManager>, MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, ConflictCountManager>, MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, ConflictCountManager>, MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, ConflictCountManager>, MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, ConflictCountManager>, MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, ConflictCountManager>, MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, ConflictCountManager>, MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, ConflictCountManager>, MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, ConflictCountManager>, HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, ConflictCountManager>, HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, ConflictCountManager>, HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, ConflictCountManager>, HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, ConflictCountManager>, HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, ConflictCountManager>, RandomSplit>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, ConflictCountManager>, RandomSplit>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, ConflictCountManager>, RandomSplit>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, ConflictCountManager>, RandomSplit>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, ConflictCountManager>, RandomSplit>(
              this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, ConflictCountManager>, RandomMinMax>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, ConflictCountManager>, RandomMinMax>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, ConflictCountManager>, RandomMinMax>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, ConflictCountManager>, RandomMinMax>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, ConflictCountManager>, RandomMinMax>(
              this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, ConflictCountManager>, MinWeightValue>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, ConflictCountManager>, MinWeightValue>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, ConflictCountManager>, MinWeightValue>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, ConflictCountManager>, MinWeightValue>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, ConflictCountManager>, MinWeightValue>(
              this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, ConflictCountManager>, Guided<MinValue>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, ConflictCountManager>, Guided<MinValue>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, ConflictCountManager>, Guided<MinValue>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, ConflictCountManager>, Guided<MinValue>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, ConflictCountManager>, Guided<MinValue>>(
              this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, ConflictCountManager>, Guided<MaxValue>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, ConflictCountManager>, Guided<MaxValue>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, ConflictCountManager>, Guided<MaxValue>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, ConflictCountManager>, Guided<MaxValue>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, ConflictCountManager>, Guided<MaxValue>>(
              this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 1, ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 2, ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 3, ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 4, ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MaxWeight, 5, ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "dom/gwdeg" || var_ordering == "dgwdeg" ||
        var_ordering == "DGWDEG" ||
        var_ordering == "Domain Over Global Weighted Degree" ||
        var_ordering == "DomainOverGlobalWeightedDegree") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, ConflictCountManager>,
              AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, ConflictCountManager>,
              MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, ConflictCountManager>,
              MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, ConflictCountManager>,
              HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, ConflictCountManager>,
              RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, ConflictCountManager>,
              RandomMinMax>(this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, ConflictCountManager>,
              MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, ConflictCountManager>,
              Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, ConflictCountManager>,
              Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "Impact" || var_ordering == "impact") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 1, ImpactManager>,
                                     AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 2, ImpactManager>,
                                     AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 3, ImpactManager>,
                                     AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 4, ImpactManager>,
                                     AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 5, ImpactManager>,
                                     AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 1, ImpactManager>,
                                     MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 2, ImpactManager>,
                                     MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 3, ImpactManager>,
                                     MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 4, ImpactManager>,
                                     MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 5, ImpactManager>,
                                     MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 1, ImpactManager>,
                                     MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 2, ImpactManager>,
                                     MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 3, ImpactManager>,
                                     MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 4, ImpactManager>,
                                     MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 5, ImpactManager>,
                                     MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 1, ImpactManager>,
                                     HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 2, ImpactManager>,
                                     HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 3, ImpactManager>,
                                     HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 4, ImpactManager>,
                                     HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 5, ImpactManager>,
                                     HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 1, ImpactManager>,
                                     RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 2, ImpactManager>,
                                     RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 3, ImpactManager>,
                                     RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 4, ImpactManager>,
                                     RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 5, ImpactManager>,
                                     RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 1, ImpactManager>,
                                     RandomMinMax>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 2, ImpactManager>,
                                     RandomMinMax>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 3, ImpactManager>,
                                     RandomMinMax>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 4, ImpactManager>,
                                     RandomMinMax>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 5, ImpactManager>,
                                     RandomMinMax>(this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 1, ImpactManager>,
                                     MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 2, ImpactManager>,
                                     MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 3, ImpactManager>,
                                     MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 4, ImpactManager>,
                                     MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 5, ImpactManager>,
                                     MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 1, ImpactManager>,
                                     Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 2, ImpactManager>,
                                     Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 3, ImpactManager>,
                                     Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 4, ImpactManager>,
                                     Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 5, ImpactManager>,
                                     Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 1, ImpactManager>,
                                     Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 2, ImpactManager>,
                                     Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 3, ImpactManager>,
                                     Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 4, ImpactManager>,
                                     Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 5, ImpactManager>,
                                     Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 1, ImpactManager>,
                                     Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 2, ImpactManager>,
                                     Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 3, ImpactManager>,
                                     Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 4, ImpactManager>,
                                     Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 5, ImpactManager>,
                                     Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 1, ImpactManager>,
                                     Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 2, ImpactManager>,
                                     Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 3, ImpactManager>,
                                     Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 4, ImpactManager>,
                                     Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 5, ImpactManager>,
                                     Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 1, ImpactManager>,
                                     Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 2, ImpactManager>,
                                     Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 3, ImpactManager>,
                                     Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 4, ImpactManager>,
                                     Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinWeight, 5, ImpactManager>,
                                     Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinWeight, 1, ImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinWeight, 2, ImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinWeight, 3, ImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinWeight, 4, ImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinWeight, 5, ImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "IBS" || var_ordering == "ibs" ||
        var_ordering == "Impact" || var_ordering == "impact") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 1, ImpactManager>, AnyValue>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 2, ImpactManager>, AnyValue>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 3, ImpactManager>, AnyValue>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 4, ImpactManager>, AnyValue>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 5, ImpactManager>, AnyValue>(
              this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 1, ImpactManager>, MinValue>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 2, ImpactManager>, MinValue>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 3, ImpactManager>, MinValue>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 4, ImpactManager>, MinValue>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 5, ImpactManager>, MinValue>(
              this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 1, ImpactManager>, MaxValue>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 2, ImpactManager>, MaxValue>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 3, ImpactManager>, MaxValue>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 4, ImpactManager>, MaxValue>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 5, ImpactManager>, MaxValue>(
              this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 1, ImpactManager>, HalfSplit>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 2, ImpactManager>, HalfSplit>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 3, ImpactManager>, HalfSplit>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 4, ImpactManager>, HalfSplit>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 5, ImpactManager>, HalfSplit>(
              this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 1, ImpactManager>, RandomSplit>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 2, ImpactManager>, RandomSplit>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 3, ImpactManager>, RandomSplit>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 4, ImpactManager>, RandomSplit>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 5, ImpactManager>, RandomSplit>(
              this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 1, ImpactManager>, RandomMinMax>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 2, ImpactManager>, RandomMinMax>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 3, ImpactManager>, RandomMinMax>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 4, ImpactManager>, RandomMinMax>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 5, ImpactManager>, RandomMinMax>(
              this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 1, ImpactManager>,
              MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 2, ImpactManager>,
              MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 3, ImpactManager>,
              MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 4, ImpactManager>,
              MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 5, ImpactManager>,
              MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 1, ImpactManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 2, ImpactManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 3, ImpactManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 4, ImpactManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 5, ImpactManager>,
              Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 1, ImpactManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 2, ImpactManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 3, ImpactManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 4, ImpactManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 5, ImpactManager>,
              Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 1, ImpactManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 2, ImpactManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 3, ImpactManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 4, ImpactManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 5, ImpactManager>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 1, ImpactManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 2, ImpactManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 3, ImpactManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 4, ImpactManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 5, ImpactManager>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 1, ImpactManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 2, ImpactManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 3, ImpactManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 4, ImpactManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 5, ImpactManager>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 1, ImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 2, ImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 3, ImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 4, ImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainTimesWeight, 5, ImpactManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "RIBS" || var_ordering == "ribs" ||
        var_ordering == "RealImpact" || var_ordering == "realimpact") {

      // std::cout << branching << std::endl;

      if (branching == "xdiv") {
        if (randomness <= 1) {
          heu = new ImpactBasedSearch<1>(this, 0.00001);
        } else if (randomness <= 2) {
          heu = new ImpactBasedSearch<2>(this, 0.00001);
        } else if (randomness <= 3) {
          heu = new ImpactBasedSearch<3>(this, 0.00001);
        } else if (randomness <= 4) {
          heu = new ImpactBasedSearch<4>(this, 0.00001);
        } else if (randomness <= 5) {
          heu = new ImpactBasedSearch<5>(this, 0.00001);
        }
      } else if (branching == "xint") {
        if (randomness <= 1) {
          heu = new ImpactBasedSearch<1>(this, 0.9);
        } else if (randomness <= 2) {
          heu = new ImpactBasedSearch<2>(this, 0.9);
        } else if (randomness <= 3) {
          heu = new ImpactBasedSearch<3>(this, 0.9);
        } else if (randomness <= 4) {
          heu = new ImpactBasedSearch<4>(this, 0.9);
        } else if (randomness <= 5) {
          heu = new ImpactBasedSearch<5>(this, 0.9);
        }
      } else if (branching == "div") {
        if (randomness <= 1) {
          heu = new ImpactBasedSearch<1>(this, 0.001);
        } else if (randomness <= 2) {
          heu = new ImpactBasedSearch<2>(this, 0.001);
        } else if (randomness <= 3) {
          heu = new ImpactBasedSearch<3>(this, 0.001);
        } else if (randomness <= 4) {
          heu = new ImpactBasedSearch<4>(this, 0.001);
        } else if (randomness <= 5) {
          heu = new ImpactBasedSearch<5>(this, 0.001);
        }
      } else if (branching == "int") {
        if (randomness <= 1) {
          heu = new ImpactBasedSearch<1>(this, 0.3);
        } else if (randomness <= 2) {
          heu = new ImpactBasedSearch<2>(this, 0.3);
        } else if (randomness <= 3) {
          heu = new ImpactBasedSearch<3>(this, 0.3);
        } else if (randomness <= 4) {
          heu = new ImpactBasedSearch<4>(this, 0.3);
        } else if (randomness <= 5) {
          heu = new ImpactBasedSearch<5>(this, 0.3);
        }
      } else if (branching == "xdivlex") {
        if (randomness <= 1) {
          heu = new ImpactBasedSearch<1>(this, 0.00001, false);
        } else if (randomness <= 2) {
          heu = new ImpactBasedSearch<2>(this, 0.00001, false);
        } else if (randomness <= 3) {
          heu = new ImpactBasedSearch<3>(this, 0.00001, false);
        } else if (randomness <= 4) {
          heu = new ImpactBasedSearch<4>(this, 0.00001, false);
        } else if (randomness <= 5) {
          heu = new ImpactBasedSearch<5>(this, 0.00001, false);
        }
      } else if (branching == "xintlex") {
        if (randomness <= 1) {
          heu = new ImpactBasedSearch<1>(this, 0.9, false);
        } else if (randomness <= 2) {
          heu = new ImpactBasedSearch<2>(this, 0.9, false);
        } else if (randomness <= 3) {
          heu = new ImpactBasedSearch<3>(this, 0.9, false);
        } else if (randomness <= 4) {
          heu = new ImpactBasedSearch<4>(this, 0.9, false);
        } else if (randomness <= 5) {
          heu = new ImpactBasedSearch<5>(this, 0.9, false);
        }
      } else if (branching == "divlex") {
        if (randomness <= 1) {
          heu = new ImpactBasedSearch<1>(this, 0.001, false);
        } else if (randomness <= 2) {
          heu = new ImpactBasedSearch<2>(this, 0.001, false);
        } else if (randomness <= 3) {
          heu = new ImpactBasedSearch<3>(this, 0.001, false);
        } else if (randomness <= 4) {
          heu = new ImpactBasedSearch<4>(this, 0.001, false);
        } else if (randomness <= 5) {
          heu = new ImpactBasedSearch<5>(this, 0.001, false);
        }
      } else if (branching == "intlex") {
        if (randomness <= 1) {
          heu = new ImpactBasedSearch<1>(this, 0.3, false);
        } else if (randomness <= 2) {
          heu = new ImpactBasedSearch<2>(this, 0.3, false);
        } else if (randomness <= 3) {
          heu = new ImpactBasedSearch<3>(this, 0.3, false);
        } else if (randomness <= 4) {
          heu = new ImpactBasedSearch<4>(this, 0.3, false);
        } else if (randomness <= 5) {
          heu = new ImpactBasedSearch<5>(this, 0.3, false);
        }
      } else if (branching == "lex") {
        if (randomness <= 1) {
          heu = new ImpactBasedSearch<1>(this, 0.1, false);
        } else if (randomness <= 2) {
          heu = new ImpactBasedSearch<2>(this, 0.1, false);
        } else if (randomness <= 3) {
          heu = new ImpactBasedSearch<3>(this, 0.1, false);
        } else if (randomness <= 4) {
          heu = new ImpactBasedSearch<4>(this, 0.1, false);
        } else if (randomness <= 5) {
          heu = new ImpactBasedSearch<5>(this, 0.1, false);
        }
      } else {
        if (randomness <= 1) {
          heu = new ImpactBasedSearch<1>(this, 0.1);
        } else if (randomness <= 2) {
          heu = new ImpactBasedSearch<2>(this, 0.1);
        } else if (randomness <= 3) {
          heu = new ImpactBasedSearch<3>(this, 0.1);
        } else if (randomness <= 4) {
          heu = new ImpactBasedSearch<4>(this, 0.1);
        } else if (randomness <= 5) {
          heu = new ImpactBasedSearch<5>(this, 0.1);
        }
      }
    }
    if (var_ordering == "ABS" || var_ordering == "abs" ||
        var_ordering == "activity based search" ||
        var_ordering == "Activity Based Search" ||
        var_ordering == "dom/activity") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, PruningCountManager>,
              AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, PruningCountManager>,
              AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, PruningCountManager>,
              AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, PruningCountManager>,
              AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, PruningCountManager>,
              AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, PruningCountManager>,
              MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, PruningCountManager>,
              MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, PruningCountManager>,
              MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, PruningCountManager>,
              MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, PruningCountManager>,
              MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, PruningCountManager>,
              MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, PruningCountManager>,
              MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, PruningCountManager>,
              MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, PruningCountManager>,
              MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, PruningCountManager>,
              MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, PruningCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, PruningCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, PruningCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, PruningCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, PruningCountManager>,
              HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, PruningCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, PruningCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, PruningCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, PruningCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, PruningCountManager>,
              RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, PruningCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, PruningCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, PruningCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, PruningCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, PruningCountManager>,
              RandomMinMax>(this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, PruningCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, PruningCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, PruningCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, PruningCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, PruningCountManager>,
              MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, PruningCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, PruningCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, PruningCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, PruningCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, PruningCountManager>,
              Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, PruningCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, PruningCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, PruningCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, PruningCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, PruningCountManager>,
              Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, PruningCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, PruningCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, PruningCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, PruningCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, PruningCountManager>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, PruningCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, PruningCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, PruningCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, PruningCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, PruningCountManager>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, PruningCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, PruningCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, PruningCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, PruningCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, PruningCountManager>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 1, PruningCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 2, PruningCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 3, PruningCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 4, PruningCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomainOverWeight, 5, PruningCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "mindomain" || var_ordering == "MinDomain" ||
        var_ordering == "min domain" || var_ordering == "minimum domain") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomMinMax>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomMinMax>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomMinMax>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomMinMax>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomMinMax>(this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MinValue>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MinValue>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MinValue>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MinValue>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MinValue>>(
              this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MaxValue>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MaxValue>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MaxValue>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MaxValue>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MaxValue>>(
              this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, Guided<RandomMinMax>>(
                  this);
        } else if (randomness <= 2) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, Guided<RandomMinMax>>(
                  this);
        } else if (randomness <= 3) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, Guided<RandomMinMax>>(
                  this);
        } else if (randomness <= 4) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, Guided<RandomMinMax>>(
                  this);
        } else if (randomness <= 5) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, Guided<RandomMinMax>>(
                  this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomain>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomain>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomain>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomain>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomain>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "Neighbor" || var_ordering == "Neighbor" ||
        var_ordering == "Neighbour" || var_ordering == "neighbour") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              RandomMinMax>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              RandomMinMax>(this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 1,
                                 ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 2,
                                 ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 3,
                                 ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 4,
                                 ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericNeighborDVO<SelfPlusAverage, MinDomainOverWeight, 5,
                                 ConflictCountManager>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "mindom->maxdeg" ||
        var_ordering == "MinDomainMaxDegree") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, HalfSplit>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, HalfSplit>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, HalfSplit>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, HalfSplit>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, HalfSplit>(
              this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, RandomSplit>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, RandomSplit>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, RandomSplit>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, RandomSplit>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, RandomSplit>(
              this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, RandomMinMax>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, RandomMinMax>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, RandomMinMax>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, RandomMinMax>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, RandomMinMax>(
              this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, MinWeightValue>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, MinWeightValue>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, MinWeightValue>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, MinWeightValue>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>, MinWeightValue>(
              this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MinValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MinValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MinValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MinValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MinValue>>(this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MaxValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MaxValue>>(this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<RandomMinMax>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              Guided<RandomMinMax>>(this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<LexCombination<MinDomain, MaxDegree>>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    }
    if (var_ordering == "mindomain" || var_ordering == "MinDomain" ||
        var_ordering == "min domain" || var_ordering == "minimum domain" ||
        var_ordering == "MinDomnMaxDeg") {
      if (branching == "No" || branching == "no" || branching == "Any" ||
          branching == "any") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, AnyValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, AnyValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, AnyValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, AnyValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, AnyValue>(this);
        }
      }
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MinValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MinValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MinValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MinValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MinValue>(this);
        }
      }
      if (branching == "AntiLex" || branching == "antilex" ||
          branching == "maxvalue" || branching == "max value" ||
          branching == "MaxValue" || branching == "Max Value" ||
          branching == "antilexicographic" ||
          branching == "Antilexicographic") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MaxValue>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MaxValue>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MaxValue>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MaxValue>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, MaxValue>(this);
        }
      }
      if (branching == "HalfSplit" || branching == "halfsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, HalfSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, HalfSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, HalfSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, HalfSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, HalfSplit>(this);
        }
      }
      if (branching == "RandomSplit" || branching == "RandSplit" ||
          branching == "randomsplit" || branching == "randsplit") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomSplit>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomSplit>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomSplit>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomSplit>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomSplit>(this);
        }
      }
      if (branching == "RandomMinMax" || branching == "randomminmax" ||
          branching == "RandMinMax" || branching == "randminmax") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomMinMax>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomMinMax>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomMinMax>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomMinMax>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, RandomMinMax>(this);
        }
      }
      if (branching == "minweight" || branching == "MinWeight") {
        if (randomness <= 1) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, MinWeightValue>(this);
        } else if (randomness <= 2) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, MinWeightValue>(this);
        } else if (randomness <= 3) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, MinWeightValue>(this);
        } else if (randomness <= 4) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, MinWeightValue>(this);
        } else if (randomness <= 5) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, MinWeightValue>(this);
        }
      }
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MinValue>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MinValue>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MinValue>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MinValue>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MinValue>>(
              this);
        }
      }
      if (branching == "MaxVal+Guided" || branching == "maxval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MaxValue>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MaxValue>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MaxValue>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MaxValue>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>, Guided<MaxValue>>(
              this);
        }
      }
      if (branching == "MinWeight+Guided" || branching == "minweight+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MinWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MinWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MinWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MinWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MinWeightValue>>(this);
        }
      }
      if (branching == "MaxWeightVal+Guided" ||
          branching == "maxweightval+guided") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MaxWeightValue>>(this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MaxWeightValue>>(this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MaxWeightValue>>(this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MaxWeightValue>>(this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<GenericDVO<MinDomain>,
                                     Guided<MaxWeightValue>>(this);
        }
      }
      if (branching == "Random+Guided" || branching == "random+guided") {
        if (randomness <= 1) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, Guided<RandomMinMax>>(
                  this);
        } else if (randomness <= 2) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, Guided<RandomMinMax>>(
                  this);
        } else if (randomness <= 3) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, Guided<RandomMinMax>>(
                  this);
        } else if (randomness <= 4) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, Guided<RandomMinMax>>(
                  this);
        } else if (randomness <= 5) {
          heu =
              new GenericHeuristic<GenericDVO<MinDomain>, Guided<RandomMinMax>>(
                  this);
        }
      }
      if (branching == "Adpated" || branching == "adapted") {
        if (randomness <= 1) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomain>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 2) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomain>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 3) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomain>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 4) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomain>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        } else if (randomness <= 5) {
          heu = new GenericHeuristic<
              GenericDVO<MinDomain>,
              ConditionalOnSize<GuidedSplit<HalfSplit>, Guided<MinValue>>>(
              this);
        }
      }
    } else if (var_ordering == "COS") {
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        heu = new ConflictOrderedSearch<
            GenericDVO<MinDomainOverWeight, 2, FailureCountManager>, MinValue,
            MinValue>(this);
      }
      if (branching == "Adpated" || branching == "adapted") {
        heu = new ConflictOrderedSearch<
            GenericDVO<MinDomainOverWeight, 2, FailureCountManager>,
            Guided<MinValue>, Guided<MinValue>>(this);
      }
    } else if (var_ordering == "ECOS") {
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        heu = new ConflictOrderedSearch<
            GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>, MinValue,
            MinValue>(this);
      }
      if (branching == "Adpated" || branching == "adapted") {
        heu = new ConflictOrderedSearch<
            GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>,
            Guided<MinValue>, Guided<MinValue>>(this);
      }
    } else if (var_ordering == "LC") {
      if (branching == "Lex" || branching == "lex" || branching == "minvalue" ||
          branching == "min value" || branching == "MinValue" ||
          branching == "Min Value" || branching == "lexicographic" ||
          branching == "Lexicographic") {
        heu = new LastConflict<
            GenericDVO<MinDomainOverWeight, 2, FailureCountManager>, MinValue,
            MinValue, 1>(this);
      }
      if (branching == "Adpated" || branching == "adapted") {
        heu = new LastConflict<
            GenericDVO<MinDomainOverWeight, 2, FailureCountManager>,
            Guided<MinValue>, Guided<MinValue>, 1>(this);
      }
      // heu = new LastConflict < Lexicographic, MinValue, MinValue, 3 > (this);
    } else if (var_ordering == "ELC") {
      if (branching == "MinVal+Guided" || branching == "minval+guided") {
        heu = new LastConflict<
            GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
            Guided<MinValue>, Guided<MinValue>, 1>(this);
      } else if (branching == "Lex" || branching == "lex" ||
                 branching == "minvalue" || branching == "min value" ||
                 branching == "MinValue" || branching == "Min Value" ||
                 branching == "lexicographic" || branching == "Lexicographic") {
        heu = new LastConflict<
            GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>, MinValue,
            MinValue, 1>(this);
      } else if (branching == "Adpated" || branching == "adapted") {
        heu = new LastConflict<
            GenericDVO<MinDomainOverWeight, 2, ConflictCountManager>,
            Guided<MinValue>, Guided<MinValue>, 1>(this);
      }
      // heu = new LastConflict < Lexicographic, MinValue, MinValue, 3 > (this);
    }
  }

  if (!heu) {
    std::cout << " " << parameters.prefix_comment
              << " Warning, there is no known heuristic \"" << var_ordering
              << "/" << branching << "\"" << std::endl;
  }

  return heu;
}

Mistral::RestartPolicy *Mistral::Solver::restart_factory(std::string rpolicy) {
  RestartPolicy *pol;
  if (rpolicy == "luby")
    pol = new Luby();
  else if (rpolicy == "geom")
    pol = new Geometric();
  else
    pol = new NoRestart();
  return pol;
}

void Mistral::Solver::initialise_random_seed(const int seed) { usrand(seed); }

void Mistral::Solver::set_learning_on() {

  parameters.backjump = true;
  if (!base) {
    base = new ConstraintClauseBase(variables);
    add(Constraint(base));
    // add(base);
  }
}

void Mistral::Solver::set_time_limit(const double limit) {
  if (limit > 0) {
    parameters.limit = 1;
    parameters.time_limit = limit;
  }
}

void Mistral::Solver::set_goal(Goal *g) {
  if (g) {
    objective = g;

    if (g->is_optimization()) {
      consolidate_manager->id_obj = g->objective.id();
    }
  }
}

#ifdef _REF_SOL_
void Mistral::Solver::check_reference() {
  for (int i = 0; i < variables.size; ++i) {
    if (reference_solution[i] != NOVAL and
        !variables[i].contain(reference_solution[i])) {
      std::cout << "value " << reference_solution[i]
                << " consistent but removed from " << variables[i] << std::endl;
      exit(1);
    }
  }
  std::cout << "REF OK\n";
}
#endif

void Mistral::Solver::check_constraint_graph_integrity() {
  // for each constraint, check if the set of active variables corresponds to
  // unbound vars
  for (unsigned int i = 0; i < constraints.size; ++i) {
    // if(constraints[i].is_active())
    constraints[i].check_active();
  }

  // for each variable, check that
  // 1/ no constraint is listed twice
  // 2/ for each constraint in the list for this var:
  //   a/ the index of the constraint for that var corresponds to the rank of
  //   var in its scope b/ the constraint has at least 2 unbound vars, or does
  //   not enforce nfc1
  BitSet cons_list(0, constraints.size, BitSet::empt);
  int trig, cons;
  Constraint c;
  Variable *scope;
  for (unsigned int i = 0; i < variables.size; ++i) {
    cons_list.clear();

    for (trig = _VALUE_; trig <= _DOMAIN_; ++trig) {

      for (cons = constraint_graph[i].on[trig].size; --cons >= 0;) {

        c = constraint_graph[i].on[trig][cons];

        if (cons_list.contain(c.id())) {

          std::cout << "Warning: " << c
                    << " is listed at least twice: " << std::endl
                    << constraint_graph[i].on[0] << std::endl
                    << constraint_graph[i].on[1] << std::endl
                    << constraint_graph[i].on[2] << std::endl;

          exit(1);

        } else {
          cons_list.add(c.id());
        }

        scope = c.get_scope();

        if (scope[c.index()].id() != (int)i) {

          std::cout << "Warning: incorrect variable indexing: "
                    << scope[c.index()] << "'s " << c << " is posted on "
                    << variables[i] << std::endl;

          exit(1);
        }

        if (c.rank() != cons) {

          std::cout << "Warning: incorrect list indexing: " << c << " is "
                    << cons << "th in the list of " << variables[i]
                    << ", but indexed " << c.rank() << "th." << std::endl;

          exit(1);
        }
      }
    }
  }
}

Mistral::SolverCmdLine::SolverCmdLine(const std::string &message,
                                      const char delimiter,
                                      const std::string &version,
                                      bool helpAndVersion)
    : CmdLine(message, delimiter, version, helpAndVersion) {
  initialise();
}

Mistral::SolverCmdLine::~SolverCmdLine() {

  delete fileArg;
  delete seedArg;
  delete timeArg;
  // delete printArg;
  //  delete printsolArg;
  delete allsolArg;
  delete printsolArg;
  delete printstaArg;
  delete printmodArg;
  delete printparArg;
  delete printinsArg;
  delete verbosityArg;
  delete randomizationArg;
  delete rewriteArg;
  delete restartArg;
  delete factorArg;
  delete baseArg;
  delete decayArg;
  delete forgetArg;
  delete incrementArg;
  delete learningArg;
  delete branchingArg;
  delete orderingArg;
  delete pcommentArg;
  delete pstatArg;
  delete pobjectiveArg;
  delete psolutionArg;
  delete poutcomeArg;

  delete r_allowed;
  delete vo_allowed;
  delete bo_allowed;
}

void Mistral::SolverCmdLine::initialise() {

  // INPUT FILE
  fileArg = new TCLAP::UnlabeledValueArg<std::string>("file", "instance file",
                                                      true, "", "string");
  add(*fileArg);

  // COMMENT INDICATOR
  pcommentArg = new TCLAP::ValueArg<std::string>(
      "", "prefix_comment", "output comments prefix", false, "c", "string");
  add(*pcommentArg);
  // STATS INDICATOR
  pstatArg = new TCLAP::ValueArg<std::string>(
      "", "prefix_stat", "output stats prefix", false, "d", "string");
  add(*pstatArg);
  // NEW OBJ INDICATOR
  pobjectiveArg = new TCLAP::ValueArg<std::string>(
      "", "prefix_obj", "output new objective prefix", false, "o", "string");
  add(*pobjectiveArg);
  // SOLUTION INDICATOR
  psolutionArg = new TCLAP::ValueArg<std::string>(
      "", "prefix_sol", "output solution prefix", false, "v", "string");
  add(*psolutionArg);
  // OUTCOME INDICATOR
  poutcomeArg = new TCLAP::ValueArg<std::string>(
      "", "prefix_res", "output search outcome prefix", false, "s", "string");
  add(*poutcomeArg);

  // TIME LIMIT
  timeArg = new TCLAP::ValueArg<double>("t", "time_limit", "time limit", false,
                                        0, "double");
  add(*timeArg);

  // // PRINTING OPTIONS
  // std::vector<std::string> pallowed;
  // pallowed.push_back("params");
  // pallowed.push_back("instance");
  // pallowed.push_back("model");
  // pallowed.push_back("sol");
  // pallowed.push_back("stat");
  // TCLAP::ValuesConstraint<std::string> p_allowed( pallowed );
  // printArg = new TCLAP::MultiArg<std::string>("p","print","objects to
  // print",false,&p_allowed); add( *printArg );

  // VERBOSITY LEVEL
  verbosityArg = new TCLAP::ValueArg<int>("v", "verbosity", "verbosity level",
                                          false, 1, "int");
  add(*verbosityArg);

  // RANDOM SEED
  seedArg =
      new TCLAP::ValueArg<int>("s", "seed", "random seed", false, 12345, "int");
  add(*seedArg);

  // HEURISTIC RANDOMIZATION
  randomizationArg = new TCLAP::ValueArg<int>(
      "z", "randomization", "randomization level", false, 0, "int");
  add(*randomizationArg);

  // WHETHER WE USE REWRITING OPTIONS
  rewriteArg = new TCLAP::SwitchArg(
      "w", "rewrite", "use rewriting during preprocessing", false);
  add(*rewriteArg);

  // WHETHER WE WANT TO FIND ALL SOLUTIONS
  allsolArg = new TCLAP::SwitchArg("a", "all", "find all solutions", false);
  add(*allsolArg);

  // RESTART POLICY
  std::vector<std::string> rallowed;
  rallowed.push_back("no");
  rallowed.push_back("geom");
  rallowed.push_back("luby");
  r_allowed = new TCLAP::ValuesConstraint<std::string>(rallowed);
  restartArg = new TCLAP::ValueArg<std::string>(
      "r", "restart", "restart policy", false, "geom", r_allowed);
  add(*restartArg);

  // RESTART FACTOR
  factorArg = new TCLAP::ValueArg<double>("m", "factor", "restart factor",
                                          false, 1.05, "double");
  add(*factorArg);

  // RESTART BASE
  baseArg =
      new TCLAP::ValueArg<int>("e", "base", "restart base", false, 200, "int");
  add(*baseArg);

  // DECAY FACTOR
  decayArg = new TCLAP::ValueArg<double>("d", "decay", "decay factor", false,
                                         .96, "double");
  add(*decayArg);

  // FORGETFULNESS
  forgetArg = new TCLAP::ValueArg<double>("g", "forget", "clause forgetfulness",
                                          false, .75, "double");
  add(*forgetArg);

  // ACTIVITY INCREMENT
  incrementArg = new TCLAP::ValueArg<double>(
      "i", "increment", "activity increment", false, .012, "double");
  add(*incrementArg);

  // USE CLAUSE LEARNING
  learningArg = new TCLAP::SwitchArg("l", "learning",
                                     "Switch on clause learning (CDCL)", false);
  add(*learningArg);

  // PRINT MODEL
  printmodArg = new TCLAP::SwitchArg("", "print_mod", "Print the model", false);
  add(*printmodArg);

  // PRINT STATISTICS
  printstaArg =
      new TCLAP::SwitchArg("", "print_sta", "Print the statistics", false);
  add(*printstaArg);

  // PRINT INSTANCE
  printinsArg =
      new TCLAP::SwitchArg("", "print_ins", "Print the instance", false);
  add(*printinsArg);

  // PRINT PARAMETERS
  printparArg =
      new TCLAP::SwitchArg("", "print_par", "Print the parameters", false);
  add(*printparArg);

  // PRINT SOLUTION
  printsolArg = new TCLAP::SwitchArg("", "print_sol",
                                     "Print the solution, if found", false);
  add(*printsolArg);

  // VARIABLE ORDERING
  std::vector<std::string> voallowed;
  voallowed.push_back("dom/deg");
  voallowed.push_back("dom/wdeg");
  voallowed.push_back("dom/gwdeg");
  voallowed.push_back("dom/pruning");
  voallowed.push_back("dom/activity");
  voallowed.push_back("ibs");
  voallowed.push_back("IBS");
  voallowed.push_back("ribs");
  voallowed.push_back("RIBS");
  voallowed.push_back("impact");
  voallowed.push_back("Impact");
  voallowed.push_back("realimpact");
  voallowed.push_back("RealImpact");
  voallowed.push_back("abs");
  voallowed.push_back("ABS");
  voallowed.push_back("activity");
  voallowed.push_back("neighbor");
  voallowed.push_back("mindomain");
  voallowed.push_back("maxdegree");
  voallowed.push_back("lexicographic");
  voallowed.push_back("input_order");
  voallowed.push_back("first_fail");
  voallowed.push_back("anti_first_fail");
  voallowed.push_back("smallest");
  voallowed.push_back("largest");
  voallowed.push_back("occurrence");
  voallowed.push_back("most_constrained");
  voallowed.push_back("max_regret");
  voallowed.push_back("impact");

  voallowed.push_back("MBAdom/deg");
  voallowed.push_back("MBAdom/wdeg");
  voallowed.push_back("MBAdom/gwdeg");
  voallowed.push_back("MBAdom/pruning");
  voallowed.push_back("MBAdom/activity");
  voallowed.push_back("MBAibs");
  voallowed.push_back("MBAactivity");
  voallowed.push_back("MBAmindomain");
  voallowed.push_back("MBAmaxdegree");
  voallowed.push_back("MBAimpact");

  voallowed.push_back("COS");
  voallowed.push_back("ECOS");
  voallowed.push_back("LC");
  voallowed.push_back("ELC");

  vo_allowed = new TCLAP::ValuesConstraint<std::string>(voallowed);
  orderingArg = new TCLAP::ValueArg<std::string>(
      "c", "choice", "variable selection", false, "dom/gwdeg", vo_allowed);
  add(*orderingArg);

  // VALUE ORDERING
  std::vector<std::string> boallowed;
  boallowed.push_back("minvalue");
  boallowed.push_back("maxvalue");
  boallowed.push_back("minweight");
  boallowed.push_back("maxweight");
  boallowed.push_back("halfsplit");
  boallowed.push_back("reversesplit");
  boallowed.push_back("randomsplit");
  boallowed.push_back("randomspivot");
  boallowed.push_back("guidedsplit");
  boallowed.push_back("random");
  boallowed.push_back("randminmax");
  boallowed.push_back("guided");
  boallowed.push_back("minweight+guided");
  boallowed.push_back("maxweight+guided");
  boallowed.push_back("minval+guided");
  boallowed.push_back("maxval+guided");
  boallowed.push_back("random+guided");
  boallowed.push_back("adapted");

  boallowed.push_back("indomain_min");
  boallowed.push_back("indomain_max");
  boallowed.push_back("indomain_middle");
  boallowed.push_back("indomain_median");
  boallowed.push_back("indomain_random");
  boallowed.push_back("indomain_split");
  boallowed.push_back("indomain_reverse_split");
  boallowed.push_back("indomain_interval");

  boallowed.push_back("int");
  boallowed.push_back("xint");
  boallowed.push_back("div");
  boallowed.push_back("xdiv");
  boallowed.push_back("intlex");
  boallowed.push_back("xintlex");
  boallowed.push_back("divlex");
  boallowed.push_back("xdivlex");
  boallowed.push_back("lex");

  bo_allowed = new TCLAP::ValuesConstraint<std::string>(boallowed);
  branchingArg = new TCLAP::ValueArg<std::string>(
      "b", "branching", "value ordering", false, "minval+guided", bo_allowed);
  add(*branchingArg);
}


void Mistral::SolverCmdLine::set_parameters(Mistral::Solver& s) {

  s.parameters.verbosity = verbosityArg->getValue();
  s.parameters.restart_factor = factorArg->getValue();
  s.parameters.restart_base = baseArg->getValue();
  s.parameters.restart_limit = s.parameters.restart_base;
  s.parameters.time_limit = timeArg->getValue();
  s.parameters.activity_decay = decayArg->getValue();
  s.parameters.backjump = learningArg->getValue();
  s.parameters.activity_increment = incrementArg->getValue();
  s.parameters.forgetfulness = forgetArg->getValue();
  s.parameters.prefix_comment = pcommentArg->getValue();
  s.parameters.prefix_statistics = pstatArg->getValue();
  s.parameters.prefix_objective = pobjectiveArg->getValue();
  s.parameters.prefix_solution = psolutionArg->getValue();
  s.parameters.prefix_outcome = poutcomeArg->getValue();
  // s.parameters.find_all = allsolArg->getValue();
}

std::string Mistral::SolverCmdLine::get_value_ordering() {
  return branchingArg->getValue();
}

std::string Mistral::SolverCmdLine::get_variable_ordering() {
  return orderingArg->getValue();
}

std::string Mistral::SolverCmdLine::get_restart_policy() {
  return restartArg->getValue();
}

std::string Mistral::SolverCmdLine::get_filename() {
  return fileArg->getValue(); //.c_str();
}

int Mistral::SolverCmdLine::get_seed() { return seedArg->getValue(); }

int Mistral::SolverCmdLine::get_randomization() {
  return randomizationArg->getValue();
}

bool Mistral::SolverCmdLine::print_model() {
  // init_print();
  return printmodArg->getValue();
  // return true;
}

bool Mistral::SolverCmdLine::print_solution() {
  // init_print();
  return printsolArg->getValue();
  // return true;
}

bool Mistral::SolverCmdLine::print_parameters() {
  // init_print();
  return printparArg->getValue();
  // return true;
}

bool Mistral::SolverCmdLine::print_instance() {
  // init_print();
  return printinsArg->getValue();
  // return true;
}

bool Mistral::SolverCmdLine::print_statistics() {
  // init_print();
  return printstaArg->getValue();
  // return true;
}

bool Mistral::SolverCmdLine::use_rewrite() { return rewriteArg->getValue(); }

bool Mistral::SolverCmdLine::enumerate_solutions() {
  return allsolArg->getValue();
}

#ifdef _CHECK_NOGOOD

void Mistral::Solver::store_reason(Explanation *expl, Atom a) {
  Explanation::iterator stop;
  Explanation::iterator lit = expl->get_reason_for(
      a, ((a != NULL_ATOM) ? assignment_level[a] : level), stop);

  Vector<Literal> ng;
  while (lit < stop) {
    ng.add(*lit);
    ++lit;
  }

  if (a != NULL_ATOM) {

    Literal l = literal(variables[a]);

    // if(l == 83) {
    //   std::cout << variables[a].get_domain() << " " << a << std::endl;
    //   exit(1);
    // }

    ng.add(l);
  }

  if (expl->is_clause())
    nogood_origin.add(NULL);
  else
    nogood_origin.add(expl);
  nogood_clause.add(ng);
  atom.add(a);
  node_num.add(statistics.num_filterings);
}

void Mistral::Solver::store_nogood(Vector< Literal >& lc) {
  atom.add(NULL_ATOM);
  node_num.add(statistics.num_filterings);
  nogood_clause.add(lc);
  nogood_origin.add(NULL);
}


void Mistral::Solver::read_solution(const char* fname) {
  std::ifstream fsol(fname);

  int val;
  int nvars = 0;
  fsol >> nvars;

  while (nvars--) {
    fsol >> val;
    solution.add(val);
  }
}

void Mistral::Solver::check_nogoods() {
  for (int i = 0; i < nogood_clause.size; ++i) {
    std::cout << "check " << node_num[i] << " " << atom[i] << " "
              << nogood_clause[i].size << " " << nogood_clause[i]
              << " produced by ";
    std::cout << " " << ((int *)(nogood_origin[i])) << " ";
    std::cout.flush();

    if (nogood_origin[i])
      std::cout << nogood_origin[i];
    else
      std::cout << "learning";

    bool ok = false;
    for (int j = 0; j < nogood_clause[i].size; ++j) {
      ok |= (solution[UNSIGNED(nogood_clause[i][j])] ==
             SIGN(nogood_clause[i][j]));
    }

    if (!ok) {
      std::cout << " WRONG NOGOOD!!\n";
      exit(1);
    }
    std::cout << " ok\n";
  }
}
#endif
