 
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


#include <limits>

#include <mistral_search.hpp>

#include <mistral_sat.hpp>
#include <mistral_variable.hpp>


//#define _COUNT_BRANCHES



Mistral::ConsolidateListener::ConsolidateListener(Mistral::Solver *s) 
  : VariableListener(), ConstraintListener(), solver(s)
{
  sequence = &(solver->sequence);
  constraints.initialise(solver->variables.size);


  // int n = solver->variables.size;
  // for(int i=0; i<n; ++i) {
  //   Vector< Constraint > neighborhood;
  //   neighborhood.initialise(solver->variables[i].get_degree());
  //   for(int k=0; k<3; ++k) {
  //     for(int j=solver->constraint_graph[i].on[k].size-1; j>=0; --j)
  //    	neighborhood.add(solver->constraint_graph[i].on[k][j]);
  //   }
  //   constraints.add(neighborhood);
  // }

  //std::cout << "Variables size: " << solver->variables.size << std::endl;

  Constraint c;
  Variable *scope;
  int i, j, arity, n=solver->variables.size;

  for(i=0; i<n; ++i) {
    constraints[i].initialise(solver->variables[i].get_degree());
  }

  n = solver->constraints.size;
  for(i=0; i<n; ++i) {
    c = solver->constraints[i];
    scope = c.get_scope();
    arity = c.arity();
    for(j=0; j<arity; ++j) {
      if(!scope[j].is_ground()) {
	//std::cout << " " << scope[j] << std::endl;
	c.set_index(j);
	constraints[scope[j].id()].add(c);
      }
    }
  }

  id_obj = -1;
  if(solver->objective) {
    id_obj = solver->objective->objective.id();
  }

  //std::cout << constraints[66] << std::endl;

  solver->add((VariableListener*)this);
  solver->add((ConstraintListener*)this);
}


Mistral::ConsolidateListener::~ConsolidateListener() {
  //std::cout << "in delete consolidate manager" << std::endl;
}

void Mistral::ConsolidateListener::notify_add_var() {
  Variable x = solver->variables.back();
  while((int)(constraints.size)<x.id()) {
    Vector< Constraint > neighborhood;
    constraints.add(neighborhood);
  }
  Vector< Constraint > neighborhood;
  neighborhood.initialise(x.get_degree());
  for(int k=0; k<3; ++k) {
    for(int j=solver->constraint_graph[x.id()].on[k].size-1; j>=0; --j)
      neighborhood.add(solver->constraint_graph[x.id()].on[k][j]);
  }
  constraints.add(neighborhood);
}

void Mistral::ConsolidateListener::notify_post (Constraint c) {};

void Mistral::ConsolidateListener::notify_relax(Constraint c) {};

void Mistral::ConsolidateListener::notify_add_con(Constraint c) {

  //std::cout << "ADD CON " << c << std::endl;

  Variable *scope = c.get_scope();
  int arity = c.arity();
  for(int i=0; i<arity; ++i) {
    constraints[scope[i].id()].add(c);
    constraints[scope[i].id()].back().set_index(i);
  }
}

void Mistral::ConsolidateListener::notify_change(const int idx) {

  //std::cout << "BEG REACT TO CHANGE ON " << solver->variables[idx] << std::endl;

  Variable X = solver->variables[idx];
  int ids = sequence->index(idx);
  if(ids>=0) sequence->list_[ids] = X;

  //if(idx==66) std::cout << constraints[idx] << std::endl;

  for(int i=constraints[idx].size; --i>=0;) {
    //std::cout << constraints[idx][i] << std::endl;
    constraints[idx][i].consolidate_var();
    //std::cout << constraints[idx][i] << std::endl;
  }

  //std::cout << "IDX=" << idx << " IDO=" << id_obj << std::endl;

  if(idx==id_obj) {
    solver->objective->objective = X;
  }

  //std::cout << "END REACT TO CHANGE ON " << solver->variables[idx] << std::endl;
}


Mistral::HeuristicPoolManager::HeuristicPoolManager(Solver *s) : solver(s) {// }

  //std::cout << " c add restart listener" << std::endl;
  
  heu_index = 1;
  solver->add((RestartListener*)this);
}

Mistral::HeuristicPoolManager::~HeuristicPoolManager() {// }
  for(unsigned int i=0; i<pool.size; ++i) {
    if(solver->heuristic != pool[i]) {
      delete [] pool[i];
    }
  }
  solver->remove((RestartListener*)this);
}


void Mistral::HeuristicPoolManager::notify_restart(const double prog) {
  //std::cout << " c notify restart (3): " << solver->statistics.num_restarts << std::endl;
  if(heu_index<(int)(pool.size)) {
    if(prog > 0.0) {
      counter = threshold;
    } else if(--counter <= 0) {
      counter = threshold;
      std::cout << " c switch heuristic!\n";
      solver->heuristic = pool[heu_index++];
    }
  }
  //std::cout << " c " << counter << std::endl;
}


// //Mistral::LiteralActivityManager::LiteralActivityManager(Solver *s, void *a) 
// Mistral::LiteralActivityManager::LiteralActivityManager(Solver *s) 
//   : solver(s) {

//   lit_activity = solver->lit_activity.stack_;
//   var_activity = solver->var_activity.stack_;
//   n_vars = solver->variables.size;


//   //   n_vars = solver->base->scope.size;

//   // if(solver->base) {
//   //   lit_activity = solver->base->lit_activity.stack_;
//   //   var_activity = solver->base->var_activity.stack_;
//   //   n_vars = solver->base->scope.size;
//   // } else {
//   //   n_vars = solver->variables.size;
//   //   lit_activity = new double[2*n_vars];
//   //   var_activity = new double[n_vars];
//   //   std::fill(lit_activity, lit_activity+2*n_vars, 0.012);
//   //   std::fill(var_activity, var_activity+n_vars, 0.024);
//   // }


//   //double activity_increment = parameters.activity_increment / (1 << clause.size);
  
//   // if(activity_increment > 0.0) {
//   //   int i=clause.size;
//    //   while(i--) {

//    //     //std::cout << clause << " " << activity_increment << std::endl;

//    //     lit_activity[clause[i]] += activity_increment;
//    //     var_activity[UNSIGNED(clause[i])] += activity_increment;
//    //   }
//    // }

   
//   decay = solver->parameters.activity_decay;
//   solver->add((DecisionListener*)this);
// }

// Mistral::LiteralActivityManager::~LiteralActivityManager() {
//   solver->remove((DecisionListener*)this);
//   // if(!solver->base) {
//   //   delete [] lit_activity;
//   //   delete [] var_activity;
//   // }
// }

// double *Mistral::LiteralActivityManager::get_weight() { return var_activity; }     

// void Mistral::LiteralActivityManager::notify_decision() {
//   int i=n_vars;
//   while(i--) {
//     //std::cout << i << " " << var_activity[i] << " -> ";
//     var_activity[i] *= decay;
//     //std::cout << var_activity[i] << std::endl;
//   }    
//   i=2*n_vars;
//   while(i--) lit_activity[i] *= decay;
// }    

Mistral::RestartPolicy::RestartPolicy(const unsigned int b) {
  base = b;
}

Mistral::RestartPolicy::~RestartPolicy() {
}

Mistral::NoRestart::NoRestart() 
  : RestartPolicy(-1)
{
}

Mistral::NoRestart::~NoRestart() {}

Mistral::Geometric::Geometric(const unsigned int b, const double f) 
  : RestartPolicy(b)
{
  increment = b;
  factor = f;
}

Mistral::Geometric::~Geometric() {}

Mistral::Luby::Luby(const unsigned int b) 
  : RestartPolicy(b)
{
  iteration = 0;
}

Mistral::Luby::~Luby() {}

Mistral::NoOrder::NoOrder(Solver *s) 
  : solver(s) {}

Mistral::NoOrder::~NoOrder() {}

Mistral::Variable Mistral::NoOrder::select() {
  return solver->sequence.back();
}

// Mistral::Lexicographic::Lexicographic(Solver *s) 
//   : solver(s) {
//   index.initialise(s->variables.size);
//   solver->add(this);
//   last.initialise(0,s);
// }

Mistral::Lexicographic::Lexicographic(Solver *s) 
  : solver(s) {
  index.initialise(s->variables.size);

  int n = solver->variables.size;
  std::fill(index.stack_, index.stack_+n, -1);

  //std::cout << n << std::endl;
 
  solver->add(this);
  last.initialise(s,0);
}

void Mistral::Lexicographic::initialise(VarStack< Variable, ReversibleNum<int> >& seq) {
  // int n = solver->variables.size;
  // std::fill(index.stack_, index.stack_+n, -1);
  //for(int i=0; i<seq.size; ++i) {
  if(order.empty())
    for(int i=seq.size; --i>=0;) {
      index[seq[i].id()] = order.size;
      order.add(seq[i]);
    }
}

//void Mistral::Lexicographic::initialise(Solver *s, void *a) {
void Mistral::Lexicographic::initialise(Solver *s) {
  solver = s;
  index.initialise(s->variables.size);

  int n = solver->variables.size;
  std::fill(index.stack_, index.stack_+n, -1);

  //std::cout << n << std::endl;

  solver->add(this);
  last.initialise(s,0);
}

Mistral::Lexicographic::~Lexicographic() {}

Mistral::Variable Mistral::Lexicographic::select() {

  // for(int i=0; i<order.size; ++i) {
  //   std::cout << order[i] << " in " << order[i].get_domain() << " ";
  // }
  // std::cout << std::endl << (int)last << std::endl;
  // for(int i=last; i<order.size; ++i) {
  //   std::cout << order[i] << " in " << order[i].get_domain() << " ";
  // }
  // std::cout << std::endl;

  while(last<(int)(order.size) && order[last].is_ground()) { 
    ++last;
  }
  return order[last];
}

void Mistral::Lexicographic::notify_change(const int idx) {
  int ido = index[idx];
  if(ido>=0) {
    order[ido] = solver->variables[idx];
  }
}

// std::ostream& operator<<(std::ostream& os, Mistral::BranchingHeuristic& x) {
//   return x.display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::BranchingHeuristic* x) {
//   return x->display(os);
// }


std::ostream& operator<<(std::ostream& os, Mistral::DecisionListener& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::RestartListener& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::SuccessListener& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::BacktrackListener& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::ConstraintListener& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::VariableListener& x) {
  return x.display(os);
}

// std::ostream& operator<<(std::ostream& os, Mistral::DecisionListener* x) {
//   return x->display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::RestartListener* x) {
//   return x->display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::SuccessListener* x) {
//   return x->display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::BacktrackListener* x) {
//   return x->display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::ConstraintListener* x) {
//   return x->display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::VariableListener* x) {
//   return x->display(os);
// }

double *weight_sorting_array;

int Mistral::decreasing_weight(const void *x, const void *y) {
    int _x = *(int*)x;
    int _y = *(int*)y;

    int r_value = 0;
    if(weight_sorting_array[_x] > weight_sorting_array[_y]) {
      r_value = -1;
    } else if(weight_sorting_array[_x] < weight_sorting_array[_y]) {
      r_value = 1;
    }

    return r_value;
  }

int mylog10(const double x){
  long int y = (long int)x;
  int log = 0;
  while(y>0) {
    ++log;
    y/=10;
  }
  return log;
}


int Mistral::ImpactManager::get_minweight_value(const Variable x) {
	
	int best_val = 0;
	int idx = x.id();
	int offset = init_min[idx];
	double *wgt = value_weight[idx];
	
	if(factor[idx]==1) {
		best_val = x.get_min();
		double min_weight = wgt[best_val-offset], aux_weight;
		int vali, vnxt=x.next(best_val);
		do {
			vali = vnxt;
			vnxt = x.next(vali);
			aux_weight = wgt[vali-offset]; 
			if(aux_weight < min_weight) {
				min_weight = aux_weight;
				best_val = vali;
			}
		} while(vali<vnxt);
	} else {
		int fact = factor[idx];
		int the_min = x.get_min();
		int the_max = x.get_max();	
		int best_group = (the_min-offset)/fact;
		int group = best_group;
		double min_weight = wgt[best_group];
		int lb = group*fact+offset;
		int ub = lb+fact-1;
		while(ub < the_max) {
			++group;
			lb += fact;
			ub += fact;
			
			if(x.intersect(lb, ub) && min_weight>wgt[group]) {
				min_weight = wgt[group];
				best_group = group;
			}
		}
		best_val = x.next(best_group*fact+offset-1);
	}
	
	return best_val;
}
int Mistral::ImpactManager::get_maxweight_value(const Variable x) {
	int best_val = 0;
	int idx = x.id();
	int offset = init_min[idx];
	double *wgt = value_weight[idx];

	if(factor[idx]==1) {
		best_val = x.get_min();
		double max_weight = wgt[best_val-offset], aux_weight;
		int vali, vnxt=x.next(best_val);
		do {
			vali = vnxt;
			vnxt = x.next(vali);
			aux_weight = wgt[vali-offset]; 
			if(aux_weight > max_weight) {
				max_weight = aux_weight;
				best_val = vali;
			}
		} while(vali<vnxt);
	} else {
		int fact = factor[idx];
		int the_min = x.get_min();
		int the_max = x.get_max();	
		int best_group = (the_min-offset)/fact;
		int group = best_group;
		double max_weight = wgt[best_group];
		int lb = group*fact+offset;
		int ub = lb+fact-1;
		while(ub < the_max) {
			++group;
			lb += fact;
			ub += fact;
		
			if(x.intersect(lb, ub) && max_weight<wgt[group]) {
				max_weight = wgt[group];
				best_group = group;
			}
		}
		best_val = x.next(best_group*fact+offset-1);
	}

	return best_val;
}


void Mistral::ImpactManager::notify_success() {
			
	// propagation went without wipe-out
	// - check if it was after a left or a right branch
	// - find out what was the decision/refutation
	int dec;
	int id;
	int i, n;
	Variable x, y;
	double residual_space;
	int size;
	if(solver->decisions.size>0) {
				
				
#ifdef _DEBUG_IMPACT
		std::cout << " left " << solver->trail_.size << std::endl;
#endif
				
				
				
		i = solver->trail_.back(5), n=solver->saved_vars.size;

		// left branch
		Decision branch = solver->decisions.back();
		x = branch.var;
		dec = x.id();
		int dec_type = branch.type();
		int dec_val = branch.value();
		
		// std::cout << solver->variables[dec] << " " << x << std::endl;
		// exit(1);

#ifdef _DEBUG_IMPACT
		std::cout << " notified of success after the " << num_probes[dec] << "th left branch on " << solver->variables[dec] << "\n";
#endif

		residual_space = 1.0;
		while(i<n) {	
			id = solver->saved_vars[i];
			y = solver->variables[id];
	    
#ifdef _DEBUG_IMPACT2
			std::cout << " -> " << y << y.get_domain() << " lost " << y.get_reduction() << " values (";
#endif
	    
			size = y.get_size();
			residual_space *= (((double)size))/((double)(size+y.get_reduction()));
	    
#ifdef _DEBUG_IMPACT2
			std::cout << residual_space << ")\n";
#endif
	    
			++i;
		} 

#ifdef _DEBUG_IMPACT
		std::cout << " ==> impact[" << solver->variables[dec] << "] was " << impact[dec]  ;
#endif
	  
		impact[dec] = (((double)(tot_probes[dec]) * impact[dec]) + residual_space)/(double)(++tot_probes[dec]);
		
#ifdef _DEBUG_IMPACT
		std::cout << " now " << impact[dec] << std::endl;
#endif
		
		int offset = init_min[dec];
		int fact = factor[dec];
		if(dec_type == Decision::ASSIGNMENT) {
			dec_val -= offset;
			dec_val /= fact;
			
#ifdef _DEBUG_IMPACT
			std::cout << "i[" << x << "=" << dec_val << "]: " << value_weight[dec][dec_val]  << " -> ";
#endif
			
#ifdef _REAL_AVG_IMPACT
			value_weight[dec][dec_val] *= ((double)(value_visit[dec][dec_val]));
			value_weight[dec][dec_val] += (1.0 - residual_space);
			value_weight[dec][dec_val] /= (double)(++value_visit[dec][dec_val]);
#else
			value_weight[dec][dec_val] *= (alpha-1);
			value_weight[dec][dec_val] += (1.0 - residual_space);
			value_weight[dec][dec_val] /= alpha;			
#endif
			
#ifdef _DEBUG_IMPACT
			std::cout << value_weight[dec][dec_val] << std::endl;
#endif
			
		} else if(fact!=1) {
			int lb = x.get_min();
			int ub = x.get_max();
			int stag = (lb-offset)/fact;
			int endg = (ub-offset)/fact;
			int decg = (dec_val-offset)/fact;
			//int numg = (endg-stag+1);
			for(int g=stag; g<=endg; ++g) if(g!=decg) {
				
#ifdef _DEBUG_IMPACT
				std::cout << "i[" << x << "=(" << (g*fact+offset) << ", " << ((g+1)*fact+offset-1) << ")" << "]: " << value_weight[dec][g]  << " -> ";
#endif
			
#ifdef _REAL_AVG_IMPACT
				value_weight[dec][g] *= ((double)(value_visit[dec][g]));// * solver->parameters.activity_decay);
				value_weight[dec][g] += (1.0 - residual_space);
				value_weight[dec][g] /= (double)(++value_visit[dec][g]);
#else
				value_weight[dec][g] *= (alpha-1);
				value_weight[dec][g] += (1.0 - residual_space);
				value_weight[dec][g] /= alpha;			
#endif

#ifdef _DEBUG_IMPACT				
				std::cout << value_weight[dec][g] << std::endl;
#endif
			}
		} else {
			int vali, vnxt=x.get_min();
			do {
				vali = vnxt;
				vnxt = x.next(vali);
				vali -= offset;

#ifdef _DEBUG_IMPACT				
				std::cout << "i[" << x << "=" << (vali+offset) << "]: " << value_weight[dec][dec_val]  << " -> ";
#endif
				
#ifdef _REAL_AVG_IMPACT
				value_weight[dec][vali] *= ((double)(value_visit[dec][vali]));// * solver->parameters.activity_decay);
				value_weight[dec][vali] += (1.0 - residual_space);
				value_weight[dec][vali] /= (double)(++value_visit[dec][vali]);		
#else
				value_weight[dec][vali] *= (alpha-1);
				value_weight[dec][vali] += (1.0 - residual_space);
				value_weight[dec][vali] /= alpha;			
#endif
				
#ifdef _DEBUG_IMPACT				
				std::cout << value_weight[dec][vali] << std::endl;		
#endif
				vali += offset;
				
			} while(vali<vnxt);
		}
	  
		++num_probes[dec];
		variable_weight[dec] = avg_branches[dec] * impact[dec];
	  
#ifdef _DEBUG_IMPACT
		std::cout << " ===> total weight of " << solver->variables[dec] << " = " << avg_branches[dec] << " * " << impact[dec] << " = " << variable_weight[dec] << " \n";
#endif
	  		
	}
      
	left = 1;
}


void Mistral::ImpactManager::notify_backtrack() {
	// propagation produced a wipe-out
	// - check if it was after a left or a right branch
	// - find out what was the decision/refutation

	int dec{0};
	Variable x;

	//if(!solver->decisions.empty()) {
	if(left==1) {
		// left branch
		x = solver->decisions.back().var;
		dec = x.id();

#ifdef _DEBUG_IMPACT
		std::cout << " notified of backtrack after the " << num_probes[dec] << " left branch on " << solver->variables[dec] << "\n";
		std::cout << " ==> left-weight[" << solver->variables[dec] << "] was " << impact[dec] ;
#endif

		impact[dec] = (((double)(tot_probes[dec]) * impact[dec]))/(double)(++tot_probes[dec]);

#ifdef _DEBUG_IMPACT
		std::cout << " now " << impact[dec] << std::endl;
#endif

		variable_weight[dec] = avg_branches[dec] * impact[dec] ;
		++num_probes[dec];

#ifdef _DEBUG_IMPACT
		std::cout << " ===> total weight of " << solver->variables[dec] << " = " << avg_branches[dec] << " * " << impact[dec] << " = " << variable_weight[dec] << " \n";
#endif

	} else if(left==0) {
		// right branch
		x = solver->decisions.back(0).var;
		dec = x.id();
		

#ifdef _DEBUG_IMPACT
		double nbr = avg_branches[dec];
#endif

#ifdef _COUNT_BRANCHES
		avg_branches[dec] = (avg_branches[dec]*(double)(tot_fails[dec]) + (double)(num_probes[dec]))/((double)(++tot_fails[dec]));
#endif		
		
#ifdef _DEBUG_IMPACT
		std::cout << " notified of backtrack after a right branch on " << solver->variables[dec] << "\n";
		std::cout << " ==> num branches on " << solver->variables[dec] << " was " << nbr
			<< " now " << avg_branches[dec] << std::endl;
#endif

		variable_weight[dec] = avg_branches[dec] * impact[dec];

#ifdef _DEBUG_IMPACT
		std::cout << " ===> total weight of " << solver->variables[dec] << " = " << avg_branches[dec] << " * " << impact[dec] << " = " << variable_weight[dec] << " \n";
#endif

		num_probes[dec] = 0;
	} 

	int offset = init_min[dec];
	int fact = factor[dec];
	
	if(fact!=1) {
		int lb = x.get_min();
		int ub = x.get_max();
		int stag = (lb-offset)/fact;
		int endg = (ub-offset)/fact;
		//int decg = (dec_val-offset)/fact;
		//int numg = (endg-stag+1);
		for(int g=stag; g<=endg; ++g) 
		//if(g!=decg) 
		{
			
#ifdef _DEBUG_IMPACT	
			std::cout << "i[" << x << "=(" << (g*fact+offset) << ", " << ((g+1)*fact+offset-1) << ")" << "]: " << value_weight[dec][g]  << " -> ";
#endif
			
#ifdef _REAL_AVG_IMPACT
			value_weight[dec][g] *= ((double)(value_visit[dec][g]));
			value_weight[dec][g] += 1.0;
			value_weight[dec][g] /= (double)(++value_visit[dec][g]);
#else
			value_weight[dec][g] *= (alpha-1);
			value_weight[dec][g] += 1.0;
			value_weight[dec][g] /= alpha;			
#endif
			// value_weight[dec][g] *= ((double)(value_visit[dec][g]) * solver->parameters.activity_decay);
			// value_weight[dec][g] += 1.0;
			// value_weight[dec][g] /= (double)(++value_visit[dec][g]);

#ifdef _DEBUG_IMPACT			
			std::cout << value_weight[dec][g] << std::endl;
#endif
			
		}
	} else {
		int vali, vnxt=x.get_min();
		do {
			vali = vnxt;
			vnxt = x.next(vali);
			vali -= offset;

#ifdef _DEBUG_IMPACT			
			std::cout << "i[" << x << "=" << vali << "]: " << value_weight[dec][vali]  << " -> ";
#endif
			
			
#ifdef _REAL_AVG_IMPACT
			value_weight[dec][vali] *= ((double)(value_visit[dec][vali]));
			value_weight[dec][vali] += 1.0;
			value_weight[dec][vali] /= (double)(++value_visit[dec][vali]);
#else
			value_weight[dec][vali] *= (alpha-1);
			value_weight[dec][vali] += 1.0;
			value_weight[dec][vali] /= alpha;			
#endif
			// value_weight[dec][vali] *= ((double)(value_visit[dec][vali]) * solver->parameters.activity_decay);
			// value_weight[dec][vali] += 1.0;
			// value_weight[dec][vali] /= (double)(++value_visit[dec][vali]);

#ifdef _DEBUG_IMPACT			
			std::cout << value_weight[dec][vali] << std::endl;				
#endif
			
		} while((vali+offset)<vnxt);
	}


	left = 0;
}


std::ostream& Mistral::ImpactManager::display(std::ostream& os, const bool all) const {
	os << "impact";
  return os;
}



int Mistral::RealImpactManager::get_minweight_value(const Variable x) {
	
	int best_val = 0;
	int idx = x.id();
	int offset = init_min[idx];
	double *wgt = value_weight[idx];
	
	if(factor[idx]==1) {
		best_val = x.get_min();
		double min_weight = wgt[best_val-offset], aux_weight;
		int vali, vnxt=x.next(best_val);
		do {
			vali = vnxt;
			vnxt = x.next(vali);
			aux_weight = wgt[vali-offset]; 
			if(aux_weight < min_weight) {
				min_weight = aux_weight;
				best_val = vali;
			}
		} while(vali<vnxt);
	} else {
		int fact = factor[idx];
		int the_min = x.get_min();
		int the_max = x.get_max();	
		int best_group = (the_min-offset)/fact;
		int group = best_group;
		double min_weight = wgt[best_group];
		int lb = group*fact+offset;
		int ub = lb+fact-1;
		while(ub < the_max) {
			++group;
			lb += fact;
			ub += fact;
			
			if(x.intersect(lb, ub) && min_weight>wgt[group]) {
				min_weight = wgt[group];
				best_group = group;
			}
		}
		best_val = x.next(best_group*fact+offset-1);
	}
	
	return best_val;
}
int Mistral::RealImpactManager::get_maxweight_value(const Variable x) {
	
	int best_val = 0;
	int idx = x.id();
	int offset = init_min[idx];
	double *wgt = value_weight[idx];
	
	if(factor[idx]==1) {
		best_val = x.get_min();
		double max_weight = wgt[best_val-offset], aux_weight;
		int vali, vnxt=x.next(best_val);
		do {
			vali = vnxt;
			vnxt = x.next(vali);
			aux_weight = wgt[vali-offset]; 
			if(aux_weight > max_weight) {
				max_weight = aux_weight;
				best_val = vali;
			}
		} while(vali<vnxt);
	} else {
		int fact = factor[idx];
		int the_min = x.get_min();
		int the_max = x.get_max();	
		int best_group = (the_min-offset)/fact;
		int group = best_group;
		double max_weight = wgt[best_group];
		int lb = group*fact+offset;
		int ub = lb+fact-1;
		while(ub < the_max) {
			++group;
			lb += fact;
			ub += fact;
			
			if(x.intersect(lb, ub) && max_weight<wgt[group]) {
				max_weight = wgt[group];
				best_group = group;
			}
		}
		best_val = x.next(best_group*fact+offset-1);
	}
	
	return best_val;
}

void Mistral::RealImpactManager::notify_success() {
			
	// propagation went without wipe-out
	// - check if it was after a left or a right branch
	// - find out what was the decision/refutation
	int dec;
	int id;
	int i, n;
	Variable x, y;
	double residual_space;
	int size;
	if(solver->decisions.size>0) {
				
		i = solver->trail_.back(5), n=solver->saved_vars.size;

		// left branch
		Decision branch = solver->decisions.back();
		x = branch.var;
		dec = x.id();
		int dec_type = branch.type();
		int dec_val = branch.value();

		residual_space = 1.0;
		while(i<n) {	
			id = solver->saved_vars[i];
			y = solver->variables[id];
	    
#ifdef _DEBUG_IMPACT2
			std::cout << " -> " << y << y.get_domain() << " lost " << y.get_reduction() << " values (";
#endif
	    
			size = y.get_size();
			residual_space *= (((double)size))/((double)(size+y.get_reduction()));
	    
#ifdef _DEBUG_IMPACT2
			std::cout << residual_space << ")\n";
#endif
	    
			++i;
		} 

#ifdef _DEBUG_IMPACT
		std::cout << " [PROP] ==> impact[" << solver->variables[dec] << "] was " << variable_weight[dec] ;
#endif
	  
		
		int offset = init_min[dec];
		int fact = factor[dec];
		if(dec_type == Decision::ASSIGNMENT) {
			dec_val -= offset;
			dec_val /= fact;
			
#ifdef _DEBUG_IMPACT
			std::cout << " i[" << x << "=" << dec_val << "]: " << value_weight[dec][dec_val]  << " -> ";
#endif
			
			variable_weight[dec] -= value_weight[dec][dec_val];
			
#ifdef _REAL_AVG_IMPACT
			value_weight[dec][dec_val] *= ((double)(value_visit[dec][dec_val]));
			value_weight[dec][dec_val] += residual_space;
			value_weight[dec][dec_val] /= (double)(++value_visit[dec][dec_val]);
#else
			value_weight[dec][dec_val] *= (alpha-1);
			value_weight[dec][dec_val] += residual_space;
			value_weight[dec][dec_val] /= alpha;			
#endif
			
			variable_weight[dec] += value_weight[dec][dec_val];
			
#ifdef _DEBUG_IMPACT
			std::cout << value_weight[dec][dec_val] ;
#endif
			
		} else if(fact!=1) {
			int lb = x.get_min();
			int ub = x.get_max();
			int stag = (lb-offset)/fact;
			int endg = (ub-offset)/fact;
			int decg = (dec_val-offset)/fact;
			//int numg = (endg-stag+1);
			for(int g=stag; g<=endg; ++g) if(g!=decg) {
				
#ifdef _DEBUG_IMPACT
				std::cout << "\n    i[" << x << "=(" << (g*fact+offset) << ", " << ((g+1)*fact+offset-1) << ")" << "]: " << value_weight[dec][g]  ;
#endif
				
				variable_weight[dec] -= value_weight[dec][g];
			
#ifdef _REAL_AVG_IMPACT
				value_weight[dec][g] *= ((double)(value_visit[dec][g]));// * solver->parameters.activity_decay);
				value_weight[dec][g] += residual_space;
				value_weight[dec][g] /= (double)(++value_visit[dec][g]);
#else
				value_weight[dec][g] *= (alpha-1);
				value_weight[dec][g] += residual_space;
				value_weight[dec][g] /= alpha;			
#endif
				
				variable_weight[dec] += value_weight[dec][dec_val];

#ifdef _DEBUG_IMPACT				
				std::cout << " -> " << value_weight[dec][g] ;
#endif
			}
			
#ifdef _DEBUG_IMPACT
		std::cout << std::endl;
#endif

		} else {
			int vali, vnxt=x.get_min();
			do {
				vali = vnxt;
				vnxt = x.next(vali);
				vali -= offset;

#ifdef _DEBUG_IMPACT				
				std::cout << "\n    i[" << x << "=" << (vali+offset) << "]: " << value_weight[dec][dec_val]  ;
#endif
				
				variable_weight[dec] -= value_weight[dec][vali];
				
#ifdef _REAL_AVG_IMPACT
				value_weight[dec][vali] *= ((double)(value_visit[dec][vali]));// * solver->parameters.activity_decay);
				value_weight[dec][vali] += residual_space;
				value_weight[dec][vali] /= (double)(++value_visit[dec][vali]);		
#else
				value_weight[dec][vali] *= (alpha-1);
				value_weight[dec][vali] += residual_space;
				value_weight[dec][vali] /= alpha;			
#endif
				
				variable_weight[dec] += value_weight[dec][vali];
				
#ifdef _DEBUG_IMPACT				
				std::cout << " -> " << value_weight[dec][vali] ;
#endif
				vali += offset;
				
			} while(vali<vnxt);

#ifdef _DEBUG_IMPACT
		std::cout << std::endl;
#endif
		
		}
			
#ifdef _DEBUG_IMPACT
		std::cout << " ---> now " << variable_weight[dec] << std::endl;
#endif
	
	}
      
	left = 1;
}


void Mistral::RealImpactManager::notify_backtrack() {
	// propagation produced a wipe-out
	// - check if it was after a left or a right branch
	// - find out what was the decision/refutation

	int dec{0};
	Variable x;

	if(left==1) {
		// left branch
		x = solver->decisions.back().var;
		dec = x.id();
		
		std::cout << "BACKTRACK FROM A LEFT BRANCH?\n";

	} else if(left==0) {
		// right branch
		x = solver->decisions.back(0).var;
		dec = x.id();
		
		std::cout << "BACKTRACK FROM A RIGHT BRANCH?\n";
		
	} 
	
#ifdef _DEBUG_IMPACT
		std::cout << " ==> impact[" << solver->variables[dec] << "] was " << variable_weight[dec] ;
#endif

	int offset = init_min[dec];
	int fact = factor[dec];
	
	if(fact!=1) {
		int lb = x.get_min();
		int ub = x.get_max();
		int stag = (lb-offset)/fact;
		int endg = (ub-offset)/fact;
		//int decg = (dec_val-offset)/fact;
		//int numg = (endg-stag+1);
		for(int g=stag; g<=endg; ++g) 
		//if(g!=decg) 
		{
			
#ifdef _DEBUG_IMPACT	
			std::cout << "\n    i[" << x << "=(" << (g*fact+offset) << ", " << ((g+1)*fact+offset-1) << ")" << "]: " << value_weight[dec][g]  << " -> ";
#endif
			
			variable_weight[dec] -= value_weight[dec][g];
			
#ifdef _REAL_AVG_IMPACT
			value_weight[dec][g] *= ((double)(value_visit[dec][g]));
			// value_weight[dec][g] += 1.0;
			value_weight[dec][g] /= (double)(++value_visit[dec][g]);
#else
			value_weight[dec][g] *= (alpha-1);
			// value_weight[dec][g] += 1.0;
			value_weight[dec][g] /= alpha;			
#endif
			
			variable_weight[dec] += value_weight[dec][g];

#ifdef _DEBUG_IMPACT			
			std::cout << value_weight[dec][g] ;
#endif
			
		}

#ifdef _DEBUG_IMPACT
		std::cout << std::endl;
#endif
		
	} else {
		int vali, vnxt=x.get_min();
		do {
			vali = vnxt;
			vnxt = x.next(vali);
			vali -= offset;

#ifdef _DEBUG_IMPACT			
			std::cout << "\n    i[" << x << "=" << vali << "]: " << value_weight[dec][vali]  << " -> ";
#endif
			
			variable_weight[dec] -= value_weight[dec][vali];
			
#ifdef _REAL_AVG_IMPACT
			value_weight[dec][vali] *= ((double)(value_visit[dec][vali]));
			// value_weight[dec][vali] += 1.0;
			value_weight[dec][vali] /= (double)(++value_visit[dec][vali]);
#else
			value_weight[dec][vali] *= (alpha-1);
			// value_weight[dec][vali] += 1.0;
			value_weight[dec][vali] /= alpha;			
#endif
			
			variable_weight[dec] += value_weight[dec][vali];
			
			// value_weight[dec][vali] *= ((double)(value_visit[dec][vali]) * solver->parameters.activity_decay);
			// value_weight[dec][vali] += 1.0;
			// value_weight[dec][vali] /= (double)(++value_visit[dec][vali]);

#ifdef _DEBUG_IMPACT			
			std::cout << value_weight[dec][vali] ;
#endif
			
		} while((vali+offset)<vnxt);

#ifdef _DEBUG_IMPACT3
		std::cout << std::endl;
#endif
		
	}
	
#ifdef _DEBUG_IMPACT
		std::cout << " ---> now " << variable_weight[dec] << std::endl;
#endif


	left = 0;
}


std::ostream& Mistral::RealImpactManager::display(std::ostream& os, const bool all) const {
	os << "impact";
  return os;
}


std::ostream& Mistral::DecisionCountManager::display(std::ostream& os, const bool all) const {
	os << "dcm";
  return os;
}

std::ostream& Mistral::FailureCountManager::display(std::ostream& os, const bool all) const {
      
	int *all_variables = new int[variable_weight.size];
	int *all_constraints = new int[constraint_weight.size];


	int w, 
	xwidth; //, // = mylog10(solver->variables[variable_weight.size-1].id()),
	//cwidth; // = mylog10(solver->constraints[constraint_weight.size-1].id());
    
	for(unsigned int i=0; i<variable_weight.size; ++i) {
		all_variables[i] = i;

		// w = mylog10(variable_weight[i]);
		// if(w>xwidth) xwidth = w;

	}

	for(unsigned int i=0; i<constraint_weight.size; ++i) {
		all_constraints[i] = i;

		// w = mylog10(constraint_weight[i]);
		// if(w>cwidth) cwidth = w;

	}


	weight_sorting_array = variable_weight.stack_;
	qsort(all_variables, variable_weight.size, sizeof(int), decreasing_weight);

	weight_sorting_array = constraint_weight.stack_;
	qsort(all_constraints, constraint_weight.size, sizeof(int), decreasing_weight);

	os << " c variable weight: \n c id: ";
	for(unsigned int i=0; i<variable_weight.size; ++i) {
		if(all || solver->sequence.contain(all_variables[i])) {
			xwidth = mylog10(variable_weight[all_variables[i]]);
			w = mylog10(all_variables[i]);
			if(w>xwidth) xwidth = w;

			os << std::setw(xwidth) << all_variables[i] << " ";
		}
	}
	os << "\n c va: ";
	for(unsigned int i=0; i<variable_weight.size; ++i) {
		if(all || solver->sequence.contain(all_variables[i])) {
			xwidth = mylog10(variable_weight[all_variables[i]]);
			w = mylog10(all_variables[i]);
			if(w>xwidth) xwidth = w;
	  
			os << std::setw(xwidth) << variable_weight[all_variables[i]] << " ";
		}
	}
	os << "\n c constraint weight: \n c id: ";
	for(unsigned int i=0; i<constraint_weight.size; ++i) {
		xwidth = mylog10(constraint_weight[all_constraints[i]]);
		w = mylog10(all_constraints[i]);
		if(w>xwidth) xwidth = w;

		os << std::setw(xwidth) << all_constraints[i] << " ";
	}
	os << "\n c va: ";
	for(unsigned int i=0; i<constraint_weight.size; ++i) {
		xwidth = mylog10(constraint_weight[all_constraints[i]]);
		w = mylog10(all_constraints[i]);
		if(w>xwidth) xwidth = w;

		os << std::setw(xwidth) << constraint_weight[all_constraints[i]] << " ";
	}
	os << std::endl;

	delete [] all_constraints;
	delete [] all_variables;

	return os;
}    



std::ostream& Mistral::ConflictCountManager::display(std::ostream& os, const bool all) const {
      
	int *all_variables = new int[variable_weight.size];

	int w, 
	xwidth; //, // = mylog10(solver->variables[variable_weight.size-1].id()),
	//cwidth; // = mylog10(solver->constraints[constraint_weight.size-1].id());
    
	for(unsigned int i=0; i<variable_weight.size; ++i) {
		all_variables[i] = i;

		// w = mylog10(variable_weight[i]);
		// if(w>xwidth) xwidth = w;

	}

	weight_sorting_array = variable_weight.stack_;
	qsort(all_variables, variable_weight.size, sizeof(int), decreasing_weight);


	os << " c variable weight: \n c id: ";
	for(unsigned int i=0; i<variable_weight.size; ++i) {
		if(all || solver->sequence.contain(all_variables[i])) {
			xwidth = mylog10(variable_weight[all_variables[i]]);
			w = mylog10(all_variables[i]);
			if(w>xwidth) xwidth = w;

			os << std::setw(xwidth) << all_variables[i] << " ";
		}
	}
	os << "\n c va: ";
	for(unsigned int i=0; i<variable_weight.size; ++i) {
		if(all || solver->sequence.contain(all_variables[i])) {
			xwidth = mylog10(variable_weight[all_variables[i]]);
			w = mylog10(all_variables[i]);
			if(w>xwidth) xwidth = w;
	  
			os << std::setw(xwidth) << variable_weight[all_variables[i]] << " ";
		}
	}
	os << std::endl;

	delete [] all_variables;

	return os;
}    

#ifdef _ABS_VAL
int Mistral::PruningCountManager::get_minweight_value(const Variable x) {
	
	
	int best_val = 0;
	int idx = x.id();
	int offset = init_min[idx];
	double *wgt = value_weight[idx];
	
	// std::cout << "get min weight " << x << " in " << x.get_domain() << " " << factor[idx] << " " << offset << std::endl;
	
	
	if(factor[idx]==1) {
		best_val = x.get_min();
		double min_weight = wgt[best_val-offset], aux_weight;
		int vali, vnxt=x.next(best_val);
		
		// std::cout << "  " << best_val << " " << min_weight << std::endl;
		
		do {
			vali = vnxt;
			vnxt = x.next(vali);
			aux_weight = wgt[vali-offset]; 
			
			// std::cout << "  " << vali << " " << aux_weight << std::endl;
			
			if(aux_weight < min_weight) {
				min_weight = aux_weight;
				best_val = vali;
			}
		} while(vali<vnxt);
	} else {
		int fact = factor[idx];
		int the_min = x.get_min();
		int the_max = x.get_max();	
		int best_group = (the_min-offset)/fact;
		int group = best_group;
		double min_weight = wgt[best_group];
		int lb = group*fact+offset;
		int ub = lb+fact-1;
		
		// std::cout << "  [" << lb << "," << ub << "] " << min_weight << std::endl;
		
		while(ub < the_max) {
			++group;
			lb += fact;
			ub += fact;
			
			if(x.intersect(lb, ub)) {
				
				// std::cout << "  [" << lb << "," << ub << "] " << wgt[group] << std::endl;
				
				if(min_weight>wgt[group]) {
					min_weight = wgt[group];
					best_group = group;
				}
			}
		}
		best_val = x.next(best_group*fact+offset-1);
	}
	
	// std::cout << " return " << best_val << std::endl;
	
	return best_val;
}
int Mistral::PruningCountManager::get_maxweight_value(const Variable x) {
	
	int best_val = 0;
	int idx = x.id();
	int offset = init_min[idx];
	double *wgt = value_weight[idx];
	
	// std::cout << "get min weight " << x << " in " << x.get_domain() << " " << factor[idx] << " " << offset << std::endl;
	
	
	if(factor[idx]==1) {
		best_val = x.get_min();
		double max_weight = wgt[best_val-offset], aux_weight;
		int vali, vnxt=x.next(best_val);
		
		// std::cout << "  " << best_val << " " << min_weight << std::endl;
		
		do {
			vali = vnxt;
			vnxt = x.next(vali);
			aux_weight = wgt[vali-offset]; 
			
			// std::cout << "  " << vali << " " << aux_weight << std::endl;
			
			if(aux_weight > max_weight) {
				max_weight = aux_weight;
				best_val = vali;
			}
		} while(vali<vnxt);
	} else {
		int fact = factor[idx];
		int the_min = x.get_min();
		int the_max = x.get_max();	
		int best_group = (the_min-offset)/fact;
		int group = best_group;
		double max_weight = wgt[best_group];
		int lb = group*fact+offset;
		int ub = lb+fact-1;
		
		// std::cout << "  [" << lb << "," << ub << "] " << min_weight << std::endl;
		
		while(ub < the_max) {
			++group;
			lb += fact;
			ub += fact;
			
			if(x.intersect(lb, ub)) {
				
				// std::cout << "  [" << lb << "," << ub << "] " << wgt[group] << std::endl;
				
				if(max_weight < wgt[group]) {
					max_weight = wgt[group];
					best_group = group;
				}
			}
		}
		best_val = x.next(best_group*fact+offset-1);
	}
	
	// std::cout << " return " << best_val << std::endl;
	
	return best_val;
}
#endif



#ifdef _ABS_VAL	
void Mistral::PruningCountManager::notify_backtrack() {
	
	int sz = solver->sequence.size;
	double wu = solver->parameters.activity_increment;
	
	if(solver->decisions.empty())
		left = 0;
	
	
	if(left) {

		Decision branch = solver->decisions.back();
		Variable x = branch.var;
		int dec = x.id();
		int dec_type = branch.type();
		int dec_val = branch.value();
		int offset = init_min[dec];
		int fact = factor[dec];
		if(dec_type == Decision::ASSIGNMENT) {
			
#ifdef _DEBUG_ABS
			std::cout << "(d) i[" << x << "=" << dec_val << "]: " << value_weight[dec][dec_val]  << " -> ";
#endif
			
			dec_val -= offset;
			dec_val /= fact;
						
#ifdef _REAL_AVG_ACTIVITY
			value_weight[dec][dec_val] *= (double)(value_visit[dec][dec_val]);
			value_weight[dec][dec_val] += wu*(double)sz;
			value_weight[dec][dec_val] /= (double)(++value_visit[dec][dec_val]);
#else
			value_weight[dec][dec_val] *= (alpha-1);
			value_weight[dec][dec_val] += wu*(double)sz;
			value_weight[dec][dec_val] /= alpha;
#endif
			
#ifdef _DEBUG_ABS
			std::cout << value_weight[dec][dec_val*fact+offset] << "(" << sz << ")" << std::endl;
#endif
			
		} else if(fact!=1) {
			int lb = x.get_min();
			int ub = x.get_max();
			int stag = (lb-offset)/fact;
			int endg = (ub-offset)/fact;
			int decg = (dec_val-offset)/fact;
			//int numg = (endg-stag+1);
			for(int g=stag; g<=endg; ++g) if(g!=decg) {
				
#ifdef _DEBUG_ABS
				std::cout << "(dr) i[" << x << "=(" << (g*fact+offset) << ", " << ((g+1)*fact+offset-1) << ")" << "]: " << value_weight[dec][g]  << " -> ";
#endif
				
#ifdef _REAL_AVG_ACTIVITY
				value_weight[dec][g] *= (double)(value_visit[dec][g]);
				value_weight[dec][g] += wu*(double)sz;
				value_weight[dec][g] /= (double)(++value_visit[dec][g]);
#else				
				value_weight[dec][g] *= (alpha-1);
				value_weight[dec][g] += wu*(double)sz;
				value_weight[dec][g] /= alpha;
#endif				
				

#ifdef _DEBUG_ABS				
				std::cout << value_weight[dec][g] << std::endl;
#endif
			}
		} else {
			int vali, vnxt=x.get_min();
			// double domsize = (double)(x.get_size());
			do {
				vali = vnxt;
				vnxt = x.next(vali);
				vali -= offset;

#ifdef _DEBUG_ABS				
				std::cout << "(dr) i[" << x << "=" << (vali+offset) << "]: " << value_weight[dec][vali]  << " -> ";
#endif
				
#ifdef _REAL_AVG_ACTIVITY
				value_weight[dec][vali] *= (double)(value_visit[dec][vali]);
				value_weight[dec][vali] += wu*(double)sz;
				value_weight[dec][vali] /= (double)(++value_visit[dec][vali]);
#else				
				value_weight[dec][vali] *= (alpha-1);
				value_weight[dec][vali] += wu*(double)sz;
				value_weight[dec][vali] /= alpha;
#endif

#ifdef _DEBUG_ABS				
				std::cout << value_weight[dec][vali] << std::endl;		
#endif
				
				vali += offset;
				
			} while(vali<vnxt);
		}
	} else {
		
		Decision branch = solver->decisions.back(0);
		Variable x = branch.var;
		int dec = x.id();
		int offset = init_min[dec];
		int fact = factor[dec];
		if(fact!=1) {
			int lb = x.get_min();
			int ub = x.get_max();
			int stag = (lb-offset)/fact;
			int endg = (ub-offset)/fact;
			for(int g=stag; g<=endg; ++g) {
				
#ifdef _DEBUG_ABS
				std::cout << "(dr) i[" << x << "=(" << (g*fact+offset) << ", " << ((g+1)*fact+offset-1) << ")" << "]: " << value_weight[dec][g]  << " -> ";
#endif
				
#ifdef _REAL_AVG_ACTIVITY
				value_weight[dec][g] *= (double)(value_visit[dec][g]);
				value_weight[dec][g] += wu*(double)sz;
				value_weight[dec][g] /= (double)(++value_visit[dec][g]);
#else				
				value_weight[dec][g] *= (alpha-1);
				value_weight[dec][g] += wu*(double)sz;
				value_weight[dec][g] /= alpha;
#endif				
				

#ifdef _DEBUG_ABS				
				std::cout << value_weight[dec][g] << std::endl;
#endif
			}
		} else {
			int vali, vnxt=x.get_min();
			// double domsize = (double)(x.get_size());
			do {
				vali = vnxt;
				vnxt = x.next(vali);
				vali -= offset;

#ifdef _DEBUG_ABS				
				std::cout << "(dr) i[" << x << "=" << (vali+offset) << "]: " << value_weight[dec][vali]  << " -> ";
#endif
				
#ifdef _REAL_AVG_ACTIVITY
				value_weight[dec][vali] *= (double)(value_visit[dec][vali]);
				value_weight[dec][vali] += wu*(double)sz;
				value_weight[dec][vali] /= (double)(++value_visit[dec][vali]);
#else				
				value_weight[dec][vali] *= (alpha-1);
				value_weight[dec][vali] += wu*(double)sz;
				value_weight[dec][vali] /= alpha;
#endif

#ifdef _DEBUG_ABS				
				std::cout << value_weight[dec][vali] << std::endl;		
#endif
				
				vali += offset;
				
			} while(vali<vnxt);
		}
	}
	
	left = 0;

}
#endif


void Mistral::PruningCountManager::notify_success() {
	
	//display(std::cout, false);
	

	
	
	int id;
	int i = solver->trail_.back(5), n=solver->saved_vars.size;
	double max_weight = 0;
	
	int sz = (n-i);
	double wu = solver->parameters.activity_increment;
	
#ifdef _DEBUG_ABS
		if(_DEBUG_ABS>1) {
			std::cout << sz << " variables have been pruned\n";
		}
#endif
	

	while(++i<n) {	
		id = solver->saved_vars[i]; 
		
#ifdef _DEBUG_ABS
		if(_DEBUG_ABS>1) {
			std::cout << solver->variables[id] << ": " << variable_weight[id]  << " -> ";
		}
#endif
		
		variable_weight[id] += weight_unit;
		
#ifdef _DEBUG_ABS
		if(_DEBUG_ABS>1) {
			std::cout << variable_weight[id] << std::endl;
		}
#endif
		
		if(variable_weight[id] > max_weight)
			max_weight = variable_weight[id];
	}
	
	
	
	if(max_weight>threshold) {
	
#ifdef _DEBUG_ABS
		std::cout << "scaling down" << std::endl;
#endif
		
		double rfactor = std::min(solver->parameters.activity_increment/weight_unit, 1/max_weight);
		weight_unit *= rfactor;
	
		n = solver->variables.size;
		for(int i=0; i<n; ++i) {

#ifdef _DEBUG_ABS
			std::cout << variable_weight[i] << " -> ";
#endif

			variable_weight[i] *= rfactor;

#ifdef _DEBUG_ABS
			std::cout << variable_weight[i] << std::endl;
#endif

		}
	
	}
	

#ifdef _ABS_VAL	

	if(n_restart==(int)(solver->statistics.num_restarts)) {
		++n_restart;
		return;
	}


	if(solver->decisions.empty())
		left = 0;

	if(left) {
		Decision branch = solver->decisions.back();
		Variable x = branch.var;
		int dec = x.id();
		int dec_type = branch.type();
		int dec_val = branch.value();
		int offset = init_min[dec];
		int fact = factor[dec];
		if(dec_type == Decision::ASSIGNMENT) {
			
#ifdef _DEBUG_ABS
			std::cout << "(a) i[" << x << "=" << dec_val << "]: <- (" << value_weight[dec][dec_val]  << " * " << (alpha-1) << " + " << (wu*(double)sz) << ") / " << alpha << " = ";
#endif

			dec_val -= offset;
			dec_val /= fact;
						
#ifdef _REAL_AVG_ACTIVITY
			value_weight[dec][dec_val] *= (double)(value_visit[dec][dec_val]);
			value_weight[dec][dec_val] += wu*double(sz);
			value_weight[dec][dec_val] /= (double)(++value_visit[dec][dec_val]);
#else			
			value_weight[dec][dec_val] *= (alpha-1);
			value_weight[dec][dec_val] += wu*(double)sz;
			value_weight[dec][dec_val] /= alpha;
#endif			
			
#ifdef _DEBUG_ABS
			std::cout << value_weight[dec][dec_val*fact+offset] << std::endl;
#endif
			
		} else if(fact!=1) {
			int lb = x.get_min();
			int ub = x.get_max();
			int stag = (lb-offset)/fact;
			int endg = (ub-offset)/fact;
			int decg = (dec_val-offset)/fact;
			//int numg = (endg-stag+1);
			for(int g=stag; g<=endg; ++g) if(g!=decg) {
				
#ifdef _DEBUG_ABS
				std::cout << "(r) i[" << x << "=(" << (g*fact+offset) << ", " << ((g+1)*fact+offset-1) << ")" << "]: " << value_weight[dec][g]  << " -> ";
#endif
				
				// value_weight[dec][g] *= ((double)(value_visit[dec][g]) * solver->parameters.activity_decay);
				// value_weight[dec][g] += weight_unit*double(sz);
				// value_weight[dec][g] /= (double)(++value_visit[dec][g]);
				
#ifdef _REAL_AVG_ACTIVITY
				value_weight[dec][g] *= (double)(value_visit[dec][g]);
				value_weight[dec][g] += wu*double(sz);
				value_weight[dec][g] /= (double)(++value_visit[dec][g]);
#else			
				value_weight[dec][g] *= (alpha-1);
				value_weight[dec][g] += wu*(double)sz;
				value_weight[dec][g] /= alpha;
#endif				
				

#ifdef _DEBUG_ABS				
				std::cout << value_weight[dec][g] << std::endl;
#endif
			}
		} else {
			int vali, vnxt=x.get_min();
			// double domsize = (double)(x.get_size());
			do {
				vali = vnxt;
				vnxt = x.next(vali);
				vali -= offset;

#ifdef _DEBUG_ABS				
				std::cout << "(r) i[" << x << "=" << (vali+offset) << "]: " << value_weight[dec][vali]  << " -> ";
#endif
				
				// value_weight[dec][vali] *= ((double)(value_visit[dec][vali]) * solver->parameters.activity_decay);
				// value_weight[dec][vali] += weight_unit*double(sz);
				// value_weight[dec][vali] /= (double)(++value_visit[dec][vali]);
				
#ifdef _REAL_AVG_ACTIVITY
				value_weight[dec][vali] *= (double)(value_visit[dec][vali]);
				value_weight[dec][vali] += wu*double(sz);
				value_weight[dec][vali] /= (double)(++value_visit[dec][vali]);
#else			
				value_weight[dec][vali] *= (alpha-1);
				value_weight[dec][vali] += wu*(double)sz;
				value_weight[dec][vali] /= alpha;
#endif

#ifdef _DEBUG_ABS				
				std::cout << value_weight[dec][vali] << std::endl;		
#endif
				
				vali += offset;
				
			} while(vali<vnxt);
		}
		
	} else {
		
		Decision branch = solver->decisions.back(0);
		Variable x = branch.var;
		int dec = x.id();
		int offset = init_min[dec];
		int fact = factor[dec];

		if(fact!=1) {
					int lb = x.get_min();
					int ub = x.get_max();
					int stag = (lb-offset)/fact;
					int endg = (ub-offset)/fact;

					for(int g=stag; g<=endg; ++g) {
				
		#ifdef _DEBUG_ABS
						std::cout << "(r) i[" << x << "=(" << (g*fact+offset) << ", " << ((g+1)*fact+offset-1) << ")" << "]: " << value_weight[dec][g]  << " -> ";
		#endif
				
						// value_weight[dec][g] *= ((double)(value_visit[dec][g]) * solver->parameters.activity_decay);
						// value_weight[dec][g] += weight_unit*double(sz);
						// value_weight[dec][g] /= (double)(++value_visit[dec][g]);
				
		#ifdef _REAL_AVG_ACTIVITY
						value_weight[dec][g] *= (double)(value_visit[dec][g]);
						value_weight[dec][g] += wu*double(sz);
						value_weight[dec][g] /= (double)(++value_visit[dec][g]);
		#else			
						value_weight[dec][g] *= (alpha-1);
						value_weight[dec][g] += wu*(double)sz;
						value_weight[dec][g] /= alpha;
		#endif				
				

		#ifdef _DEBUG_ABS				
						std::cout << value_weight[dec][g] << std::endl;
		#endif
					}
				} else {
					int vali, vnxt=x.get_min();
					// double domsize = (double)(x.get_size());
					do {
						vali = vnxt;
						vnxt = x.next(vali);
						vali -= offset;

		#ifdef _DEBUG_ABS				
						std::cout << "(r) i[" << x << "=" << (vali+offset) << "]: " << value_weight[dec][vali]  << " -> ";
		#endif
				
						// value_weight[dec][vali] *= ((double)(value_visit[dec][vali]) * solver->parameters.activity_decay);
						// value_weight[dec][vali] += weight_unit*double(sz);
						// value_weight[dec][vali] /= (double)(++value_visit[dec][vali]);
				
		#ifdef _REAL_AVG_ACTIVITY
						value_weight[dec][vali] *= (double)(value_visit[dec][vali]);
						value_weight[dec][vali] += wu*double(sz);
						value_weight[dec][vali] /= (double)(++value_visit[dec][vali]);
		#else			
						value_weight[dec][vali] *= (alpha-1);
						value_weight[dec][vali] += wu*(double)sz;
						value_weight[dec][vali] /= alpha;
		#endif

		#ifdef _DEBUG_ABS				
						std::cout << value_weight[dec][vali] << std::endl;		
		#endif
				
						vali += offset;
				
					} while(vali<vnxt);
				}
		
		
	}
#endif	
		
	if(solver->parameters.activity_decay<1 && solver->parameters.activity_decay>0) 
		weight_unit /= solver->parameters.activity_decay;
	
	left = 1;
	
	// std::cout << 22 << std::endl;
}

std::ostream& Mistral::PruningCountManager::display(std::ostream& os, const bool all) const {
      
	int *all_variables = new int[variable_weight.size];


	int w, 
	xwidth; //, // = mylog10(solver->variables[variable_weight.size-1].id()),
	//cwidth; // = mylog10(solver->constraints[constraint_weight.size-1].id());
    
	for(unsigned int i=0; i<variable_weight.size; ++i) {
		all_variables[i] = i;

		// w = mylog10(variable_weight[i]);
		// if(w>xwidth) xwidth = w;

	}


	weight_sorting_array = variable_weight.stack_;
	qsort(all_variables, variable_weight.size, sizeof(int), decreasing_weight);

	os << " c variable weight: \n c id: ";
	for(unsigned int i=0; i<variable_weight.size; ++i) {
		if(all || solver->sequence.contain(all_variables[i])) {
			xwidth = mylog10(variable_weight[all_variables[i]]);
			w = mylog10(all_variables[i]);
			if(w>xwidth) xwidth = w;
	  
			os << std::setw(xwidth) << all_variables[i] << " ";
		}
	}
	os << "\n c va: ";
	for(unsigned int i=0; i<variable_weight.size; ++i) {
		if(all || solver->sequence.contain(all_variables[i])) {
			xwidth = mylog10(variable_weight[all_variables[i]]);
			w = mylog10(all_variables[i]);
			if(w>xwidth) xwidth = w;
	  
			os << std::setw(xwidth) << variable_weight[all_variables[i]] << " ";
		}
	}
	os << std::endl;

	delete [] all_variables;

	return os;
}    


Mistral::LearningActivityManager::LearningActivityManager(Solver *s) : solver(s) {
  weight_unit = solver->parameters.activity_increment;
  decay = solver->parameters.activity_decay;
  max_weight = std::numeric_limits<int>::max();
  
  var_activity.initialise(solver->variables.size, solver->variables.size, 0);
  lit_activity.initialise(2*solver->variables.size, 2*solver->variables.size, 0);
  
  int i = solver->constraints.size;
  Constraint *cons = solver->constraints.stack_;
  while(i--) {
    cons[i].initialise_activity(lit_activity.stack_, var_activity.stack_, weight_unit);
  }

  max_activity = 0;
  i = var_activity.size;
  while(i--) {
    if(var_activity[i]>max_activity)
      max_activity = var_activity[i];
  }

#ifdef _DEBUG_ACTIVITY
  for(int a=0; a<solver->variables.size; ++a) {
    std::cout << "init x" << a << " (" << lit_activity[NEG(a)] << "/" << lit_activity[POS(a)] << ")/" << var_activity[a]
	       << std::endl; 
  }
#endif
  
  //solver->lit_activity = lit_activity.stack_;
  //solver->var_activity = var_activity.stack_;
  
  solver->add((BacktrackListener*)this);
}

Mistral::LearningActivityManager::~LearningActivityManager() {
  solver->remove((BacktrackListener*)this);
}


void Mistral::LearningActivityManager::notify_backtrack() {

      // //std::cout << "d " << lit_activity.stack_ << " " << lit_activity[0] << " " << lit_activity[1] << std::endl;


#ifdef _DEBUG_ACTIVITY
  std::cout << "NOTIFY BACKTRACK!" << std::endl;
  std::cout << std::endl;
#endif

  int i;
  Literal q;

  weight_unit /= decay;
  if(max_weight - weight_unit <= max_activity) {
    // risk of double overflow

#ifdef _DEBUG_ACTIVITY
    std::cout << "\n RISK OF OVERFLOW (" << max_activity << " + " << weight_unit << " >= " << max_weight << ")\n";
#endif

    i=lit_activity.size;
    while(i--) lit_activity[i] /= max_activity;

    i=var_activity.size;
    while(i--) {
#ifdef _DEBUG_ACTIVITY
      std::cout << "x" << i << " (" << lit_activity[NEG(i)] << "/" << lit_activity[POS(i)] << ")/" << var_activity[i] << " -> " 
		<< (lit_activity[POS(i)]+lit_activity[NEG(i)]) << std::endl;
#endif
      var_activity[i] = lit_activity[POS(i)]+lit_activity[NEG(i)];
    }
   


    weight_unit = 1.0/decay;
    max_activity = 1.0;
  }
  i = solver->visited_literals.size;
  while(i--) {

#ifdef _DEBUG_ACTIVITY
    std::cout << "\n UPDATE ACTIVITY BECAUSE OF LITERAL " << q << ":\n";
#endif

    q = solver->visited_literals[i];
    Atom a = UNSIGNED(q);

#ifdef _DEBUG_ACTIVITY
    std::cout << "x" << a << " (" << lit_activity[NEG(a)] << "/" << lit_activity[POS(a)] << ")/" << var_activity[a]
     	      << " -> " ;     
#endif

    lit_activity[q] += weight_unit;
    var_activity[a] += weight_unit;
    if(var_activity[a] > max_activity)
      max_activity = var_activity[a];

    
#ifdef _DEBUG_ACTIVITY
    std::cout << " (" << lit_activity[NEG(a)] << "/" << lit_activity[POS(a)] << ")/" << var_activity[a]
     	      << std::endl ; 
#endif

  }
//   //std::cout << std::endl ; 


// #else

//       if(decay > 0 && decay < 1) {
//       	int i=var_activity.size;
//       	while(i--) {

// // #ifdef _DEBUG_ACTIVITY
// // 	  int a = i;
// // 	  std::cout << "decay x" << a << " (" << lit_activity[NEG(a)] << "/" << lit_activity[POS(a)] << ")/" << var_activity[a]
// // 		 << " -> ";
// // #endif

// 	  var_activity[i] *= decay;

// // #ifdef _DEBUG_ACTIVITY
// // 	  std::cout << " (" << lit_activity[NEG(a)] << "/" << lit_activity[POS(a)] << ")/" << var_activity[a]
// // 		    << std::endl;
// // #endif

// 	}
//       	// i=lit_activity.size;
//       	// while(i--) lit_activity[i] *= decay;
//       }

//   //solver->parameters.activity_increment *= (2.0-decay);

// #endif

}


std::ostream& Mistral::LearningActivityManager::display(std::ostream& os, const bool all) const {

  
  int *all_variables = new int[var_activity.size];


      int w, 
	xwidth; //, // = mylog10(solver->variables[var_activity.size-1].id()),
	//cwidth; // = mylog10(solver->constraints[constraint_weight.size-1].id());
    
      for(unsigned int i=0; i<var_activity.size; ++i) {
	all_variables[i] = i;

	// w = mylog10(var_activity[i]);
	// if(w>xwidth) xwidth = w;

      }


      weight_sorting_array = var_activity.stack_;
      qsort(all_variables, var_activity.size, sizeof(int), decreasing_weight);

      os << " c variable weight: \n c id: ";
      for(unsigned int i=0; i<var_activity.size; ++i) {

	if(!(i%1000)) {

	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = mylog10(var_activity[all_variables[i]])+7;
	  w = mylog10(all_variables[i]);
	  if(w>xwidth) xwidth = w;
	  
	  os << std::setw(xwidth) << all_variables[i] << " ";
	}

	}

      }
      os << "\n c va: ";
      for(unsigned int i=0; i<var_activity.size; ++i) {

	if(!(i%1000)) {

	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = mylog10(var_activity[all_variables[i]])+7;
	  w = mylog10(all_variables[i]);
	  if(w>xwidth) xwidth = w;
	  
	  long long int intval = (long long int)(100000 * var_activity[all_variables[i]]);
	  double outputval = ((double)intval)/100000;

	  os << std::setw(xwidth) << outputval << " ";
	}

	}

      }
     os << "\n c  0: ";
      for(unsigned int i=0; i<var_activity.size; ++i) {

	if(!(i%1000)) {

	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = mylog10(var_activity[2*all_variables[i]])+7;
	  w = mylog10(all_variables[i]);
	  if(w>xwidth) xwidth = w;
	  
	  long long int intval = (long long int)(100000 * lit_activity[2*all_variables[i]]);
	  double outputval = ((double)intval)/100000;

	  os << std::setw(xwidth) << outputval << " ";
	}

	}

      }
     os << "\n c  1: ";
      for(unsigned int i=0; i<var_activity.size; ++i) {

	if(!(i%1000)) {

	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = mylog10(lit_activity[2*all_variables[i]+1])+7;
	  w = mylog10(all_variables[i]);
	  if(w>xwidth) xwidth = w;
	  
	  long long int intval = (long long int)(100000 * lit_activity[2*all_variables[i]+1]);
	  double outputval = ((double)intval)/100000;

	  os << std::setw(xwidth) << outputval << " ";
	}

	}

      }
      os << std::endl;

    delete [] all_variables;


  return os;
}

Decision Mistral::MaxWeightValue::make(Variable x) {
	
	int best_val = 0;

	if(weight) {
		best_val = x.get_min();
		double *wgt = weight[x.id()];
		//double min_weight = weight[id_x][best_val]// [best_val][id_x]
		double max_weight = wgt[best_val]// [best_val][id_x]
			, aux_weight;
		int vali, vnxt=x.next(best_val);
		do {
			vali = vnxt;
			vnxt = x.next(vali);
			aux_weight = wgt[vali]; //weight[id_x][vali]; //weight[vali][id_x];
			if(aux_weight > max_weight) {
				max_weight = aux_weight;
				best_val = vali;
			}
		} while(vali<vnxt);
		
	} else {
		
		best_val = w_function->get_maxweight_value(x);
		
	}
  
	Decision d(x, Decision::ASSIGNMENT, best_val);
	
	return d;
}

Decision Mistral::MinWeightValue::make(Variable x) {
	
	int best_val = 0;

	if(weight) {
		best_val = x.get_min();

		double *wgt = weight[x.id()];
		//double min_weight = weight[id_x][best_val]// [best_val][id_x]
		double min_weight = wgt[best_val]// [best_val][id_x]
			, aux_weight;
		int vali, vnxt=x.next(best_val);
		do {
			vali = vnxt;
			vnxt = x.next(vali);
			aux_weight = wgt[vali]; //weight[id_x][vali]; //weight[vali][id_x];
			if(aux_weight < min_weight) {
				min_weight = aux_weight;
				best_val = vali;
			}
		} while(vali<vnxt);
		
	} else {
		
		best_val = w_function->get_minweight_value(x);
		
	}
  
	Decision d(x, Decision::ASSIGNMENT, best_val);
	
	return d;
}


std::ostream& operator<<(std::ostream& os, Mistral::MinDomain& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinMin& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxMax& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxRegret& x) {
  return x.display(os);
}

// std::ostream& operator<<(std::ostream& os, Mistral::MinDomainMaxDegree& x) {
//   return x.display(os);
// }

std::ostream& operator<<(std::ostream& os, Mistral::MinDomainOverDegree& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinDomainOverWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinDomainTimesWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinNeighborDomainOverWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinNeighborDomainOverNeighborWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::AnyValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MiddleValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MedianValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::RandomValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::HalfSplit& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::ReverseSplit& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::RandomSplit& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::RandomMinMax& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinWeightValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinWeightBound& x) {
  return x.display(os);
}

// std::ostream& operator<<(std::ostream& os, Mistral::Guided& x) {
//   return x.display(os);
// }

std::ostream& operator<<(std::ostream& os, Mistral::BoolMinWeightValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::BoolMaxWeightValue& x) {
  return x.display(os);
}



std::ostream& operator<<(std::ostream& os, Mistral::MinDomain* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinMin* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxMax* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxRegret* x) {
  return x->display(os);
}

// std::ostream& operator<<(std::ostream& os, Mistral::MinDomainMaxDegree* x) {
//   return x->display(os);
// }

std::ostream& operator<<(std::ostream& os, Mistral::MinDomainOverDegree* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinDomainOverWeight* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinDomainTimesWeight* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinNeighborDomainOverWeight* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinNeighborDomainOverNeighborWeight* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxWeight* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinWeight* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::AnyValue* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinValue* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxValue* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MiddleValue* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MedianValue* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::RandomValue* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::HalfSplit* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::ReverseSplit* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::RandomSplit* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::RandomMinMax* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinWeightValue* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinWeightBound* x) {
  return x->display(os);
}

// std::ostream& operator<<(std::ostream& os, Mistral::Guided* x) {
//   return x->display(os);
// }

std::ostream& operator<<(std::ostream& os, Mistral::BoolMinWeightValue* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::BoolMaxWeightValue* x) {
  return x->display(os);
}
