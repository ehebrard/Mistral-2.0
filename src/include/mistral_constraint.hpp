
/*
  Mistral 2.0 is a constraint satisfaction and optimisation library
  Copyright (C) 2009  Emmanuel Hebrard
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  The author can be contacted electronically at emmanuel.hebrard@gmail.com.
*/


/** \file mistral_constraint.hpp
    \brief Header for the Constraint library.
*/


#include <string>
#include <vector>

#include <mistral_variable.hpp>


#ifndef __CONSTRAINT_HPP
#define __CONSTRAINT_HPP


#define FILTER1( var, method )				    \
  event_type[(var)] = scope[(var)].method ;		     \
  if(event_type[(var)] == FAIL_EVENT) wiped = FAILURE(var);		\
  else if(event_type[(var)] != NO_EVENT && !changes.contain(var)) changes.add(var); 


#define FILTER2( evt, var, method )		\
  evt = scope[(var)].method ; \
  if(evt == FAIL_EVENT) wiped = FAILURE(var); \
  else if(evt != NO_EVENT && !changes.contain(var)) { \
  changes.add(var); \
  event_type[(var)] |= evt; \
  }



/**
 *Constraint:
 -propagator         (pointer to the constraint)
 -postponed          (should it be propagated on changed and/or postponed)


 *ConstraintImplementation: 
 -scope              (the variables)
 -id                 (an id)
 -trigger / index    (post/relax operations)

 +virtual propagate()
 +virtual propagate(const int changedIdx, const Event evt)
 +virtual check()

 *BinaryConstraintImplementation:

 *TernaryConstraintImplementation:
 
 *GlobalConstraintImplementation:
 */




namespace Mistral {

  /********************************************
   * ConstraintWrapper 
   ********************************************/
  /*! \class ConstraintWrapper
    \brief Wrapper for constraints. 
    
    Holds a pointer to a constraint and the index 
    (in the scope) of the variable that owns the list 
  */  

  
  class Variable;
  class ConstraintImplementation {
    
  public:
    
    /*!@name Parameters*/
    //@{
    /// Pointer to its solver
    Environment *solver;
    /// An unique id
    int id;
  
    // used to post/relax the constraints on the right triggers
    // the list on wich it should be added [among the 3 triggers of the variable]
    Vector< Trigger* > on;

    // the position of that constraint in the list
    int            *index;
    
    // what to put on the list: a pointer to the constraint and the rank of the variable 
    // (so that it knows who triggered it)
    Constraint      *self;

    // the variable to which the trigger applies to 
    Vector< Variable > _scope;

    
    //int arity;
    //int num_triggers;
    unsigned int type;
    //@}
    

    /*!@name Constructors*/
    //@{
    /// The _scope is build by copying an array "scp" containing "l" variables
    ConstraintImplementation();
    //ConstraintImplementation(const int a);
    virtual Constraint clone() = 0;
    virtual ~ConstraintImplementation();

    virtual void initialise() { type = get_type(); }
    virtual void initialise_vars(Solver*) = 0;

    virtual int get_type() = 0;
    virtual int postponed() = 0;
    virtual int idempotent() = 0;

    virtual void mark_domain() {}    

    virtual void desactivate(const int var) = 0;


   void un_post_from(const int var) {
     //std::cout << "relax [" << id << "] from " << _scope[var] << ": " << on[var] << " -> " ;      
      
      on[var]->relax(index[var]);
      index[var] = -1;    
      
      //std::cout << on[var]  << std::endl;
    }
    
    void un_relax_from(const int var) {
      // std::cout << var << std::endl;
      // std::cout << " urf post [" << id << "] on " << _scope[var] << ": " ;
      // std::cout.flush();
      // std::cout << on[var] << " -> " ;      
      
      index[var] = on[var]->post(self[var]);
      
      //std::cout << on[var]  << std::endl;
    }


    void trigger_on(const int t, Variable x) ;//{
    //   trigger.add(solver->constraint_trigger_graph[x.id()][t]);
    // }

    bool is_triggered_on(const int i, const int t);

    void initial_post(Solver *s);
   
    int get_trigger_type(const int i) ;
    void set_scope(const int i, Variable x);

    /*!@name Propagators*/
    //@{
    /*!
     * This methods returns 0 if the constraint is satisfied
     * by the combination and any positive value otherwise.  
     */
    virtual int check(const int*) const = 0;
    /*!
     *  This method is called when the domain of at least one   
     *  variable in _scope has been modified. The list of 
     *  changed variables is in 'modified', and the type
     *  of event in 'event_type'.
     *  Some domains in _scope may be reduced, and NULL is
     *  returned on success (there is still at least one consistent
     *  assignment) or a ptr to the wiped-out variable otherwise. 
     */
    virtual PropagationOutcome propagate() = 0; // { return NULL; }
    virtual PropagationOutcome checker_propagate() = 0;
    virtual PropagationOutcome bound_checker_propagate() = 0;
    virtual PropagationOutcome propagate(const int changed_idx, 
					 const Event evt) = 0;
    virtual PropagationOutcome checker_propagate(const int changed_idx, 
					 const Event evt) = 0;
    virtual PropagationOutcome bound_checker_propagate(const int changed_idx, 
					 const Event evt) = 0;
    //  { return CONSISTENT; }
    virtual bool absorb_negation(const int var) { return false; }
    virtual Constraint get_negation(const int var, Variable x) { return Constraint(); }
    virtual bool rewritable() { return false; }
    virtual RewritingOutcome rewrite() { return NO_EVENT; }
    virtual void consolidate() = 0; 
    virtual void consolidate_var(const int idx) = 0; 

    //virtual int get_backtrack_level();// {return solver->level-1;}
    //virtual Decision get_decision();// { return solver->decisions.back(0); }
    //@}


    inline Solver* get_solver() { return (Solver*)solver; }

    /*!@name Miscellanous methods*/
    //@{
    /// Print the constraint
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "c"; }
    //@}
  };





  class ConstraintTriggerArray {
    
  public:
    
    Trigger on[3];
    //Vector< ConstraintWrapper > range_trigger;
    //Vector< ConstraintWrapper > domain_trigger;

    ConstraintTriggerArray() ;
    ConstraintTriggerArray(const int size) ;
    void initialise(const int size) ;

    virtual ~ConstraintTriggerArray() ;

    inline int size() { return (on[0].size+on[1].size+on[2].size); }

    // each constraint keeps its index in the array it appears in
    // to remove: trigger.remove(index)
    // to add: index[i] = trigger[i].size; trigger.add(self[i]);

    std::ostream& display(std::ostream& os) const ;


    // class friend Iterator {
    
    // public :

    //   int T;
    //   int i;

    // };

  };


  /**
     
     POST/RELAX

     relax -> remove the constraint from all the active variables' lists
     relax_from(var) -> remove the constraint from var
     -> flag == active

     un_relax -> add the constraint to all the active vars' lists
     


   */

  template< int ARITY >
  class ValTuple {
  public:

    ValTuple() { std::fill(data, data+ARITY, NOVAL); }

    int data[ARITY];
    
    inline int& operator[](const int i) { return data[i]; }
    inline int operator[](const int i) const { return data[i]; }
  };



  template< int ARITY >
  class FixedArityConstraint : public ConstraintImplementation {

  public:

    int solution[ARITY];
    ValTuple<ARITY> *support[ARITY];
    Variable scope[ARITY];
    int active;


    FixedArityConstraint() : ConstraintImplementation() { 
      for(int i=0; i<ARITY; ++i) support[i] = NULL;
      //std::cout << "there" << std::endl;
      active = (1 << ARITY)-1; 
      //init_type(); //type = get_type();
    }

    FixedArityConstraint(Variable x, Variable y) : ConstraintImplementation() { 
      for(int i=0; i<ARITY; ++i) support[i] = NULL;
      //std::cout << "there" << std::endl;
      active = (1 << ARITY)-1; 
      //init_type(); //type = get_type();
      scope[0] = x;
      scope[1] = y;
    }

    FixedArityConstraint(Variable x, Variable y, Variable z) : ConstraintImplementation() { 
      for(int i=0; i<ARITY; ++i) support[i] = NULL;
      //std::cout << "there" << std::endl;
      active = (1 << ARITY)-1; 
      //init_type(); //type = get_type();
      scope[0] = x;
      scope[1] = y;
      scope[2] = z;
    }

    FixedArityConstraint(Vector< Variable >& scp) : ConstraintImplementation() { 
      for(int i=0; i<ARITY; ++i) support[i] = NULL;
      //std::cout << "there" << std::endl;
      active = (1 << ARITY)-1; 

      //std::cout << "there " << active << std::endl;

      //init_type(); //type = get_type();
      for(int i=0; i<ARITY; ++i) {
	scope[i] = scp[i];
      }
    }

    FixedArityConstraint(std::vector< Variable >& scp) : ConstraintImplementation() { 
      for(int i=0; i<ARITY; ++i) support[i] = NULL;
      //std::cout << "there" << std::endl;
      active = (1 << ARITY)-1; 
      //init_type(); //type = get_type();
      for(int i=0; i<ARITY; ++i) {
	scope[i] = scp[i];
      }
    }


    virtual void initialise_vars(Solver *s) {
      for(unsigned int i=0; i<ARITY; ++i) {
	scope[i].initialise(s, false);
      }
    }


    //int init_type() { std::cout << "here" << std::endl; type=get_type(); }
    virtual int get_type() {return 0;}
    virtual int postponed() {return 0;}
    //virtual void initialise() {type = get_type();}
    //virtual int idempotent() {return 0;}

    virtual void desactivate(const int var) { 
      int var_elt = (1 << var);
      if(active&var_elt) active ^= var_elt;
    }

    void initialise_supports() {
      int vari, vali, vnext, min_vari;
      for(vari=0; vari<ARITY; ++vari) {
	min_vari = scope[vari].get_min();
	
	support[vari] = new ValTuple<ARITY>[scope[vari].get_max() - min_vari + 1];
	support[vari] -= min_vari;
	
	vnext = scope[vari].get_min();
	do {
	  vali = vnext;
	  //support[vari][vali] = new int[ARITY];
	  support[vari][vali][vari] = vali;
	  vnext = scope[vari].next(vali);

	} while( vali < vnext );
      }
    }

    virtual void consolidate() {
      
      // std::cout << "consolidate " ;
      // display(std::cout);
      // std::cout << std::endl;

      for(unsigned int i=0; i<_scope.size; ++i) {
	_scope[i] = _scope[i].get_var();
      }
      
      for(unsigned int i=0; i<ARITY; ++i) {
	scope[i] = scope[i].get_var();
      }

      //std::cout << _scope << std::endl << scope[0] << " " << scope[1] << std::endl;

    }

   virtual void consolidate_var( const int idx ) {
      
     // std::cout << "consolidate " ;
     // display(std::cout);
     // std::cout << std::endl;
     
     _scope[idx] = _scope[idx].get_var();
     scope[idx] = scope[idx].get_var();
     
     //std::cout << _scope << std::endl << scope[0] << " " << scope[1] << std::endl;
   }


    void check_active() {
      for(int i=on.size; i--;) {
	if(active & (1 << i)) {
	  if(scope[i].is_ground()) {
	    std::cout << "Warning: " << scope[i] << " = " << scope[i].get_domain()
		      << " is ground and active!!" << std::endl;
	  }
	} else {
	  if(!(scope[i].is_ground())) {
	    std::cout << "Warning: " << scope[i] << " in " << scope[i].get_domain()
		      << " is not ground and no active!!" << std::endl;
	  }
	}
      }
    }
      
    void print_active() {
      int k=0, x = (active & (7<<k)) >> k;
      while(x) {
	
	if(x==7) std::cout << "{" << scope[0] << "," << scope[1] << "," << scope[2] << "}" ;
	else if(x==6) std::cout << "{"  << scope[1] << "," << scope[2] << "}" ;
	else if(x==5) std::cout << "{" << scope[0] << "," << scope[2] << "}" ;
	else if(x==3) std::cout << "{" << scope[0] << "," << scope[1] << "}" ;
	else if(x==4) std::cout << "{" << scope[2] << "}" ;
	else if(x==2) std::cout << "{" << scope[1] << "}" ;
	else if(x==1) std::cout << "{" << scope[0] << "}" ;
	else std::cout << "{}" ;
	
	k+=3;
	x = (active & (7<<k)) >> k;

      }
    }


    void post_on(const int var) {

      // std::cout << "po post " ;
      // display(std::cout);
      // std::cout << " on " << scope[var] << std::endl;


      if(scope[var].is_ground()) {
	active^=(1<<var);
	index[var] = -1;
      } else {

	//std::cout << "(yep)" << std::endl;

	solver->save( self[var] );
	un_relax_from(var);
      }
    }
  
    void post() {

      // std::cout << "p post " ;
      // display(std::cout);
      // std::cout << std::endl;

      active = (1 << ARITY)-1;
      for(int i=on.size; i--;) {
	post_on(i);
      }
    }

    void relax() {
      
      // std::cout << "active variables: ";
      // print_bitset(active, 0, std::cout);
      // std::cout << "( ";
      // for(int i=on.size; i--;) {
      // 	if(active & (1 << i)) {
      // 	  std::cout << scope[i] << " ";
      // 	}
      // }
      // std::cout << ")" << std::endl;

      for(int i=on.size; i--;) {
	relax_from(i);
      }    
    }

    void relax_from(const int var) {
      if(active & (1 << var)) {
      //if(index[var] >= 0) {
	solver->save( self[var] );
	un_post_from(var);
      }
    }    

    std::ostream& display(std::ostream& os) const {
      os << name() << "(" << scope[0] << ", " << scope[1];
      if(ARITY > 2)
	os << ", " << scope[2];
      os << ")";
      return os;
    }
  
  };


  class BinaryConstraint : public FixedArityConstraint<2> {
    
  public:

    BinaryConstraint() : FixedArityConstraint<2>() {}
    BinaryConstraint(Variable x, Variable y)  : FixedArityConstraint<2>(x,y) {}
    BinaryConstraint(Vector< Variable >& scp) : FixedArityConstraint<2>(scp) {}
    BinaryConstraint(std::vector< Variable >& scp) : FixedArityConstraint<2>(scp) {}
    virtual int get_type() { return BINARY|(IDEMPOTENT*idempotent()); }

    virtual PropagationOutcome propagate(); // { return NULL; }
    virtual PropagationOutcome bound_propagate(); // { return NULL; }
    virtual PropagationOutcome checker_propagate() { return BinaryConstraint::propagate(); }
    virtual PropagationOutcome bound_checker_propagate() { return BinaryConstraint::bound_propagate(); }
    virtual PropagationOutcome propagate(const int changed_idx, 
					 const Event evt);
    virtual PropagationOutcome bound_propagate(const int changed_idx, 
					 const Event evt);
    virtual PropagationOutcome checker_propagate(const int changed_idx, 
						 const Event evt) { return BinaryConstraint::propagate(changed_idx, evt); }
    virtual PropagationOutcome bound_checker_propagate(const int changed_idx, 
						       const Event evt) { return BinaryConstraint::propagate(changed_idx, evt); }
    bool find_support(const int revise_idx, const int vli);
    bool find_bound_support(const int revise_idx, const int vli);

   void restore(const int rtype) {
      int var = rtype&CTYPE;
      if(index[var]<0) {
	un_relax_from(var);
      } else {
	un_post_from(var);
      }
      active = 3;
    }  

    void trigger();
  };


  class TernaryConstraint : public FixedArityConstraint<3> {

  public:

    int lvl;

    TernaryConstraint() : FixedArityConstraint<3>() { lvl = 3; }
    TernaryConstraint(Variable x, Variable y, Variable z) : FixedArityConstraint<3>(x,y,z) { lvl = 3; }
    TernaryConstraint(Vector< Variable >& scp) : FixedArityConstraint<3>(scp) { lvl = 3; }
    TernaryConstraint(std::vector< Variable >& scp) : FixedArityConstraint<3>(scp) { lvl = 3; }
    virtual int get_type() {return TERNARY|(IDEMPOTENT*idempotent());}

    virtual PropagationOutcome checker_propagate() 
    { return TernaryConstraint::propagate(); }
    virtual PropagationOutcome bound_checker_propagate() { 
      return TernaryConstraint::bound_propagate(); 
    }
    virtual PropagationOutcome propagate(); // { return NULL; }
    virtual PropagationOutcome bound_propagate(); // { return NULL; }
    virtual PropagationOutcome propagate(const int changed_idx, 
					 const Event evt);
    virtual PropagationOutcome bound_propagate(const int changed_idx, 
					 const Event evt);
    virtual PropagationOutcome checker_propagate(const int changed_idx, 
						 const Event evt) 
    { return TernaryConstraint::propagate(changed_idx, evt); }
    virtual PropagationOutcome bound_checker_propagate(const int changed_idx, 
						       const Event evt) { 

      //std::cout << "constraint.cpp: bound checker propagate(evt)" << std::endl; 

      return TernaryConstraint::bound_propagate(changed_idx, evt); }
    bool find_support(const int revise_idx, const int vli);
    bool find_bound_support(const int revise_idx, const int vli);


    /*
    // the 3 first bits stand for the current set of active vars, 
    // the 3 next bits stand for the intermediate state 
    // lvl stands for the level at wich the intermediate


    -------- FIRST POSSIBILITY -----------
    // when assigning a new variable, 
    -- if lvl = -1 :  do the remove and change lvl (save)
    -- if lvl = s->level : do the remove
    -- if -1 < lvl < s->level : store the current state, do the remove (save)

    // when backtracking
    -- if s->level > lvl : copy the stored state
    -- otherwise : set to {0,1,2}

              -1|a    b
    {0,1,2} | {1,2} | {2}

             -1|a|a
    {0,1,2} | {1}  

    ------- SECOND POSSIBILITY -----------
    lvl store the current position where values should eb stored
    
    Whenever a variable is assigned, store the current state into the next slot and save.

    {1,2}{0,1,2}

    */



    bool assign(const int var) {

      // //if(id==36) {
      // std::cout << std::endl;
      // print_active();
      // std::cout << " Assign " << scope[var] << " " ;
      // //}
	
      bool ret_val = false;
      int elt = (1 << var);
      if(active & elt){
	int tmp = active&7;
	active <<= 3;
	active |= tmp;
	active ^= elt;

	if(active&448) ret_val = true;
	solver->save(Constraint(this, type|ACTIVITY));

	// if(active&64) ret_val = true;
	// else solver->save(Constraint(this, type|ACTIVITY));
      }

      // //if(id==36) {
      // print_active();
      // std::cout << std::endl;
      // //}

      return ret_val;
    }


    void restore(const int rtype) {

      // //if(id==36) {
      // //std::cout << std::endl;

      // // std::cout << scope[0] << " " << on[0] << std::endl;
      // // std::cout << scope[1] << " " << on[1] << std::endl;
      // // std::cout << scope[2] << " " << on[2] << std::endl;
      // print_active();
      // std::cout << " restore " ;
      // //}

      if(rtype&ACTIVITY) active >>= 3;	
      else {
      int var = rtype&CTYPE;

      //std::cout << var << " " << index[var] << " ";


      if(index[var]<0) {

	// std::cout << on[var] << " repost " ;
	// display(std::cout);
	// //std::cout << std::endl; 

	un_relax_from(var);

	// std::cout << " " << on[var] << " ";

      } else {
	un_post_from(var);
      }
      }
      // //if(id==36) {
      // print_active();
      // std::cout << std::endl ;

      // // std::cout << scope[0] << " " << on[0] << std::endl;
      // // std::cout << scope[1] << " " << on[1] << std::endl;
      // // std::cout << scope[2] << " " << on[2] << std::endl;
      // std::cout << std::endl;

      // //}

    }

    void update(const int changed_idx, const Event evt) {
      if(ASSIGNED(evt) && assign(changed_idx)) {

	
	// std::cout << " (" << ((active&7)/2) << ") ";
	// std::cout.flush();
	// std::cout << on[((active&7)/2)] << " relax " ;
	// display(std::cout) ;
	// std::cout.flush();


	relax_from((active&7)/2);

	// std::cout << " from " << scope[((active&7)/2)] 
	// 	  << on[((active&7)/2)] << std::endl;

      }
    }
    
    void trigger();
  };

  class GlobalConstraint : public ConstraintImplementation {

  public:

    Vector< Variable > scope;
    /// We use two lists so that the active constraint can add events to its own 
    /// list (events) without changing the one used during propagation (changes)
    IntStack changes; // this is the list that is accessible from a propagator
    IntStack events; // this is the list that collects the events
    /// The type of event for each modified variable
    Event *event_type;
    /// Set of non-ground variables
    ReversibleSet active;

    int priority;

    ////
    int   *solution;
    int ***supports;
    ////

  
    GlobalConstraint() : ConstraintImplementation() {}
    GlobalConstraint(Vector< Variable > scp);
    GlobalConstraint(std::vector< Variable > scp);
    GlobalConstraint(Variable* scp, const int n);
    void initialise();
    virtual int get_type() {return (PUSHED*pushed())|(POSTPONED*postponed())|(IDEMPOTENT*idempotent());};
    virtual int pushed() = 0;
    virtual int postponed() = 0;
    virtual int idempotent() = 0;
    virtual ~GlobalConstraint();

    virtual void initialise_vars(Solver*);

   virtual void desactivate(const int var) { 
     active.remove(var);
   }

    void trigger();

    virtual void consolidate() {
      for(unsigned int i=0; i<scope.size; ++i) {
	scope[i] = scope[i].get_var();
      }
     for(unsigned int i=0; i<_scope.size; ++i) {
	_scope[i] = _scope[i].get_var();
      }
    }

    virtual void consolidate_var(const int idx) {
      scope[idx] = scope[idx].get_var();
      _scope[idx] = _scope[idx].get_var();
    }
  
    /// An idempotent constraint should not be called on events triggered by itself.
    /// To forbid that, the lists 'events' and 'changes' are merged
    /// during its propagation, events will be added to the events list of the constraint 
    /// after the propagation, the lists events and changes are swapped back
    /// and the change list is cleared. When idempotent, since the two lists
    /// point to the same object, they are both cleared.
    inline void set_idempotent() {
      if(idempotent()) {
	events.size = 0;
	events.capacity = changes.capacity;
	events.list_ = changes.list_;
	events.index_ = changes.index_;
	events.start_ = NULL;
      } else {
	events.initialise(0, changes.capacity-1, false);
      }
    }


    // called before propagation, the events strored in the list 'events'
    // are copied onto the list 'changes'.
    inline void freeze() {
      // if the constraint is idempotent, the two lists point to the same
      // elements, we just set the size up to what it should be.
      changes.size = events.size;

      // otherwise,
      // before each propagation, the lists events and changes are swapped 
      if(changes.list_ != events.list_) {

	// make the changes list points to the events list, and clear the other
	events.size = 0;
	
	int *iaux = events.list_;
	events.list_ = changes.list_;
	changes.list_ = iaux;
	
	unsigned int *uaux = events.index_;
	events.index_ = changes.index_;
	changes.index_ = uaux;
      }
    }

    inline void defrost() {
      //      if(is_posted && active.size <= stress) {

      // #ifdef _DEBUG_CGRAPH
      // 	std::cout << " ---> relax " ;
      // 	display(std::cout);
      // 	std::cout << std::endl;
      // #endif

      // 	relax();
      //       }
      if(changes.list_ == events.list_)
	// if idempotent, clear the events list
	events.size = 0;	
    }

    inline void un_relax() {
      for(int n=active.size; --n;) {
	index[active[n]] = on[active[n]]->post(self[active[n]]);
      }
    }

    inline void post() {

      solver->save( Constraint(this, type|2) );
    
      active.fill();
      for(int i=on.size; --i;) {
	if(scope[i].is_ground()) {
	  active.reversible_remove(i);
	  index[i] = -1;
	} else index[i] = on[i]->post(self[i]);
      }
    }

    inline void un_post() {
      for(int i=on.size; --i;) {
	if(index[i]>=0) {
	  un_post_from(i);
	  // on[i]->relax(index[i]);
	  // index[i] = -1;
	}
      }
    }

    inline void relax() {
      solver->save( Constraint(this, type|1) );

      for(int i=on.size; --i;) {
	if(active.contain(i)) {
	  un_post_from(i);
	  // on[i]->relax(index[i]);
	  // index[i] = -1;
	}
      }    
    }

    virtual PropagationOutcome checker_propagate() { return GlobalConstraint::propagate(); }
    virtual PropagationOutcome bound_checker_propagate() { return GlobalConstraint::bound_propagate(); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome bound_propagate();
    virtual PropagationOutcome bound_propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome checker_propagate(const int changed_idx, 
						 const Event evt) 
    { return GlobalConstraint::propagate(changed_idx, evt); } 
    virtual PropagationOutcome bound_checker_propagate(const int changed_idx, 
						       const Event evt) 
    { return GlobalConstraint::bound_propagate(changed_idx, evt); } 
    virtual int get_backtrack_level();
    virtual Decision get_decision();// { return solver->decisions.back(0); }

    bool first_support(const int vri, const int vli);
    bool find_support(const int vri, const int vli);
    bool find_bound_support(const int vri, const int vli);

    inline void notify_other_event(const int var, 
				   const Mistral::Event evt) {

      

      if(ASSIGNED(evt) && active.contain(var)) {

	// std::cout << active << " -> " << scope[var] << " in " << scope[var].get_domain() 
	// 	  << " is assigned -> " ;
	
	active.reversible_remove(var);
	
	// std::cout << active << std::endl;

      }

      if(events.contain(var)) {
	event_type[var] |= evt;
      } else {
	events.add(var);
	event_type[var] = evt;
      }
    }

    inline void notify_first_event(const int var, 
				   const Mistral::Event evt) {

      if(ASSIGNED(evt) && active.contain(var)) {

	// std::cout << active << " -> " << scope[var] << " in " << scope[var].get_domain() 
	// 	  << " is assigned -> " ;

	active.reversible_remove(var);

	// std::cout << active << std::endl;

      }

      events.set_to(var);
      event_type[var] = evt;
    }


    std::ostream& display(std::ostream& os) const;  
  };




  // class ConstraintImplementation;
  // class ConstraintWrapper {




  /**********************************************
   * == Constraint
   **********************************************/ 
  /*! \class ConstraintEqual
    \brief  Binary equal Constraint (x0 == x1).
  */
  class ConstraintEqual : public BinaryConstraint {

  public:  
    /**@name Constructors*/
    //@{
    ConstraintEqual() : BinaryConstraint() {}
    ConstraintEqual(Variable x, Variable y) 
      : BinaryConstraint(x, y) {}
    ConstraintEqual(Vector< Variable >& scp) 
      : BinaryConstraint(scp[0], scp[1]) {}
    ConstraintEqual(std::vector< Variable >& scp) 
      : BinaryConstraint(scp[0], scp[1]) {}
    virtual Constraint clone() { return Constraint(new ConstraintEqual(scope[0], scope[1])// , type
						   ); }
    virtual void initialise();
    virtual int idempotent() { return 1;}
    virtual ~ConstraintEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check(const int* sol) const { return (sol[0] != sol[1]); }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    virtual bool rewritable() { return true; }
    virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "="; }
    //@}
  };

  /**********************************************
   * =/= Constraint
   **********************************************/ 
  /*! \class ConstraintNotEqual
    \brief  Binary not-equal Constraint (x0 =/= x1).
  */
  class ConstraintNotEqual : public BinaryConstraint {

  public:  
    /**@name Constructors*/
    //@{
    ConstraintNotEqual() : BinaryConstraint() {}
    ConstraintNotEqual(Variable x, Variable y)
      : BinaryConstraint(x, y) {}
    ConstraintNotEqual(Vector< Variable >& scp) 
      : BinaryConstraint(scp// [0], scp[1]
			 ) {}
    ConstraintNotEqual(std::vector< Variable >& scp) 
      : BinaryConstraint(scp[0], scp[1]) {}
    virtual Constraint clone() { return Constraint(new ConstraintNotEqual(scope[0], scope[1])// , type
						   ); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual ~ConstraintNotEqual() {}
    //@}


    /**@name Solving*/
    //@{
    virtual int check(const int* sol) const { return (sol[0] == sol[1]); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);  
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new ConstraintEqual( (var?scope[0]:x), (var?x:scope[1]) ) );
    }
    virtual bool absorb_negation(const int var) { 
      return (scope[0].get_min()==0 &&
	      scope[1].get_min()==0 &&
	      scope[0].get_max()==1 &&
	      scope[1].get_max()==1);
    }
    virtual bool rewritable() { return true; }
    virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "=/="; }
    //@}
  };




  /**********************************************
   * Equality Predicate
   **********************************************/
  /*! \class PredicateEqual
    \brief  Truth value of a binary equality ((x0 = x1) <-> y)
  */
  class PredicateEqual : public TernaryConstraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int spin;
    //@}

    /**@name Constructors*/
    //@{
    PredicateEqual() : TernaryConstraint() {}
    PredicateEqual(Variable x, Variable y, Variable z, const int sp=1) 
      : TernaryConstraint(x, y, z) { spin = sp; }
    PredicateEqual(Vector< Variable >& scp, const int sp=1) 
      : TernaryConstraint(scp) { spin = sp; }
    PredicateEqual(std::vector< Variable >& scp, const int sp=1) 
      : TernaryConstraint(scp) { spin = sp; }
    virtual Constraint clone() { return Constraint(new PredicateEqual(scope[0], scope[1], scope[2], spin)// , type
						   ); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual ~PredicateEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return((sol[0] == sol[1]) == (sol[2] ^ spin)); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual bool rewritable() { return true; }
    virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "=?="; }
    //@}
  };


  /**********************************************
   * ConstantEquality Predicate
   **********************************************/
  /*! \class PredicateConstantEqual
    \brief  Truth value of a binary equality ((x0 = x1) <-> y)
  */
  class PredicateConstantEqual : public BinaryConstraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int spin;
    int value;
    //@}

    /**@name Constructors*/
    //@{
    PredicateConstantEqual() : BinaryConstraint() {}
    PredicateConstantEqual(Variable x, Variable y, const int val=0, const int sp=1) 
      : BinaryConstraint(x,y) { value = val; spin = sp; }
    PredicateConstantEqual(Vector< Variable >& scp, const int val=0, const int sp=1) 
      : BinaryConstraint(scp) { value = val; spin = sp; }
    PredicateConstantEqual(std::vector< Variable >& scp, const int val=0, const int sp=1) 
      : BinaryConstraint(scp) { value = val; spin = sp; }
    virtual Constraint clone() { return Constraint(new PredicateConstantEqual(scope[0], scope[1], value, spin)// , type
						   ); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new PredicateConstantEqual( scope[0], x, value, !spin ) );
    }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual bool absorb_negation(const int var) { return var==1; }
    virtual ~PredicateConstantEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return((sol[0] == value) == (sol[1] ^ spin)); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "=k?"; }
    //@}
  };


  /**********************************************
   * Not Predicate
   **********************************************/
  /*! \class PredicateNot
    \brief  Truth value of a conjunction (!x0 <-> y)
  */
  class PredicateNot : public BinaryConstraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateNot() : BinaryConstraint() {}
    PredicateNot(Variable x, Variable y) 
      : BinaryConstraint(x,y) {}
    PredicateNot(Vector< Variable >& scp) 
      : BinaryConstraint(scp) {}
    PredicateNot(std::vector< Variable >& scp) 
      : BinaryConstraint(scp) {}
    virtual Constraint clone() { return Constraint(new PredicateNot(scope[0], scope[1])); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new ConstraintEqual( (var?scope[0]:x), (var?x:scope[1]) ) );
    }
    virtual void initialise();
    virtual int idempotent() { return 1; }
    virtual bool absorb_negation(const int var) { return true; }
    virtual ~PredicateNot() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0]>0) == sol[1]);
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual bool rewritable() { return true; }
    virtual RewritingOutcome rewrite();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "not?"; }
    //@}
  };


  /**********************************************
   * Neg Predicate
   **********************************************/
  /*! \class PredicateNeg
    \brief  Truth value of a conjunction (!x0 <-> y)
  */
  class PredicateNeg : public BinaryConstraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateNeg() : BinaryConstraint() {}
    PredicateNeg(Variable x, Variable y) 
      : BinaryConstraint(x,y) {}
    PredicateNeg(Vector< Variable >& scp) 
      : BinaryConstraint(scp) {}
    PredicateNeg(std::vector< Variable >& scp) 
      : BinaryConstraint(scp) {}
    virtual Constraint clone() { return Constraint(new PredicateNeg(scope[0], scope[1])); }
    virtual void initialise();
    virtual int idempotent() { return 1; }
    virtual ~PredicateNeg() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return(sol[0] != -sol[1]);
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "neg?"; }
    //@}
  };


  /**********************************************
   * Interval Membership Predicate
   **********************************************/
  /*! \class PredicateIntervalMember
    \brief  Truth value of a binary equality ((x0 = x1) <-> y)
  */
  class PredicateIntervalMember : public BinaryConstraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int spin;
    int lower_bound;
    int upper_bound;
    //@}

    /**@name Constructors*/
    //@{
    PredicateIntervalMember() : BinaryConstraint() {}
    PredicateIntervalMember(Variable x, Variable y, const int lb=-INFTY, const int ub=+INFTY, const int sp=1) 
      : BinaryConstraint(x,y) { spin = sp; lower_bound=lb; upper_bound=ub; }
    PredicateIntervalMember(Vector< Variable >& scp, const int lb=-INFTY, const int ub=+INFTY, const int sp=1) 
      : BinaryConstraint(scp) { spin = sp; lower_bound=lb; upper_bound=ub; }
    PredicateIntervalMember(std::vector< Variable >& scp, const int lb, const int ub, const int sp=1) 
      : BinaryConstraint(scp) { spin = sp; lower_bound=lb; upper_bound=ub; }
    virtual Constraint clone() { return Constraint(new PredicateIntervalMember(scope[0], scope[1], lower_bound, upper_bound, spin)); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new PredicateIntervalMember( scope[0], x, lower_bound, upper_bound, !spin ) );
    }
    virtual void initialise();
    virtual int idempotent() { return 1; }
    virtual bool absorb_negation(const int var) { return var==1; }
    virtual void mark_domain();
    virtual ~PredicateIntervalMember() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return(((sol[0]>=lower_bound) && (sol[0]<=upper_bound)) 
						       == (sol[1] ^ spin)); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "in[]?"; }
    //@}
  };


  /**********************************************
   * Set Membership Predicate
   **********************************************/
  /*! \class PredicateSetMember
    \brief  Truth value of a binary equality ((x0 = x1) <-> y)
  */
  class PredicateSetMember : public BinaryConstraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int spin;
    BitSet values;
    BitSet non_values;
    //@}

    /**@name Constructors*/
    //@{
    PredicateSetMember() : BinaryConstraint() {}
    PredicateSetMember(Vector< Variable >& scp, const int sp=1) 
      : BinaryConstraint(scp) { spin = sp; }
    PredicateSetMember(Variable x, Variable y, const BitSet& vals, const int sp=1) 
      : BinaryConstraint(x,y) { spin = sp; values=vals; }
    PredicateSetMember(Vector< Variable >& scp, const BitSet& vals, const int sp=1) 
      : BinaryConstraint(scp) { spin = sp; values=vals; }
    PredicateSetMember(std::vector< Variable >& scp, const BitSet& vals, const int sp=1) 
      : BinaryConstraint(scp) { spin = sp; values=vals; }
    virtual Constraint clone() { return Constraint(new PredicateSetMember(scope[0], scope[1], values, spin)); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new PredicateSetMember( scope[0], x, values, !spin ) );
    }
    virtual void initialise();
    virtual int idempotent() { return 1; }
    virtual bool absorb_negation(const int var) { return var==1; }
    virtual void mark_domain();
    virtual ~PredicateSetMember() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return((values.contain(sol[0]) == (sol[1] ^ spin))); }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "in{}?"; }
    //@}
  };


  /**********************************************
   * <= Constraint
   **********************************************/ 
  /*! \class ConstraintLess
    \brief  Binary Less Than Constraint (x0 + k <= x1).
  */
  class ConstraintLess : public BinaryConstraint {

  public: 
    /**@name Parameters*/
    //@{  
    int offset;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintLess() : BinaryConstraint() {}
    ConstraintLess(const int ofs=0) 
      : BinaryConstraint() { offset = ofs; }
    ConstraintLess(Variable x, Variable y, const int ofs=0) 
      : BinaryConstraint(x, y) { offset = ofs; }
    ConstraintLess(Vector< Variable >& scp, const int ofs=0) 
      : BinaryConstraint(scp) { offset = ofs; }
    ConstraintLess(std::vector< Variable >& scp, const int ofs=0) 
      : BinaryConstraint(scp) { offset = ofs; }
    virtual Constraint clone() { return Constraint(new ConstraintLess(scope[0], scope[1], offset)// , type
						   ); }
    virtual void initialise();
    virtual int idempotent() { return 1;}
    // virtual bool absorb_negation(const int var) { 
    //   return (offset = 0 &&
    // 	      scope[0].get_min()==0 &&
    // 	      scope[1].get_min()==0 &&
    // 	      scope[0].get_max()==1 &&
    // 	      scope[1].get_max()==1);
    // }
    virtual ~ConstraintLess() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return (sol[0]+offset > sol[1]); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "<"; }
    //@}
  };

  /**********************************************
   * <= Predicate
   **********************************************/ 
  /*! \class PredicateLess

    \brief  Truth value of a precedence ((x0 + k <= x1) <-> y)
  */
  class PredicateLess : public TernaryConstraint {

  public: 
    /**@name Parameters*/
    //@{  
    int offset;
    //@}

    /**@name Constructors*/
    //@{
    PredicateLess() : TernaryConstraint(), offset(0) {}
    PredicateLess(Variable x, Variable y, Variable z, const int ofs=0)
      : TernaryConstraint(x, y, z), offset(ofs) {}
    PredicateLess(Vector< Variable >& scp, const int ofs=0) 
      : TernaryConstraint(scp), offset(ofs) {}
    PredicateLess(std::vector< Variable >& scp, const int ofs=0) 
      : TernaryConstraint(scp), offset(ofs) {}
    virtual Constraint clone() { return Constraint(new PredicateLess(scope[0], scope[1], scope[2], offset)// , type
						   ); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new PredicateLess( scope[1], scope[0], x, 1-offset ) );
    }
    virtual void initialise();
    virtual int idempotent() { return 1;}
    virtual bool absorb_negation(const int var) { return var==2; }
    virtual ~PredicateLess() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return ((sol[0]+offset <= sol[1]) != sol[2]); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //@}
    
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "<?"; }
    //@}
  };



  /**********************************************
   * <=K Predicate
   **********************************************/ 
  /*! \class PredicateUpperBound

    \brief  Truth value of a precedence ((x0 <= k) <-> y)
  */
  class PredicateUpperBound : public BinaryConstraint {

  public: 

    int bound;

    /**@name Constructors*/
    //@{
    PredicateUpperBound() : BinaryConstraint() {}
    PredicateUpperBound(Variable x, Variable y, const int b=0) 
      : BinaryConstraint(x,y) { bound = b; }
    PredicateUpperBound(Vector< Variable >& scp, const int b=0) 
      : BinaryConstraint(scp) { bound = b; }
    PredicateUpperBound(std::vector< Variable >& scp, const int b=0) 
      : BinaryConstraint(scp) { bound = b; }
    virtual Constraint clone() { return Constraint(new PredicateUpperBound(scope[0], scope[1], bound)// , type
						   ); }
    virtual Constraint get_negation(const int var, Variable x);
    // { 
    //   return Constraint( new PredicateLowerBound( scope[0], x, bound-1 ) );
    // }
    virtual void initialise();
    virtual int idempotent() { return 1;}
    virtual bool absorb_negation(const int var) { return var==1; }
    virtual ~PredicateUpperBound() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return ((sol[0] <= bound) != sol[1]); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "<u?"; }
    //@}
  };


  /**********************************************
   * >=K Predicate
   **********************************************/ 
  /*! \class PredicateLowerBound

    \brief  Truth value of a precedence ((x0 >= k) <-> y)
  */
  class PredicateLowerBound : public BinaryConstraint {

  public: 

    int bound;

    /**@name Constructors*/
    //@{
    PredicateLowerBound() : BinaryConstraint() {}
    PredicateLowerBound(Variable x, Variable y, const int b=0) 
      : BinaryConstraint(x,y) { bound = b; }
    PredicateLowerBound(Vector< Variable >& scp, const int b=0) 
      : BinaryConstraint(scp) { bound = b; }
    PredicateLowerBound(std::vector< Variable >& scp, const int b=0) 
      : BinaryConstraint(scp) { bound = b; }
    virtual Constraint clone() { return Constraint(new PredicateLowerBound(scope[0], scope[1], bound)// , type
						   ); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new PredicateUpperBound( scope[0], x, bound+1 ) );
    }
    virtual void initialise();
    virtual int idempotent() { return 1;}
    virtual bool absorb_negation(const int var) { return var==1; }
    virtual ~PredicateLowerBound() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return ((sol[0] >= bound) != sol[1]); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return ">l?"; }
    //@}
  };



  /**********************************************
   * And Predicate
   **********************************************/
  /*! \class PredicateAnd
    \brief  Truth value of a conjunction ((x0 and x1) <-> y)
  */
  class PredicateAnd : public TernaryConstraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateAnd() : TernaryConstraint() {}
    PredicateAnd(Vector< Variable >& scp) 
      : TernaryConstraint(scp) {}
    PredicateAnd(Variable x, Variable y, Variable z) 
      : TernaryConstraint(x,y,z) {}
    PredicateAnd(std::vector< Variable >& scp) 
      : TernaryConstraint(scp) {}
    virtual Constraint clone() { return Constraint(new PredicateAnd(scope[0], scope[1], scope[2])); }
    // virtual Constraint get_negation(const int var, Variable x) { 
    //   return Constraint( new PredicateUpperBound( scope[0], x, bound+1 ) );
    // }
    virtual void initialise();
    virtual int idempotent() { return 1;}
    //virtual bool absorb_negation(const int var) { return true; }
    virtual ~PredicateAnd() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0] && sol[1]) != (sol[2])); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "and?"; }
    //@}
  };

  /**********************************************
   * And Constraint
   **********************************************/
  /*! \class ConstraintAnd
    \brief  Conjunction (x0 and x1)
  */
  class ConstraintAnd : public BinaryConstraint
  {

  public:
    /**@name Constructors*/
    //@{
    ConstraintAnd() : BinaryConstraint() {}
    ConstraintAnd(Vector< Variable >& scp) 
      : BinaryConstraint(scp) {}
    ConstraintAnd(Variable x, Variable y) 
      : BinaryConstraint(x,y) {}
    ConstraintAnd(std::vector< Variable >& scp) 
      : BinaryConstraint(scp) {}
    virtual Constraint clone() { return Constraint(new ConstraintAnd(scope[0], scope[1])); }
    virtual void initialise();
    virtual int idempotent() { return 1;}
    //virtual bool absorb_negation(const int var) { return true; }
    virtual ~ConstraintAnd() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return(!(sol[0] && sol[1])); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "and"; }
    //@}
  };


  /**********************************************
   * Or Predicate
   **********************************************/
  /*! \class PredicateOr
    \brief  Truth value of a disjunction ((x0 or x1) <-> y)
  */
  class PredicateOr : public TernaryConstraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateOr() : TernaryConstraint() {}
    PredicateOr(Vector< Variable >& scp) 
      : TernaryConstraint(scp) {}
    PredicateOr(Variable x, Variable y, Variable z) 
      : TernaryConstraint(x,y,z) {}
    PredicateOr(std::vector< Variable >& scp) 
      : TernaryConstraint(scp) {}
    virtual Constraint clone() { return Constraint(new PredicateOr(scope[0], scope[1], scope[2])); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( (var<2 ? 
			  new PredicateLess( x, (var?scope[0]:scope[1]), scope[2], 0 ) : 
			  NULL
			  //new PredicateNotOr(scope[0], scope[1], scope[2]) 
			  ) );
    }
    virtual void initialise();
    virtual int idempotent() { return 1;}
    virtual bool absorb_negation(const int var) { return true; }
    virtual ~PredicateOr() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0] || sol[1]) != (sol[2])); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "or?"; }
    //@}
  };

  /**********************************************
   * Or Constraint
   **********************************************/
  /*! \class ConstraintOr
    \brief  Disjunction (x0 or x1)
  */
  class ConstraintOr : public BinaryConstraint
  {

  public:
    /**@name Constructors*/
    //@{
    ConstraintOr(Vector< Variable >& scp) 
      : BinaryConstraint(scp) { }
    ConstraintOr(Variable x, Variable y) 
      : BinaryConstraint(x,y) {}
    ConstraintOr(std::vector< Variable >& scp) 
      : BinaryConstraint(scp) { }
    virtual Constraint clone() { return Constraint(new ConstraintOr(scope[0], scope[1])); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new ConstraintLess( x, (var?scope[0]:scope[1]), 0 ) );
    }
    virtual void initialise();
    virtual int idempotent() { return 1;}
    virtual bool absorb_negation(const int var) { return true; }
    virtual ~ConstraintOr() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return(!(sol[0] || sol[1])); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "or"; }
    //@}
  };


  /**********************************************
   * Lex Constraint
   **********************************************/
  /*! \class ConstraintLex
    \brief  Basic bloc of a lex lt/leq constraint

    let x0 and x1 being two cells with same rank on two rows/columns 
    and b0, b1 be two Boolean variables.
    This constraint ensures that 
       - x0 =/= x1 => b1
       - b0 < b1 => x0 < x1
       - b0 <= b1 
  */
  class ConstraintLex : public GlobalConstraint
  {

  public:
    /**@name Constructors*/
    //@{
    ConstraintLex() : GlobalConstraint() { priority=1; }
    ConstraintLex(Vector< Variable >& scp) 
      : GlobalConstraint(scp) { priority=1; }
    ConstraintLex(std::vector< Variable >& scp) 
      : GlobalConstraint(scp) { priority=1; }
    virtual Constraint clone() { return Constraint(new ConstraintLex(scope)); }
    virtual void initialise();
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual ~ConstraintLex() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return( ((!sol[2] && !sol[3]) > (sol[0] == sol[1])
	       || (sol[2] < sol[3]) > (sol[0] <  sol[1])
	       || sol[2] > sol[3])
	      ); 
    }
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "lex"; }
    //@}
  };


  /**********************************************
   * Offset Predicate
   **********************************************/
  /*! \class PredicateOffset
    \brief  Offset (x+k = y)
  */
  class PredicateOffset : public BinaryConstraint
  {

  public:

    /**@name Parameters*/
    //@{  
    int offset;
    //@}

    /**@name Constructors*/
    //@{
    PredicateOffset() : BinaryConstraint() {}
    PredicateOffset(Vector< Variable >& scp, const int ofs) 
      : BinaryConstraint(scp) { offset=ofs; }
    PredicateOffset(Variable x, Variable y, const int ofs) 
      : BinaryConstraint(x,y) { offset=ofs; }
    PredicateOffset(std::vector< Variable >& scp, const int ofs) 
      : BinaryConstraint(scp) { offset=ofs; }
    virtual Constraint clone() { return Constraint(new PredicateOffset(scope[0], scope[1], offset)); }
    virtual void initialise();
    virtual int idempotent() { return 1;}
    virtual ~PredicateOffset() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0]+offset) != sol[1]);
    }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "+k"; }
    //@}
  };


  /**********************************************
   * Addition Predicate
   **********************************************/
  /*! \class PredicateAdd
    \brief  Binary addition predicate (x0 + x1 = y)
  */
  class PredicateAdd : public TernaryConstraint
  {
    
  public:
    /**@name Constructors*/
    //@{
    PredicateAdd() : TernaryConstraint() {}
    PredicateAdd(Variable x, Variable y, Variable z)
      : TernaryConstraint(x, y, z) {}
    PredicateAdd(Vector< Variable >& scp)
      : TernaryConstraint(scp) {}
    PredicateAdd(std::vector< Variable >& scp)
      : TernaryConstraint(scp) {}
    virtual Constraint clone() { return Constraint(new PredicateAdd(scope[0], scope[1], scope[2])// , type
						   ); }
    virtual void initialise();
    virtual int idempotent() { return 0;}
    virtual ~PredicateAdd() {}
    //@}
    
    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return (sol[2] != (sol[0]+sol[1])); }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    virtual RewritingOutcome rewrite();
    //@}
    
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "+="; }
    //@}
  };


  /**********************************************
   * Factor Predicate
   **********************************************/
  /*! \class PredicateFactor
    \brief  Truth value of a conjunction (!x0 <-> y)
  */
  class PredicateFactor : public BinaryConstraint
  {

  public:

    /**@name Parameters*/
    //@{  
    int factor;
    //@}

    /**@name Constructors*/
    //@{
    PredicateFactor() : BinaryConstraint() {}
    PredicateFactor(Vector< Variable >& scp, const int fct) 
      : BinaryConstraint(scp) { factor=fct; }
    PredicateFactor(Variable x, Variable y, const int fct) 
      : BinaryConstraint(x,y) { factor=fct; }
    PredicateFactor(std::vector< Variable >& scp, const int fct) 
      : BinaryConstraint(scp) { factor=fct; }
    virtual Constraint clone() { return Constraint(new PredicateFactor(scope[0], scope[1], factor)); }
    virtual void initialise();
    virtual int idempotent() { return 1;}
    virtual ~PredicateFactor() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0]*factor) != sol[1]);
    }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "*k"; }
    //@}
  };



  /**********************************************
   * Multiplication Predicate
   **********************************************/
  /*! \class PredicateMul
    \brief  Binary mulition predicate (x0 + x1 = y)
  */
  class PredicateMul : public GlobalConstraint
  {
    
  public:

    int min_pos[3];
    int max_pos[3];
    int min_neg[3];
    int max_neg[3];
    int zero[3];

    /**@name Constructors*/
    //@{
    PredicateMul() : GlobalConstraint() { priority = 1; }
    PredicateMul(Vector< Variable >& scp) 
      : GlobalConstraint(scp) { priority = 1; }
    // PredicateMul(Variable x, Variable y, Variable z) 
    //   : GlobalConstraint(x,y,z) {}
    PredicateMul(std::vector< Variable >& scp) 
      : GlobalConstraint(scp) { priority = 1; }
    //virtual Constraint clone() { return Constraint(new PredicateMul(scope[0], scope[1], scope[2])); }
    virtual Constraint clone() { return Constraint(new PredicateMul(scope)); }
    virtual void initialise();
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual ~PredicateMul() {}
    //@}
    
    /**@name Solving*/
    //@{
    //
    PropagationOutcome revise_division(const int X, const int Y, const int Z);
    //
    PropagationOutcome revise_multiplication(const int X, const int Y, const int Z);
    //
    PropagationOutcome prune(const int lb_neg, 
			     const int ub_neg, 
			     const int lb_pos, 
			     const int ub_pos,
			     const bool pzero,
			     const int Z);

    virtual int check( const int* sol ) const { return (sol[2] != (sol[0]*sol[1])); }
    virtual PropagationOutcome propagate();
    virtual RewritingOutcome rewrite();
    //@}
    
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "*="; }
    //@}
  };




  /**********************************************
   * Disjunctive Constraint (implemented with two precedence constraints)
   **********************************************/ 
  /*! \class ConstraintDisjunctive
    \brief  Binary Disjunctive Constraint (x0 + p0 <= x1 || x1 + p1 <= x0).
  */
  class ConstraintDisjunctive : public BinaryConstraint {
    
  public: 
    /**@name Parameters*/
    //@{
    int processing_time[2];
    Constraint precedence[2];
    //@}
    
    /**@name Constructors*/
    //@{
    ConstraintDisjunctive() : BinaryConstraint() {}
    ConstraintDisjunctive(Variable x, Variable y, const int p0, const int p1); 
    ConstraintDisjunctive(Vector< Variable >& scp, const int p0, const int p1); 
    ConstraintDisjunctive(std::vector< Variable >& scp, const int p0, const int p1); 
    virtual Constraint clone() { 
      return Constraint(new ConstraintDisjunctive(scope[0], scope[1], 
						  processing_time[0], processing_time[1])); }
    virtual void initialise();
    virtual int idempotent() { return 1; }
    virtual ~ConstraintDisjunctive() {}
    //@}

    /**@name Solving*/
    //@{
    void decide(const int choice);
    virtual int check( const int* sol ) const { 
      return ((sol[0]+processing_time[0] > sol[1])
	      &&
	      (sol[1]+processing_time[1] > sol[0])); 
    }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    virtual void consolidate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "<>"; }
    //@}
  };


  /**********************************************
   * ReifiedDisjunctive Constraint (implemented with two precedence constraints)
   **********************************************/ 
  /*! \class ConstraintReifiedDisjunctive
    \brief  Binary ReifiedDisjunctive Constraint (x0 + p0 <= x1 || x1 + p1 <= x0).
  */
  class ConstraintReifiedDisjunctive : public TernaryConstraint {

  public: 
    /**@name Parameters*/
    //@{
    int processing_time[2];
    int *min_t0_ptr;
    int *max_t0_ptr;
    int *min_t1_ptr;
    int *max_t1_ptr;
    int *state;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintReifiedDisjunctive() : TernaryConstraint() {}
    ConstraintReifiedDisjunctive(Variable x, Variable y, Variable z, const int p0, const int p1); 
    ConstraintReifiedDisjunctive(Vector< Variable >& scp, const int p0, const int p1); 
    ConstraintReifiedDisjunctive(std::vector< Variable >& scp, const int p0, const int p1); 

    virtual Constraint clone() { return Constraint(new ConstraintReifiedDisjunctive(scope[0], scope[1], scope[2], processing_time[0], processing_time[1])); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new ConstraintReifiedDisjunctive( scope[0], scope[1], x, processing_time[0], processing_time[1] ) );
    }
    virtual void initialise();
    virtual ~ConstraintReifiedDisjunctive() {}
    virtual int idempotent() { return 1; }
    virtual bool absorb_negation(const int var) { return var==2; }
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      int ret_value = 0;

      

      if(sol[2]) { // sol[0] + p[0] <= sol[1]

	//std::cout << sol[0] << " + " << processing_time[0] << " <= " << sol[1] << std::endl;

	ret_value = (sol[0] + processing_time[0] > sol[1]);
      } else {

	//std::cout << sol[1] << " + " << processing_time[1] << " <= " << sol[0] << std::endl;

	ret_value = (sol[1] + processing_time[1] > sol[0]);
      }
      return ret_value;

      // return ( (!sol[2] && (sol[1] + processing_time[1] > sol[0])) 
      // 	       ||
      // 	       (sol[2] && (sol[0] + processing_time[0] > sol[1])) );
    }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //virtual void consolidate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "<>="; }
    //@}
  };



  /**********************************************
   * BoolSum Equal Constraint
   **********************************************/
  //  lb <= x1 + ... + xn <= ub
  /// Constraint on the value of the sum of a set of variables.
  class ConstraintBoolSumEqual : public GlobalConstraint {
 
  public:
    /**@name Parameters*/
    //@{
    int total;
    ReversibleNum<int> min_;
    ReversibleNum<int> max_;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintBoolSumEqual() : GlobalConstraint() { priority = 1; total = 0; }
    ConstraintBoolSumEqual(Vector< Variable >& scp, const int t);
    ConstraintBoolSumEqual(std::vector< Variable >& scp, const int t);
    virtual Constraint clone() { return Constraint(new ConstraintBoolSumEqual(scope, total)); }
    virtual void initialise();
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    //virtual bool absorb_negation(const int var) { return true; }
    virtual ~ConstraintBoolSumEqual();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "bsum=k"; }
    //@}
  };


  /**********************************************
   * BoolSum Interval Constraint
   **********************************************/
  //  lb <= x1 + ... + xn <= ub
  /// Constraint on the value of the sum of a set of variables.
  class ConstraintBoolSumInterval : public GlobalConstraint {
 
  public:
    /**@name Parameters*/
    //@{
    int lb;
    int ub;
    ReversibleNum<int> min_;
    ReversibleNum<int> max_;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintBoolSumInterval() : GlobalConstraint() { priority = 1; }
    ConstraintBoolSumInterval(Vector< Variable >& scp, const int l, const int u);
    ConstraintBoolSumInterval(std::vector< Variable >& scp, const int l, const int u);
    virtual Constraint clone() { return Constraint(new ConstraintBoolSumInterval(scope, lb, ub)); }
    virtual void initialise();
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    //virtual bool absorb_negation(const int var) { return true; }
    virtual ~ConstraintBoolSumInterval();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "bsum[]"; }
    //@}
  };


  /**********************************************
   * BoolSum  Predicate
   **********************************************/
  //  x1 + ... + xn-1 = xn
  /// predicate on the value of the sum of a set of variables.
  class PredicateBoolSum : public GlobalConstraint {
 
  public:
    /**@name Parameters*/
    //@{
    int lb;
    int ub;
    ReversibleNum<int> min_;
    ReversibleNum<int> max_;
    //@}

    /**@name Constructors*/
    //@{
    PredicateBoolSum() : GlobalConstraint() { priority = 1; }
    PredicateBoolSum(Vector< Variable >& scp);
    PredicateBoolSum(std::vector< Variable >& scp);
    PredicateBoolSum(Vector< Variable >& scp, Variable tot);
    PredicateBoolSum(std::vector< Variable >& scp, Variable tot);
    virtual Constraint clone() { return Constraint(new PredicateBoolSum(scope)); }
    virtual void initialise();
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual bool absorb_negation(const int var) { return true; }
    virtual ~PredicateBoolSum();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "bsum="; }
    //@}
  };


  /**********************************************
   * Min Predicate
   **********************************************/
  /*! \class PredicateMin
    \brief  Predicate on the mininum of a set of variables z = (x1 , ... , xn)
  */
  class PredicateMin : public GlobalConstraint {

  public:
    /**@name Parameters*/
    //@{ 
    int last_min;
    ReversibleSet candidates;
    //@}

    /**@name Constructors*/
    //@{
    PredicateMin(Vector< Variable >& scp);
    virtual ~PredicateMin();
    virtual Constraint clone() { return Constraint(new PredicateMin(scope)); }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual void initialise();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "min="; }
    //@}
  };


  /**********************************************
   * Max Predicate
   **********************************************/
  /*! \class PredicateMax
    \brief  Predicate on the maxinum of a set of variables z = (x1 , ... , xn)
  */
  class PredicateMax : public GlobalConstraint {

  public:
    /**@name Parameters*/
    //@{ 
    int last_max;
    ReversibleSet candidates;
    //@}

    /**@name Constructors*/
    //@{
    PredicateMax(Vector< Variable >& scp);
    virtual ~PredicateMax();
    virtual Constraint clone() { return Constraint(new PredicateMax(scope)); }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual void initialise();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "max="; }
    //@}
  };


  /**********************************************
   * Element Predicate
   **********************************************/
  /*! \class PredicateElement
    \brief  
  */
  class PredicateElement : public GlobalConstraint {

  public:
    /**@name Parameters*/
    //@{ 
    BitSet aux_dom;
    int offset;
    //@}

    /**@name Constructors*/
    //@{
    PredicateElement(Vector< Variable >& scp, const int o=0);
    PredicateElement(std::vector< Variable >& scp, const int o=0);
    virtual ~PredicateElement();
    virtual Constraint clone() { return Constraint(new PredicateElement(scope, offset)); }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual void initialise();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "[x]="; }
    //@}
  };


  /**********************************************
   * WeightedSum Constraint
   **********************************************/
  /*! \class ConstraintWeightedSum
    \brief  Constraint on a sum of variables (a1 * x1 + ... + an * xn = total)
  */
  class PredicateWeightedSum : public GlobalConstraint {

  public:
    /**@name Parameters*/
    //@{ 
    // Lower bound of the linear experssion
    int lower_bound;

    // Upper bound of the linear experssion
    int upper_bound;

    Vector< int > weight;
    // from index 0 to wpos (not included), the coefficients are all 1s
    int wpos;
    // from index wpos to wneg (not included), the coefficients are all >0
    int wneg;
    // from index wneg to size (not included), the coefficients are all <0

    // utils for the propagation
    int *lo_bound;
    int *up_bound;
    int *span;
    ReversibleNum<int> parity;
    //ReversibleIntStack unknown_parity;
    ReversibleSet unknown_parity;
    
    //@}

    /**@name Constructors*/
    //@{
    PredicateWeightedSum(Vector< Variable >& scp,
			 const int L=0, const int U=0);
    PredicateWeightedSum(Vector< Variable >& scp,
			 Vector< int >& coefs,
			 const int L=0, const int U=0);
    PredicateWeightedSum(std::vector< Variable >& scp,
			 std::vector< int >& coefs,
			 const int L=0, const int U=0);
    virtual ~PredicateWeightedSum();
    virtual Constraint clone() { return Constraint(new PredicateWeightedSum(scope, weight)); }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual void initialise();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "sum="; }
    //@}
  };


  /***********************************************
   * AllDifferent Constraint (forward checking).
   ***********************************************/
  /*! \class ConstraintCliqueNotEqual
    \brief  Clique of NotEqual Constraint.
  */
  class ConstraintCliqueNotEqual : public GlobalConstraint {

  public:
    
    /**@name Constructors*/
    //@{
    ConstraintCliqueNotEqual() : GlobalConstraint() { priority = 2; }
    ConstraintCliqueNotEqual(Vector< Variable >& scp);
    ConstraintCliqueNotEqual(std::vector< Variable >& scp);
    ConstraintCliqueNotEqual(Variable* scp, const int n);
    virtual Constraint clone() { return Constraint(new ConstraintCliqueNotEqual(scope)// , type
						   ); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual ~ConstraintCliqueNotEqual();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "{=/=}"; }
    //@}
    
  };



  /***********************************************
   * All Different Constraint (bounds consistency).
   ***********************************************/
  typedef struct {
    int min, max;		// start, end of interval
    int minrank, maxrank; // rank of min & max in bounds[] of an adcsp
  } Interval;
  
  /*! \class ConstraintAllDiff
    \brief  AllDifferent Constraint.

    Constraint of difference on a set of variables.
    Ensures that a set of variables are assigned distinct 
    values (only Bounds Consistency is implemented)
    The code is from Lopez-Ortiz, Quimper, Tromp and van Beek
  */
  class ConstraintAllDiff : public GlobalConstraint {
    
  private:
    /**@name Parameters*/
    //@{  
    //int *level;
    int lastLevel;
    int *t;		// tree links
    int *d;		// diffs between critical capacities
    int *h;		// hall interval links
    Interval *iv;
    Interval **minsorted;
    Interval **maxsorted;
    int *bounds;  // bounds[1..nb] hold set of min & max in the niv intervals
                  // while bounds[0] and bounds[nb+1] allow sentinels
    int nb;
    void sortit();
    int filterlower();
    int filterupper();
    void propagateValue();
    //@}

  public:
    /**@name Constructors*/
    //@{
    ConstraintAllDiff() : GlobalConstraint() { priority = 0; }
    ConstraintAllDiff(Vector< Variable >& scp);
    ConstraintAllDiff(std::vector< Variable >& scp);
    ConstraintAllDiff(Variable* scp, const int n);
    virtual void mark_domain();
    virtual Constraint clone() { return Constraint(new ConstraintAllDiff(scope)// , type
						   ); }
    virtual void initialise();
    virtual ~ConstraintAllDiff();
    //void init();
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    //virtual std::string getString() const ;
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "alldiff"; }
    //@}
  };




  std::ostream& operator<< (std::ostream& os, const Constraint& x);
  std::ostream& operator<< (std::ostream& os, const Constraint* x);

  std::ostream& operator<< (std::ostream& os, const ConstraintImplementation& x);
  std::ostream& operator<< (std::ostream& os, const ConstraintImplementation* x);


  std::ostream& operator<< (std::ostream& os, const ConstraintTriggerArray& x);
  std::ostream& operator<< (std::ostream& os, const ConstraintTriggerArray* x);


}

#endif //__CONSTRAINT_HPP




