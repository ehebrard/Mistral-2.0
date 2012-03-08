
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


  //class Constraint;
  // class Trigger : public Vector< Constraint > {

  // public:

  //   int post(Constraint ct) ;

  //   void relax(const int idx) ;

  // };


  
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
    virtual void initialise_vars(Solver*);

    virtual int get_type() = 0;
    virtual int postponed() = 0;
    virtual int idempotent() = 0;

    virtual void mark_domain() {}    


   void un_post_from(const int var) {
     //std::cout << "relax [" << id << "] from " << _scope[var] << ": " << on[var] << " -> " ;      
      
      on[var]->relax(index[var]);
      index[var] = -1;    
      
      //std::cout << on[var]  << std::endl;
    }
    
    void un_relax_from(const int var) {
      //std::cout << "post [" << id << "] on " << _scope[var] << ": " << on[var] << " -> " ;      
      
      index[var] = on[var]->post(self[var]);
      
      //std::cout << on[var]  << std::endl;
    }


    void trigger_on(const int t, Variable x) ;//{
    //   trigger.add(solver->constraint_trigger_graph[x.id()][t]);
    // }

    bool is_triggered_on(const int i, const int t);

    void initial_post(Solver *s);

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
    virtual PropagationOutcome propagate(const int changed_idx, 
					 const Event evt) = 0;
    //  { return CONSISTENT; }
    virtual PropagationOutcome rewrite() {}
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
  class FixedArityConstraint : public ConstraintImplementation {

  public:

    Variable scope[ARITY];
    int active;


    FixedArityConstraint() : ConstraintImplementation() { 
      //std::cout << "there" << std::endl;
      active = (1 << ARITY)-1; 
      //init_type(); //type = get_type();
    }

    FixedArityConstraint(Variable x, Variable y) : ConstraintImplementation() { 
      //std::cout << "there" << std::endl;
      active = (1 << ARITY)-1; 
      //init_type(); //type = get_type();
      scope[0] = x;
      scope[1] = y;
    }

    FixedArityConstraint(Variable x, Variable y, Variable z) : ConstraintImplementation() { 
      //std::cout << "there" << std::endl;
      active = (1 << ARITY)-1; 
      //init_type(); //type = get_type();
      scope[0] = x;
      scope[1] = y;
      scope[2] = z;
    }

    FixedArityConstraint(Vector< Variable >& scp) : ConstraintImplementation() { 
      //std::cout << "there" << std::endl;
      active = (1 << ARITY)-1; 

      //std::cout << "there " << active << std::endl;

      //init_type(); //type = get_type();
      for(int i=0; i<ARITY; ++i) {
	scope[i] = scp[i];
      }
    }

    FixedArityConstraint(std::vector< Variable >& scp) : ConstraintImplementation() { 
      //std::cout << "there" << std::endl;
      active = (1 << ARITY)-1; 
      //init_type(); //type = get_type();
      for(int i=0; i<ARITY; ++i) {
	scope[i] = scp[i];
      }
    }

    //int init_type() { std::cout << "here" << std::endl; type=get_type(); }
    virtual int get_type() {return 0;}
    virtual int postponed() {return 0;}
    //virtual void initialise() {type = get_type();}
    //virtual int idempotent() {return 0;}


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

    // void un_relax() {
    //   // std::cout << "UNRELAX: " ;
    //   // print_active();
    //   // std::cout << std::endl;

    //   int i=0, x = (active&7);
    //   while(x) {
    // 	if(x&1) index[i] = on[i]->post(self[i]);
    // 	++i;
    // 	x>>=1;
    //   }
    // }


    void post_on(const int var) {

      if(scope[var].is_ground()) {
	active^=(1<<var);
	index[var] = -1;
      } else {
	solver->save( self[var] );
	//index[var] = on[var]->post(self[var]);
	un_relax_from(var);
      }

    }
  
    void post() {

      // std::cout << "post [" << id << "]";
      // //display(std::cout);
      // std::cout << " (domains: " << scope[0].get_domain()
      // 		<< ", " << scope[1].get_domain() << " - " 
      // 		<< active << ")" << std::endl; 


    
      active = (1 << ARITY)-1;
      for(int i=on.size; i--;) {
	//for(int i=0; i<ARITY; ++i) {

	//std::cout << "add " << index[i] << " to " << on[i] << " -> ";
	post_on(i);

	
      }
    }

    // void un_post() {

    //   // std::cout << "un_post [" << id << "]";
    //   // //display(std::cout);
    //   // std::cout << " (domains: " << scope[0].get_domain()
    //   // 		<< ", " << scope[1].get_domain() << " - " 
    //   // 		<< active << ")" << std::endl; 

    //   for(int i=on.size; i--;) {
    // 	//for(int i=0; i<ARITY; ++i) {

    // 	//std::cout << "remove " << index[i] << " from " << on[i] << " -> ";
	
    // 	if(index[i]>=0) {
    // 	  on[i]->relax(index[i]);
    // 	  index[i] = -1;
    // 	}

    // 	//std::cout << on[i]  << std::endl;
    //   }

    // }


 

    void relax() {
      
      // std::cout << "relax [" << id << "]";
      // //display(std::cout);
      // std::cout << " (domains: " << scope[0].get_domain()
      // 		<< ", " << scope[1].get_domain() << " - " 
      // 		<< active << ")" << std::endl; 


      //solver->save( Constraint(this, type|1) );
      
      for(int i=on.size; i--;) {

	//std::cout << active << "&" << (1 << i) << "?" << std::endl;

	relax_from(i);
	// //for(int i=0; i<ARITY; ++i) {
	// if(active & (1 << i)) {

	//   std::cout << "remove " << index[i] << " from " << on[i] << " -> ";
	//   solver->save( Constraint(this, type|i) );

	//   //on[i]->relax(index[i]);
	//   //index[i] = -1;
	//   un_post_from(i);
	  

	//   std::cout << on[i]  << std::endl;

	// }

	//std::cout << std::endl;
      }    
    }

    void relax_from(const int var) {

   

      // std::cout << "relax " ;
      // display(std::cout) ;
      // std::cout << " from " << _scope[var] << " in " << _scope[var].get_domain() 
      // 		<< " (" << _scope[1-var] << " in " << _scope[1-var].get_domain()
      // 		<< ") - " << active << std::endl;

      // std::cout << on[var]->size << std::endl;

      // std::cout << on[var] << std::endl;

      if(active & (1 << var)) {

	solver->save( self[var] );

	un_post_from(var);
	// on[var]->relax(index[var]);
	// index[var] = -1;
      }

      //std::cout << on[var] << std::endl;

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


   void restore(const int rtype) {

     // std::cout << "restore [" << id << "]" ;
     // print_active();
     // std::cout << std::endl
     // 	       << scope[0] << " in " << scope[0].get_domain() 
     // 	       << ": " << on[0]
     // 	       << std::endl
     // 	       << scope[1] << " in " << scope[1].get_domain() 
     // 	       << ": " << on[1]
     // 	       << std::endl;
     
      int var = rtype&CTYPE;
      if(index[var]<0) {
	un_relax_from(var);
	//active |= (1 << var);
      } else {
	un_post_from(var);
	//active |= (1 << var);
      }

      // if(rtype&1) un_relax();
      // else un_post();
      // // a binary constraint is relaxed as soon as the first variable get assigned
      // // so active is to be set back to {0,1} anyway
      active = 3;

      // print_active();
      // std::cout << std::endl << scope[0] << " in " << scope[0].get_domain() 
      // 		<< ": " << on[0]
      // 		<< std::endl
      // 		<< scope[1] << " in " << scope[1].get_domain() 
      // 		<< ": " << on[1]
      // 		<< std::endl << std::endl;
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
      
      // std::cout << std::endl;
      // print_active();
      // std::cout << " Assign " << scope[var] << " " ;

	
      bool ret_val = false;
      int elt = (1 << var);
      if(active & elt){
	int tmp = active&7;
	active <<= 3;
	active |= tmp;
	active ^= elt;
	if(active&64) ret_val = true;
	else solver->save(Constraint(this, type));
      }

      // print_active();
      // std::cout << std::endl;

      return ret_val;
    }


    void restore(const int rtype) {


      if(rtype&ACTIVITY) active <<= 3;	

      int var = rtype&CTYPE;
      if(index[var]<0) {
	un_relax_from(var);
      } else {
	un_post_from(var);
      }


      // //std::cout << "Restore " ;

      // if(rtype&2) {
	
      // 	//std::cout << " un-post" << std::endl;

      // 	un_post();
      // } else {
      // 	active <<= 3;	

      // 	//print_active();

      // 	if(rtype&1) {

      // 	  //std::cout << " un-relax" ;
      // 	  un_relax();
      // 	}

      // 	//std::cout << std::endl;
      // }
      // // a binary constraint is relaxed as soon as the first variable get assigned
      // // so active is to be set back to {0,1} anyway
    }

    void update(const int changed_idx, const Event evt) {
      if(ASSIGNED(evt) && assign(changed_idx)) {

  // std::cout << "RELAXFROM: " ;
  //     print_active() ;
  //     std::cout << std::endl;


	relax_from(active/2);
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
	  active.remove(i);
	  index[i] = -1;
	} else index[i] = on[i]->post(self[i]);
      }
    }

    inline void un_post() {
      for(int i=on.size; --i;) {
	if(index[i]>=0) {
	  on[i]->relax(index[i]);
	  index[i] = -1;
	}
      }
    }

    inline void relax() {
      solver->save( Constraint(this, type|1) );

      for(int i=on.size; --i;) {
	if(active.contain(i)) {
	  on[i]->relax(index[i]);
	  index[i] = -1;
	}
      }    
    }

    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt) {}
    virtual int get_backtrack_level();
    virtual Decision get_decision();// { return solver->decisions.back(0); }

    bool first_support(const int vri, const int vli);
    bool find_support(const int vri, const int vli);

    inline void notify_other_event(const int var, 
				   const Mistral::Event evt) {

      if(ASSIGNED(evt) && active.contain(var)) {
	active.remove(var);
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
	active.remove(var);
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
    virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "not_equal"; }
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
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "not_equal"; }
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
    virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "are_equal"; }
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
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual ~PredicateConstantEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return((sol[0] == value) == (sol[1] ^ spin)); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "equal_to"; }
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
    virtual ~ConstraintLess() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return (sol[0]+offset > sol[1]); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "less_than"; }
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
    virtual void initialise();
    virtual int idempotent() { return 1;}
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
    virtual std::string name() const { return "less_than"; }
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
    virtual void initialise();
    virtual int idempotent() { return 1;}
    virtual ~PredicateUpperBound() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return ((sol[0] <= bound) != sol[1]); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "upper_bound"; }
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
    virtual void initialise();
    virtual int idempotent() { return 1;}
    virtual ~PredicateLowerBound() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return ((sol[0] >= bound) != sol[1]); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "lower_bound"; }
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
    virtual int idempotent() { return 1;}
    virtual ~PredicateAdd() {}
    //@}
    
    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return (sol[2] != (sol[0]+sol[1])); }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome rewrite();
    //@}
    
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "plus"; }
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
    //virtual PropagationOutcome rewrite();
    virtual void consolidate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "disjunctive"; }
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
    virtual void initialise();
    virtual ~ConstraintReifiedDisjunctive() {}
    virtual int idempotent() { return 1; }
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
    //virtual PropagationOutcome rewrite();
    //virtual void consolidate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "r-disjunctive"; }
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
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "clique_ne"; }
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
    //virtual PropagationOutcome rewrite();
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




