
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


/*! \file mistral_backtrack.hpp
    \brief Header for the reversible structures
*/


#ifndef _MISTRAL_BACKTRACK_HPP
#define _MISTRAL_BACKTRACK_HPP


#include <mistral_global.hpp>


namespace Mistral {


  /*!
  A reversible structure keeps a pointer to an
  "Environment" (i.e., Solver) that manages the
  backtracking process.
  Second, it implements a virtual function restore
  that "undo" the last change. 
  Most reversible structures work this way:
  A Vector<int> trail_ is used to encode the
  delta-information used to undo. The last
  integer stored on the trail_ Vector is the
  value of solver->level when the last change
  occured. It is used, when changing the object,
  to decide if one should "save" the current state,
  or simply replace the current value.

  It is still undecided if this process can be 
  made generic enough so that we can make 
  restore() a static method.
 */


  /*! \class Environment
    \brief The minimal structures used to control the backtracking process
  */
  /********************************************
   * Environement Objects
   ********************************************/

  class Constraint;
  class Trigger : public Vector< Constraint > {

  public:
    
    int post(Constraint ct) ;

    void relax(const int idx) ;

  };


  template<class T> class ReversibleNum;
  class ReversibleSet;
  class Reversible;
  class Variable;
  class Decision;
  class Solver;
  class ConstraintImplementation;
  class Constraint  {
    
  public:
    
    ConstraintImplementation* propagator;
    unsigned int data;
    
    Constraint() {propagator = NULL; data=0;}
    Constraint(ConstraintImplementation* p) ;
    Constraint(ConstraintImplementation* p, const int t) { propagator = p; data = t; }
    void initialise(Solver*);
    virtual ~Constraint() {}
    
    
    inline int  index()      const { return data&CTYPE; }
    inline bool postponed()  const { return (data&POSTPONED); }
    inline bool pushed()     const { return (data&PUSHED); }
    inline bool global()     const { return !(data&0xc0000000); }
    inline bool binary()     const { return (data&BINARY); }
    inline bool ternary()    const { return (data&TERNARY); }
    inline bool idempotent() const { return (data&TERNARY); }

    inline void clear() { propagator = NULL; data = 0; }
    inline bool empty() const { return !propagator; }

    int id() const ;
    int priority() const ;
    void set_index(const int idx);
    void set_rank(const int idx);
    void set_id(const int idx);
    void post(Solver*);
    void awaken();
    void trigger();
    int check(int* sol);
    void consolidate();
    void consolidate_var();
    void re_link_to(Trigger* t);

    ConstraintImplementation *freeze();
    ConstraintImplementation *defrost();


    int arity() const ;
    int num_active() const ;
    int get_active(const int i) const ;
    Variable* get_scope() ;

    PropagationOutcome propagate();
    PropagationOutcome propagate(const Event evt);
    void restore();

    int get_backtrack_level();
    Decision get_decision();
    
    bool operator==(Constraint c) const {return propagator == c.propagator;}
    bool operator!=(Constraint c) const {return propagator != c.propagator;}

    std::ostream& display(std::ostream& os) const;

};



  typedef TwoWayStack< Triplet < int, Event, ConstraintImplementation*> > VariableQueue;

  class Environment {
  public:
    
    /*!@name Parameters*/
    //@{
    int level;

    Vector< int >                 saved_vars;
    Vector< Constraint >          saved_cons;
    Vector< int* >                saved_bools;
    Vector< ReversibleSet* >      saved_lists;
    Vector< ReversibleNum<int>* > saved_ints;

    /// The delimitation between different levels is kept by this vector of integers
    Vector< int > trail_;

    VariableQueue active_variables;

    ConstraintImplementation *taboo_constraint;
    //@}

    /*!@name Constructors*/
    //@{
    Environment() { 
      level = 0;
      taboo_constraint = NULL;
    }
    virtual ~Environment() {}
    //@}


    /*!@name Backtrack method*/
    //@{
    inline void save() {

      trail_.add(saved_vars.size);
      trail_.add(saved_cons.size);
      trail_.add(saved_bools.size);
      trail_.add(saved_lists.size);
      trail_.add(saved_ints.size);

      ++level;

    }


    void trigger_event(const int var, const Event evt) {
      Triplet< int, int, ConstraintImplementation* > t(var, evt, taboo_constraint);
      active_variables.push_back(t);
    }

    void _restore_();

    inline void save(ReversibleNum<int> *r) {saved_ints.add(r);}
    inline void save(ReversibleSet *r) {saved_lists.add(r);}
    inline void save(int *r) {saved_bools.add(r);}
    
    inline void save(int r) {saved_vars.add(r);}
    inline void save(Constraint r) {saved_cons.add(r);}
    //@}

  };

 


  /********************************************
   * Reversible Objects
   ********************************************/
  /*! \class Reversible
    \brief Backtrackable data structures.

    All structures that need to be restored to previous
    state during search implement the methods 
    save() and restore(). 
  */
  //class Environment;
  class Reversible {
  public:
    /*!@name Parameters*/
    //@{
    Environment *env;
    //@}

    /*!@name Constructors*/
    //@{
    Reversible() { env=NULL; }
    Reversible(Environment *s) : env(s) {;}
    void initialise(Environment *s) {env=s;}
    virtual ~Reversible() {}
    //@}

  };



  /********************************************
   * Reversible Boolean Domain
   ********************************************/
  /*! \class ReversibleBool
    \brief Backtrackable Primitive
  */
  class ReversibleBool : public Reversible
  {
  
  public:
    /*!@name Parameters*/
    //@{  
    int value;
    //@}

    /*!@name Constructors*/
    //@{ 
    ReversibleBool() : Reversible() {
      value = 3;
    }
    ReversibleBool(Environment *s) 
      : Reversible(s)
    {
      value = 3;
    }

    void initialise(Environment *s) {
      Reversible::initialise(s);
    }
    
    virtual ~ReversibleBool() {}
    //@}


    /*!@name Backtrack method*/
    //@{
    inline void save() { 
      if(value == 3) env->save(&value); 
    }
    inline void restore() { 
      value = 3;
    }
    //@}

    /*!@name Manipulation*/
    //@{  
    inline void remove(const int v) { 
      save();
      value = 2-v;
    }

    inline void operator=(const int v) { 
      save();
      value = 1+v;
    }

    inline bool contain(const int v) {
      return value&(v+1);
    }

    inline bool equal(const int v) {
      return (value&(v+1)) == value;
    }
    
    inline int size() { return (1+(value==3)); }
    //@}
    
    /*!@name Printing*/
    //@{
    std::ostream& display(std::ostream& os) const { os << (value == 3 ? "{0,1}" : (value == 1 ? "0" : "1")) ; return os; }
    //@}
  };

  /********************************************
   * Reversible Primitive Type
   ********************************************/
  /*! \class ReversibleNum
    \brief Backtrackable Primitive
  */
  template < class PRIMITIVE_TYPE >
  class ReversibleNum : public Reversible
  {
  
  public:
    /*!@name Parameters*/
    //@{  
    /// value trail
    Vector< PRIMITIVE_TYPE > trail_;
    /// current value
    PRIMITIVE_TYPE value;
    //@}

    /*!@name Constructors*/
    //@{ 
    ReversibleNum() : Reversible() {
    }
    ReversibleNum(const PRIMITIVE_TYPE v) 
    {
      initialise(v);
    }
    ReversibleNum(Environment *s) 
      : Reversible(s)
    {
      Reversible::initialise(s);
    }
    ReversibleNum(Environment *s, const PRIMITIVE_TYPE v) 
      : Reversible(s)
    {
      Reversible::initialise(s);
      initialise(v);
    }

    void initialise(Environment *s, const PRIMITIVE_TYPE v) 
    {
      Reversible::initialise(s);
      initialise(v);
    }

    void initialise(Environment *s) 
    {
      Reversible::initialise(s);
    }

    void initialise(const PRIMITIVE_TYPE v) 
    {
      value = v;
      trail_.initialise(0, 16);
      trail_.add(value);
      trail_.add(-1);
    }
    virtual ~ReversibleNum() {}
    //@}

    /*!@name Accessors*/
    //@{  
    inline operator const PRIMITIVE_TYPE() const { return value; }
    //@}  

    /*!@name Backtrack method*/
    //@{
    inline void save() { 
      if((int)trail_.back() != env->level) { 
	env->save(this); 
	trail_.add((int)value); 
	trail_.add(env->level); 
      } 
    }
    inline void restore() { 
      trail_.pop();//current_level); 
      value = (PRIMITIVE_TYPE)(trail_.pop()); 
    }
    //@}

    /*!@name Manipulation*/
    //@{  
    inline bool operator!() { return !value; }
    inline PRIMITIVE_TYPE operator-() { return -value; }
    inline void operator=  ( const PRIMITIVE_TYPE x ) { save(); value  = x; }
    inline void operator+= ( const PRIMITIVE_TYPE x ) { save(); value += x; }
    inline void operator-= ( const PRIMITIVE_TYPE x ) { save(); value -= x; }
    inline void operator*= ( const PRIMITIVE_TYPE x ) { save(); value *= x; }
    inline void operator/= ( const PRIMITIVE_TYPE x ) { save(); value /= x; }
    inline void operator|= ( const PRIMITIVE_TYPE x ) { save(); value |= x; }
    inline void operator&= ( const PRIMITIVE_TYPE x ) { save(); value &= x; }
    inline void operator^= ( const PRIMITIVE_TYPE x ) { save(); value ^= x; }
    inline ReversibleNum< PRIMITIVE_TYPE >& operator++ () { save(); ++value; return *this; }
    inline ReversibleNum< PRIMITIVE_TYPE >& operator-- () { save(); --value; return *this; } 
    //@}
  };


  /********************************************
   * Reversible Set
   ********************************************/
  /*! \class ReversibleSet
    \brief Backtrackable IntStack
  */
  class ReversibleSet : public Reversible, public IntStack {
    
  public:

    /*!@name Parameters*/
    //@{  
    /// value trail
    Vector< int > trail_;
    //@}

    /*!@name Constructors*/
    //@{ 
    ReversibleSet() : Reversible(), IntStack() {}
    ReversibleSet(Environment *s, const int lb=0, const int ub=0, const bool full=true)
      : Reversible(s)
    {
      initialise(lb, ub, full);
    }

    virtual ~ReversibleSet()
    {
    }

    virtual void initialise(Environment *s, const int lb, const int ub, const bool full)
    {
      Reversible::initialise(s);
      initialise(lb, ub, full);
    }

    virtual void initialise(const int lb, const int ub, const bool full)
    {
      int l = lb;
      int u = ub;
      if(l>u) {
	u = l-1;
	l = 0;
      }
      IntStack::initialise(l, u, full);
      trail_.initialise(0, 2*(u-l+1));
      trail_.add(-1);
    }
    //@}

    /*!@name Backtrack method*/
    //@{    
    inline void restore() { 
      trail_.pop(); size = trail_.pop(); 
    } 
    inline void save() { 
      if(trail_.back() != env->level) {
	trail_.add(size);
	trail_.add(env->level);
	env->save(this);
      }
    }
    //@}

    /*!@name Manipulation*/
    //@{  
    // it's either 'add's...
    inline void reversible_add(const int elt) {
      save();
      add(elt);
    }

    // ...or remove, but not both!!
    inline void reversible_remove(const int elt) {
      save();
      remove(elt);
    }

    inline int reversible_pop()
    {
      save();
      return pop();
    }

    inline int reversible_popHead()
    {
      save();
      return popHead();
    }
    //@}    

  };


  std::ostream& operator<< (std::ostream& os, const ReversibleBool& x);
  std::ostream& operator<< (std::ostream& os, const ReversibleBool* x);

}

#endif //__BACKTRACK_HPP
