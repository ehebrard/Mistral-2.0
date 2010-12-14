
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


//#include <mistral_solver.hpp>



#ifndef __BACKTRACK_HPP
#define __BACKTRACK_HPP

#include <mistral_global.hpp>


/*!
  A reversible structure implements:
  - a pointer to the store of reversible objects
  - a pointer to the current level
  - the level of last modification
  - a function save
  - a function restore

  Upon the first modification, it should:
  1/ put itself on the store for the current level
  
  When a new node is open, the save() function will be called on every stored object
  
  Upon backtrack over the current level, the restore() function will be called
 */

namespace Mistral {


  class Reversible;
  class Environment {
    public:
    
    int level;
    
    Vector< Reversible* > saved_objs;

    void save(Reversible *r);

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
    //bool first_change;
    //@}

    /*!@name Constructors*/
    //@{
    Reversible() {}
    Reversible(Environment *s) {env=s;}
    virtual ~Reversible() {}
    //@}

    /*!@name Modification method*/
    //@{
    //@}

    /*!@name Backtrack method*/
    //@{
    virtual void restore() = 0; 
    //virtual void save() = 0; 
    //@}

    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    virtual void debug_print() const {
      std::cout << "r" << std::endl;
    }
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
    ReversibleNum(PRIMITIVE_TYPE& v, Environment *s) : Reversible(s)
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

    /*!@name Backtrack*/
    //@{
    inline void save() { 
//       int lvl = env->level;
//       if(current_level != lvl) { 
// 	env->save(this); 
// 	trail_.add(value); 
// 	trail_.add(lvl); 
//       } 
//       if(current_level != env->level) { 
// 	env->save(this); 
// 	trail_.add(value); 
// 	trail_.add(env->level); 
//       } 

      if(trail_.back() != env->level) { 
	env->save(this); 
	trail_.add(value); 
	trail_.add(env->level); 
      } 
    }
    void restore() { 
      trail_.pop();//current_level); 
      trail_.pop(value); 
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
   * Reversible IntStack
   ********************************************/
  /*! \class ReversibleIntStack
    \brief Backtrackable IntStack
  */
  //typedef IntegerStack< ReversibleNum< unsigned int > > ReversibleIntStack;
  class ReversibleIntStack : public Reversible, public IntStack {
    
  public:
    Vector< unsigned int > trail_;
    Vector< int > lvl_;
    
    ReversibleIntStack(const int lb, const int ub, bool full=true)
    {
      initialise(lb, ub, full);
    }

    virtual ~ReversibleIntStack()
    {
    }

    virtual void initialise(const int lb, const int ub, const bool full=true)
    {
      IntStack::initialise(lb, ub, full);
      trail_.initialise(0, (ub-lb+1));
    }
    
    virtual void restore() { trail_.pop(size); } 
    virtual void save() {}

    void _save_() {

    }
 

  };


  /********************************************
   * Reversible MultiList
   ********************************************/
  /*! \class ReversibleMultiList
    \brief Backtrackable MultiList
  */
  template < class DATA_TYPE, int NUM_HEAD >
  class ReversibleMultiList : public MultiList< DATA_TYPE, NUM_HEAD >, public Reversible
  {
  
  public:
    /*!@name Parameters*/
    //@{  
    Vector< int > trail_;
    //@}

    /*!@name Constructors*/
    //@{
    ReversibleMultiList(Environment *env) 
      : MultiList< DATA_TYPE, NUM_HEAD >(), Reversible(env) {
      trail_.initialise(0, 4+4*NUM_HEAD);
      trail_.add(-1);
    }
    //@}

    /*!@name List Manipulation*/
    //@{
    inline void reversible_erase(const int idx, const int k=0) {
      _notify_();
      MultiList< DATA_TYPE, NUM_HEAD >::erase(idx, k);      
      _save_(idx, k);
    }
//     inline int reversible_add(const int idx, const int k=0) {
//       _notify_();
//       return MultiList< DATA_TYPE, NUM_HEAD >::create(idx, k);
//     }
    inline int reversible_add(DATA_TYPE elt, const int k=0) {
      _notify_();
      return MultiList< DATA_TYPE, NUM_HEAD >::create(elt, k);
    }
    //@}


    /*!@name Backtrack method*/
    //@{
    virtual void restore() {
      _restore_();
    }

    inline void _restore_() {
      //first_change = true;
      int k, idx_start, idx_end, succ;

      // pop the level flag
      trail_.pop();

      // restore the deletions
      for(k=0; k<NUM_HEAD; ++k) {
	// pop the first element of the circular list
	trail_.pop(idx_start);

	if(idx_start) {
	  idx_end   = MultiList< DATA_TYPE, NUM_HEAD >::data.stack_[idx_start].prev;
	  
	  // insert [idx_start -> idx_end] right at the start of the list
	  succ = MultiList< DATA_TYPE, NUM_HEAD >::data.stack_[MultiList< DATA_TYPE, NUM_HEAD >::head[k]].next;
	  MultiList< DATA_TYPE, NUM_HEAD >::data.stack_[succ].prev = idx_end;
	  MultiList< DATA_TYPE, NUM_HEAD >::data.stack_[MultiList< DATA_TYPE, NUM_HEAD >::head[k]].next = idx_start;
	  MultiList< DATA_TYPE, NUM_HEAD >::data.stack_[idx_end].next = succ;
	  MultiList< DATA_TYPE, NUM_HEAD >::data.stack_[idx_start].prev = MultiList< DATA_TYPE, NUM_HEAD >::head[k];
	}
      }

      // restore the degree
      MultiList< DATA_TYPE, NUM_HEAD >::degree = trail_.pop();

      // restore the insertions
      idx_start = MultiList< DATA_TYPE, NUM_HEAD >::data.size;
      trail_.pop(idx_end);
      while(idx_start > idx_end) {
	k = --MultiList< DATA_TYPE, NUM_HEAD >::data.size;
	MultiList< DATA_TYPE, NUM_HEAD >::erase(k);
	++MultiList< DATA_TYPE, NUM_HEAD >::degree;
	--idx_start;
      }
    } 

    virtual void save() {}

    inline void _notify_() {
      if(trail_.back() != env->level) {
	env->save(this);
	trail_.add(MultiList< DATA_TYPE, NUM_HEAD >::data.size);
	trail_.add(MultiList< DATA_TYPE, NUM_HEAD >::degree);
	for(int k=0; k<NUM_HEAD; ++k)
	  trail_.add(0);
	trail_.add(env->level);
      }
    }

    inline void _save_(const int idx, const int k) {
      int succ, pred;
      pred = trail_.back(k+2);
      if(pred) {
	// set idx as the successor of the head of the circular list	
	succ = MultiList< DATA_TYPE, NUM_HEAD >::data.stack_[pred].next;
	
	MultiList< DATA_TYPE, NUM_HEAD >::data.stack_[pred].next = idx;
	MultiList< DATA_TYPE, NUM_HEAD >::data.stack_[succ].prev = idx;
	MultiList< DATA_TYPE, NUM_HEAD >::data.stack_[idx].next = succ;
	MultiList< DATA_TYPE, NUM_HEAD >::data.stack_[idx].prev = pred;
      } else {
	// set up a circular list with only 'idx' in
	trail_.setBack(idx, k+2);
	MultiList< DATA_TYPE, NUM_HEAD >::data.stack_[idx].next = idx;
	MultiList< DATA_TYPE, NUM_HEAD >::data.stack_[idx].prev = idx;
      }
    }
    //@}


    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    void reversible_debug_print(std::ostream& o) const 
    {
      //for(unsigned int i=MultiList< DATA_TYPE, NUM_HEAD >::head[0]; i<MultiList< DATA_TYPE, NUM_HEAD >::data.size; ++i)
      //o << MultiList< DATA_TYPE, NUM_HEAD >::data.stack_[i];
      MultiList<int,NUM_HEAD>::debug_print(o);
      o << " / " << trail_
	<< " ("  << MultiList< DATA_TYPE, NUM_HEAD >::degree << ")";
    }
    //@}

  };

  
  class Constraint;
  class ConstraintTrigger {
  public:
    Constraint *constraint;
    int index;

    ConstraintTrigger(Constraint* con, const int id) {
      constraint = con;
      index = id;
    }

    ConstraintTrigger(const int c=0) {
      constraint = (Constraint*)c;
      index = -1;
    }

    bool operator!() { return !constraint; }
    ConstraintTrigger& operator=(Constraint *con) { 
      constraint = con; 
      return *this;
    }
    ConstraintTrigger& operator=(const int id) { 
      index = id; 
      return *this;
    }

//     std::string getString() const {
//       return (constraint ? constraint->getString() : ".");
//     }
    virtual std::ostream& display(std::ostream& os) const {
      if(constraint) 
	os << constraint;
      else
	os << ".";
      return os;
    }
  };

  typedef Node< ConstraintTrigger > ConstraintNode;
  typedef ReversibleMultiList< ConstraintTrigger, 3 > ConstraintList;

  std::ostream& operator<< (std::ostream& os, const ConstraintTrigger& x);

  std::ostream& operator<< (std::ostream& os, const ConstraintTrigger* x);



}

#endif //__BACKTRACK_HPP
