/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Guido Tack <tack@gecode.org>
 *
 *  Copyright:
 *     Guido Tack, 2007
 *
 *  Last modified:
 *     $Date: 2010-07-21 11:42:47 +0200 (Wed, 21 Jul 2010) $ by $Author: tack $
 *     $Revision: 11243 $
 *
 *  This file is part of Gecode, the generic constraint
 *  development environment:
 *     http://www.gecode.org
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the
 *  "Software"), to deal in the Software without restriction, including
 *  without limitation the rights to use, copy, modify, merge, publish,
 *  distribute, sublicense, and/or sell copies of the Software, and to
 *  permit persons to whom the Software is furnished to do so, subject to
 *  the following conditions:
 *
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#ifndef __GECODE_FLATZINC_HH__
#define __GECODE_FLATZINC_HH__

#include <iostream>
#include <string>
#include <map>

#include "conexpr.hpp"
#include "ast.hpp"
#include "varspec.hpp"

#include <mistral_solver.hpp>

using namespace Mistral;

typedef Vector<Variable> IntVarArray;
typedef Vector<Variable> BoolVarArray;
typedef Vector<Variable> SetVarArray;


/**
 * \namespace Gecode::FlatZinc
 * \brief Interpreter for the %FlatZinc language
 *
 * The Gecode::FlatZinc namespace contains all functionality required
 * to parse and solve constraint models written in the %FlatZinc language.
 *
 */

namespace FlatZinc {

  /**
   * \brief Output support class for %FlatZinc interpreter
   *
   */
  class Printer {
  private:
    AST::Array* _output;
    void printElem(std::ostream& out,
                   Solver& solver,
                   AST::Node* ai,
                   const IntVarArray& iv,
                   const BoolVarArray& bv,
                   const SetVarArray& sv
                   ) const;

  public:
    Printer(void) : _output(NULL) {}
    void init(AST::Array* output);

    void print(std::ostream& out,
               Solver& solver,
               const IntVarArray& iv,
               const BoolVarArray& bv,
               const SetVarArray& sv
               ) const;

    ~Printer(void);
  private:
    Printer(const Printer&);
    Printer& operator=(const Printer&);
  };

  /**
   * \brief A space that can be initialized with a %FlatZinc model
   *
   */
  class FlatZincModel {
  public:
    enum Meth {
      SATISFACTION, //< Solve as satisfaction problem
      MINIMIZATION, //< Solve as minimization problem
      MAXIMIZATION  //< Solve as maximization problem
    };
  protected:
    /// Mistral stuff
    Solver &solver;
    BranchingHeuristic *heuristic;
    RestartPolicy *policy;


    /// Number of integer variables
    int intVarCount;
    /// Number of Boolean variables
    int boolVarCount;
    /// Number of set variables
    int setVarCount;

    /// Index of the integer variable to optimize
    int _optVar;

    /// Whether to solve as satisfaction or optimization problem
    Meth _method;

    /// Annotations on the solve item
    AST::Array* _solveAnnotations;
  public:
    /// The integer variables
    IntVarArray iv;
    /// Indicates whether an integer variable is introduced by mzn2fzn
    std::vector<bool> iv_introduced;
    /// Indicates whether an integer variable aliases a Boolean variable
    std::vector<int> iv_boolalias;
    /// The Boolean variables
    BoolVarArray bv;
    /// Indicates whether a Boolean variable is introduced by mzn2fzn
    std::vector<bool> bv_introduced;
    /// The Set variables
    SetVarArray sv;
    /// Indicates whether a set variable is introduced by mzn2fzn
    std::vector<bool> sv_introduced;

    /// vars fixed to true and false, in case they are encountered often
    //Variable vartrue(1,1);
    //Variable varfalse(0,0);
    //std::map<int, Variable> constants;

    /// Construct empty space
    FlatZincModel(Solver& s);

    /// Destructor
    ~FlatZincModel(void);

    /// Initialize space with given number of variables
    void init(int intVars, int boolVars, int setVars);

    /// Create new integer variable from specification
    void newIntVar(IntVarSpec* vs);
    /// Link integer variable \a iv to Boolean variable \a bv
    void aliasBool2Int(int iv, int bv);
    /// Return linked Boolean variable for integer variable \a iv
    int aliasBool2Int(int iv);
    /// Create new Boolean variable from specification
    void newBoolVar(BoolVarSpec* vs);
    /// Create new set variable from specification
    void newSetVar(SetVarSpec* vs);

    /// Post a constraint specified by \a ce
    void postConstraint(const ConExpr& ce, AST::Node* annotation);

    /// Post the solve item
    void solve(AST::Array* annotation);
    /// Post that integer variable \a var should be minimized
    void minimize(int var, AST::Array* annotation);
    /// Post that integer variable \a var should be maximized
    void maximize(int var, AST::Array* annotation);

    /// setup parameters from the command line
    void set_parameters(SolverParameters& p);

    /// setup parameters from the command line
    void set_strategy(std::string var_o, std::string val_o, std::string r_pol);

    /// Run the search
    void run(std::ostream& out, const Printer& p);

    /// Produce output on \a out using \a p
    void print(std::ostream& out, const Printer& p) const;

    /**
     * \brief Remove all variables not needed for output
     *
     * After calling this function, no new constraints can be posted through
     * FlatZinc variable references, and the createBranchers method must
     * not be called again.
     *
     */
    void shrinkArrays(Printer& p);

    /// Return whether to solve a satisfaction or optimization problem
    Meth method(void) const;

    /// Return index of variable used for optimization
    int optVar(void) const;

    /**
     * \brief Create branchers corresponding to the solve item annotations
     *
     * If \a ignoreUnknown is true, unknown solve item annotations will be
     * ignored, otherwise a warning is written to \a err.
     */
    void createBranchers(AST::Node* ann, bool ignoreUnknown,
                         std::ostream& err = std::cerr);

    /// Return the solve item annotations
    AST::Array* solveAnnotations(void) const;

    /// Implement optimization
    //void constrain();

    /// options
    bool findall; // find all solutions
  };

  /// %Exception class for %FlatZinc errors
  class Error {
  private:
    const std::string msg;
  public:
    Error(const std::string& where, const std::string& what)
    : msg(where+": "+what) {}
    const std::string& toString(void) const { return msg; }
  };

  /**
   * \brief Parse FlatZinc file \a fileName into \a fzs and return it.
   *
   * Creates a new empty FlatZincSpace if \a fzs is NULL.
   */
  FlatZincModel* parse(const std::string& fileName,
                       Solver& solver,
                       Printer& p, std::ostream& err = std::cerr,
                       FlatZincModel* fzs=NULL);

  /**
   * \brief Parse FlatZinc from \a is into \a fzs and return it.
   *
   * Creates a new empty FlatZincSpace if \a fzs is NULL.
   */
  FlatZincModel* parse(std::istream& is,
                       Solver& solver,
                       Printer& p, std::ostream& err = std::cerr,
                       FlatZincModel* fzs=NULL);

}

#endif

// STATISTICS: flatzinc-any
