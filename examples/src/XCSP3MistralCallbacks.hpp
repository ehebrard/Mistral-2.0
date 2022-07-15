/*=============================================================================
* parser for CSP instances represented in XCSP3 Format
*
* Copyright (c) 2015 xcp3.org (contact <at> xcsp3.org)
* Copyright (c) 2008 Olivier ROUSSEL (olivier.roussel <at> cril.univ-artois.fr)
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
* THE SOFTWARE.
*=============================================================================
*/
#ifndef COSOCO_XCSP3PRINTCALLBACKS_H
#define COSOCO_XCSP3PRINTCALLBACKS_H

#include "XCSP3CoreCallbacks.h"
#include "XCSP3Variable.h"

#include "tree.hpp"

#include "mistral_solver.hpp"
#include <map>

/**
* This is an example that prints useful informations of a XCSP3 instance.
* You need to create your own class and to override functions of the callback.
* We suggest to make a map between XVariable and your own variables in order to
* facilitate the constructions of constraints.
*
* see main.cc to show declaration of the parser
*
*/

// #define _VERBOSE_ true
// #define _DEBUG_CUMULATIVE
// #define _DEBUG_MDD
// #define _DEBUG_REG true

// #define _REF_SOL_ true

#ifdef _VERBOSE_
#define _ID_(e) e
#else
#define _ID_(e) " "
#endif



//#define _OLD_TREE



namespace XCSP3Core {

class XCSP3MistralCallbacks : public XCSP3CoreCallbacks {
public:
  Mistral::Solver &solver;

  map<string, Mistral::VarArray> array;
  map<string, Mistral::Variable> variable;
  map<string, int> id_map;
  // BitSet var_set;

  Mistral::VarArray variables;
  std::vector<string> var_ids;
  std::vector<string> declared_var_ids;

  Mistral::BitSet *last_domain;
  Mistral::Variable last_var;
  Mistral::Vector<const int *> *last_table; // THIS WAS BUGGY BECAUSE LAST TABLE
                                            // WAS CREATED BY EXPENDING STARS
                                            // WRT THE VARIABLES' DOMAINS, WHICH
                                            // MAY BE DIFFERENT IN THE NEXT
                                            // CONSTRAINT
  vector<vector<int>> *last_tuple_list;

  Mistral::Goal *goal;

  Mistral::Vector<Mistral::Variable> _demand_;

  vector<int> util;

  vector<int> initial_degree;
// vector<int> initial_size;

#ifdef _REF_SOL_
  vector<int> reference;
#endif

  XCSP3MistralCallbacks(Mistral::Solver &s);

  virtual void beginInstance(InstanceType type) override;

  virtual void endInstance() override;

  virtual void beginVariables() override;

  virtual void endVariables() override;

  virtual void beginVariableArray(string id) override;

  virtual void endVariableArray() override;

  virtual void beginConstraints() override;

  virtual void endConstraints() override;

  virtual void beginGroup(string id) override;

  virtual void endGroup() override;

  virtual void beginBlock(string classes) override;

  virtual void endBlock() override;

  virtual void beginSlide(string id, bool circular) override;

  virtual void endSlide() override;

  virtual void beginObjectives() override;

  virtual void endObjectives() override;

  virtual void buildVariableInteger(string id, int minValue,
                                    int maxValue) override;

  virtual void buildVariableInteger(string id, vector<int> &values) override;

  virtual void buildConstraintExtension(string id, vector<XVariable *> list,
                                        vector<vector<int>> &tuples,
                                        bool support, bool hasStar) override;

  virtual void buildConstraintExtension(string id, XVariable *variable,
                                        vector<int> &tuples, bool support,
                                        bool hasStar) override;

  virtual void buildConstraintExtensionAs(string id, vector<XVariable *> list,
                                          bool support, bool hasStar) override;

  virtual void buildConstraintIntension(string id, string expr) override;

  virtual void buildConstraintIntension(string id, Tree *tree) override;

  virtual void buildConstraintPrimitive(string id, OrderType op, XVariable *x,
                                        int k, XVariable *y) override;

  virtual void buildConstraintPrimitive(string id, OrderType op, XVariable *x,
                                        int k) override;

  virtual void buildConstraintPrimitive(string id, XVariable *x, bool in,
                                        int min, int max) override;

  virtual void buildConstraintPrecedence(string id, vector<XVariable *> &list,
                                         vector<int> values) override;

  virtual void buildConstraintClause(string id, vector<XVariable *> &positive,
                                     vector<XVariable *> &negative) override;

  virtual void
  buildConstraintRegular(string id, vector<XVariable *> &list, string st,
                         vector<string> & final,
                         vector<XTransition> &transitions) override;

  virtual void buildConstraintMDD(string id, vector<XVariable *> &list,
                                  vector<XTransition> &transitions) override;

  virtual void buildConstraintAlldifferent(string id,
                                           vector<XVariable *> &list) override;

  virtual void buildConstraintAlldifferentExcept(string id,
                                                 vector<XVariable *> &list,
                                                 vector<int> &except) override;

  virtual void
  buildConstraintAlldifferentList(string id,
                                  vector<vector<XVariable *>> &lists) override;

  virtual void buildConstraintAlldifferentMatrix(
      string id, vector<vector<XVariable *>> &matrix) override;

  virtual void buildConstraintAllEqual(string id,
                                       vector<XVariable *> &list) override;

  virtual void buildConstraintNotAllEqual(string id,
                                          vector<XVariable *> &list) override;

  virtual void buildConstraintOrdered(string id, vector<XVariable *> &list,
                                      OrderType order) override;

  virtual void buildConstraintLex(string id, vector<vector<XVariable *>> &lists,
                                  OrderType order) override;

  virtual void buildConstraintLexMatrix(string id,
                                        vector<vector<XVariable *>> &matrix,
                                        OrderType order) override;

  virtual void buildConstraintSum(string id, vector<XVariable *> &list,
                                  vector<int> &coeffs,
                                  XCondition &cond) override;

  virtual void buildConstraintSum(string id, vector<XVariable *> &list,
                                  XCondition &cond) override;

  virtual void buildConstraintSum(string id, vector<XVariable *> &list,
                                  vector<XVariable *> &coeffs,
                                  XCondition &cond) override;

  virtual void buildConstraintSum(string id, vector<Tree *> &trees,
                                  XCondition &cond) override;

  virtual void buildConstraintSum(string id, vector<Tree *> &trees,
                                  vector<int> &coefs,
                                  XCondition &cond) override;

  virtual void buildConstraintAtMost(string id, vector<XVariable *> &list,
                                     int value, int k) override;

  virtual void buildConstraintAtLeast(string id, vector<XVariable *> &list,
                                      int value, int k) override;

  virtual void buildConstraintExactlyK(string id, vector<XVariable *> &list,
                                       int value, int k) override;

  virtual void buildConstraintAmong(string id, vector<XVariable *> &list,
                                    vector<int> &values, int k) override;

  virtual void buildConstraintExactlyVariable(string id,
                                              vector<XVariable *> &list,
                                              int value, XVariable *x) override;

  virtual void buildConstraintCount(string id, vector<XVariable *> &list,
                                    vector<int> &values,
                                    XCondition &xc) override;

  virtual void buildConstraintCount(string id, vector<XVariable *> &list,
                                    vector<XVariable *> &values,
                                    XCondition &xc) override;

  virtual void buildConstraintCount(string id, vector<Tree *> &trees,
                                    vector<int> &values,
                                    XCondition &xc) override;

  virtual void buildConstraintCount(string id, vector<Tree *> &trees,
                                    vector<XVariable *> &values,
                                    XCondition &xc) override;

  virtual void buildConstraintNValues(string id, vector<XVariable *> &list,
                                      vector<int> &except,
                                      XCondition &xc) override;

  virtual void buildConstraintNValues(string id, vector<XVariable *> &list,
                                      XCondition &xc) override;

  virtual void buildConstraintCardinality(string id, vector<XVariable *> &list,
                                          vector<int> values,
                                          vector<int> &occurs,
                                          bool closed) override;

  virtual void buildConstraintCardinality(string id, vector<XVariable *> &list,
                                          vector<int> values,
                                          vector<XVariable *> &occurs,
                                          bool closed) override;

  virtual void buildConstraintCardinality(string id, vector<XVariable *> &list,
                                          vector<int> values,
                                          vector<XInterval> &occurs,
                                          bool closed) override;

  virtual void buildConstraintCardinality(string id, vector<XVariable *> &list,
                                          vector<XVariable *> values,
                                          vector<int> &occurs,
                                          bool closed) override;

  virtual void buildConstraintCardinality(string id, vector<XVariable *> &list,
                                          vector<XVariable *> values,
                                          vector<XVariable *> &occurs,
                                          bool closed) override;

  virtual void buildConstraintCardinality(string id, vector<XVariable *> &list,
                                          vector<XVariable *> values,
                                          vector<XInterval> &occurs,
                                          bool closed) override;

  virtual void buildConstraintMinimum(string id, vector<XVariable *> &list,
                                      XCondition &xc) override;

  virtual void buildConstraintMinimum(string id, vector<XVariable *> &list,
                                      XVariable *index, int startIndex,
                                      RankType rank, XCondition &xc) override;

  virtual void buildConstraintMinimum(string id, vector<Tree *> &list,
                                      XCondition &xc) override;

  virtual void buildConstraintMaximum(string id, vector<XVariable *> &list,
                                      XCondition &xc) override;

  virtual void buildConstraintMaximum(string id, vector<XVariable *> &list,
                                      XVariable *index, int startIndex,
                                      RankType rank, XCondition &xc) override;

  virtual void buildConstraintMaximum(string id, vector<Tree *> &list,
                                      XCondition &xc) override;

  virtual void buildConstraintElement(string id, vector<XVariable *> &list,
                                      int value) override;

  virtual void buildConstraintElement(string id, vector<XVariable *> &list,
                                      XVariable *value) override;

  virtual void buildConstraintElement(string id, vector<XVariable *> &list,
                                      int startIndex, XVariable *index,
                                      RankType rank, int value) override;

  virtual void buildConstraintElement(string id, vector<XVariable *> &list,
                                      int startIndex, XVariable *index,
                                      RankType rank, XVariable *value) override;

  virtual void buildConstraintElement(string id, vector<int> &list,
                                      int startIndex, XVariable *index,
                                      RankType rank, XVariable *value) override;

  virtual void buildConstraintElement(string id, vector<XVariable *> &list,
                                      XVariable *index, int startIndex,
                                      XCondition &xc) override;

  virtual void buildConstraintElement(string id,
                                      vector<vector<XVariable *>> &matrix,
                                      int startRowIndex, XVariable *rowIndex,
                                      int startColIndex, XVariable *colIndex,
                                      XVariable *value) override;

  virtual void buildConstraintElement(string id,
                                      vector<vector<XVariable *>> &matrix,
                                      int startRowIndex, XVariable *rowIndex,
                                      int startColIndex, XVariable *colIndex,
                                      int value) override;

  virtual void buildConstraintElement(string id, vector<vector<int>> &matrix,
                                      int startRowIndex, XVariable *rowIndex,
                                      int startColIndex, XVariable *colIndex,
                                      XVariable *value) override;

  virtual void buildConstraintChannel(string id, vector<XVariable *> &list,
                                      int startIndex) override;

  virtual void buildConstraintChannel(string id, vector<XVariable *> &list1,
                                      int startIndex1,
                                      vector<XVariable *> &list2,
                                      int startIndex2) override;

  virtual void buildConstraintChannel(string id, vector<XVariable *> &list,
                                      int startIndex,
                                      XVariable *value) override;

  virtual void buildConstraintStretch(string id, vector<XVariable *> &list,
                                      vector<int> &values,
                                      vector<XInterval> &widths) override;

  virtual void buildConstraintStretch(string id, vector<XVariable *> &list,
                                      vector<int> &values,
                                      vector<XInterval> &widths,
                                      vector<vector<int>> &patterns) override;

  virtual void buildConstraintNoOverlap(string id, vector<XVariable *> &origins,
                                        vector<int> &lengths,
                                        bool zeroIgnored) override;

  virtual void buildConstraintNoOverlap(string id, vector<XVariable *> &origins,
                                        vector<XVariable *> &lengths,
                                        bool zeroIgnored) override;

  virtual void buildConstraintNoOverlap(string id,
                                        vector<vector<XVariable *>> &origins,
                                        vector<vector<int>> &lengths,
                                        bool zeroIgnored) override;

  virtual void buildConstraintNoOverlap(string id,
                                        vector<vector<XVariable *>> &origins,
                                        vector<vector<XVariable *>> &lengths,
                                        bool zeroIgnored) override;

  virtual void buildConstraintInstantiation(string id,
                                            vector<XVariable *> &list,
                                            vector<int> &values) override;

  void p_encode_cumulative(Mistral::Solver &s,
                           Mistral::Vector<Mistral::Variable> &start,
                           Mistral::Vector<Mistral::Variable> &dur,
                           Mistral::Vector<Mistral::Variable> &req,
                           Mistral::Variable cap);
  void p_cumulative_flow(Mistral::Solver &s,
                         Mistral::Vector<Mistral::Variable> &start,
                         Mistral::Vector<Mistral::Variable> &dur,
                         Mistral::Vector<Mistral::Variable> &req,
                         Mistral::Variable cap, int horizon);
  void p_cumulative_discretization(Mistral::Solver &s,
                                   Mistral::Vector<Mistral::Variable> &start,
                                   Mistral::Vector<Mistral::Variable> &dur,
                                   Mistral::Vector<Mistral::Variable> &req,
                                   Mistral::Variable cap, int horizon);

  virtual void buildConstraintCumulative(string id,
                                         vector<XVariable *> &origins,
                                         vector<int> &lengths,
                                         vector<int> &heights,
                                         XCondition &xc) override;

  virtual void buildConstraintCumulative(string id,
                                         vector<XVariable *> &origins,
                                         vector<int> &lengths,
                                         vector<XVariable *> &varHeights,
                                         XCondition &xc) override;

  virtual void buildConstraintCumulative(string id,
                                         vector<XVariable *> &origins,
                                         vector<XVariable *> &lengths,
                                         vector<int> &heights,
                                         XCondition &xc) override;

  virtual void buildConstraintCumulative(string id,
                                         vector<XVariable *> &origins,
                                         vector<XVariable *> &lengths,
                                         vector<XVariable *> &heights,
                                         XCondition &xc) override;

  virtual void buildConstraintCircuit(string id, vector<XVariable *> &list,
                                      int startIndex) override;

  virtual void buildConstraintCircuit(string id, vector<XVariable *> &list,
                                      int startIndex, int size) override;

  virtual void buildConstraintCircuit(string id, vector<XVariable *> &list,
                                      int startIndex, XVariable *size) override;

  virtual void buildObjectiveMinimizeExpression(string expr) override;

  virtual void buildObjectiveMaximizeExpression(string expr) override;

  virtual void buildObjectiveMinimizeVariable(XVariable *x) override;

  virtual void buildObjectiveMaximizeVariable(XVariable *x) override;

  virtual void buildObjectiveMinimize(ExpressionObjective type,
                                      vector<Tree *> &trees,
                                      vector<int> &coefs) override;

  virtual void buildObjectiveMaximize(ExpressionObjective type,
                                      vector<Tree *> &trees,
                                      vector<int> &coefs) override;

  virtual void buildObjectiveMinimize(ExpressionObjective type,
                                      vector<XVariable *> &list,
                                      vector<int> &coefs) override;

  virtual void buildObjectiveMaximize(ExpressionObjective type,
                                      vector<XVariable *> &list,
                                      vector<int> &coefs) override;

  virtual void buildObjectiveMinimize(ExpressionObjective type,
                                      vector<XVariable *> &list) override;

  virtual void buildObjectiveMaximize(ExpressionObjective type,
                                      vector<XVariable *> &list) override;

  void getVariables(vector<XVariable *> &list,
                    Mistral::Vector<Mistral::Variable> &scope);

  Mistral::Variable postExpression(XCSP3Mistral::Node *n, int level = 1);

  Mistral::Variable postExpression(Node *n, int level);
};
}

using namespace XCSP3Core;
using namespace Mistral;

XCSP3MistralCallbacks::XCSP3MistralCallbacks(Solver &s)
    : XCSP3CoreCallbacks(), solver(s) {
  goal = NULL;

  intensionUsingString = true;
  recognizeSpecialIntensionCases = false;
}

void XCSP3MistralCallbacks::getVariables(vector<XVariable *> &list,
                                         Vector<Variable> &scope) {
  for (size_t i = 0; i < list.size();
       i++) { // not necessary to check if the variable exists
    // if(variable.count( list[i]->id )) {
    scope.add(variable[list[i]->id]);
    ++initial_degree[id_map[list[i]->id]];
    // } else if(list[i]->domain == NULL) {
    // 	// std::cout << list[i]->id << " has an empty domain " <<
    // std::endl;
    // 	int v = stoi(list[i]->id);
    // 	scope.add(Variable(v,v));
    // } else if(list[i]->domain->values.size()==1) {
    // 	int v = list[i]->domain->minimum();
    // 	scope.add(Variable(v,v));
    // } else {
    // 	cout << "s UNSUPPORTED\n";
    // 	exit(1);
    // }
  }
}

template <class T> void displayList(vector<T> &list, string separator = " ") {
  if (list.size() > 8) {
    for (int i = 0; i < 3; i++)
      cout << list[i] << separator;
    cout << " ... ";
    for (size_t i = list.size() - 4; i < list.size(); i++)
      cout << list[i] << separator;
    cout << endl;
    return;
  }
  for (size_t i = 0; i < list.size(); i++)
    cout << list[i] << separator;
  cout << endl;
}

void displayList(vector<XVariable *> &list, string separator = " ") {
  if (list.size() > 8) {
    for (int i = 0; i < 3; i++)
      cout << list[i]->id << separator;
    cout << " ... ";
    for (unsigned int i = list.size() - 4; i < list.size(); i++)
      cout << list[i]->id << separator;
    cout << endl;
    return;
  }
  for (unsigned int i = 0; i < list.size(); i++)
    cout << list[i]->id << separator;
  cout << endl;
}

void XCSP3MistralCallbacks::beginInstance(InstanceType type) {
#ifdef _VERBOSE_
    cout << "Start Instance - type=" << type << endl;
#endif

    last_domain = NULL;
}


void XCSP3MistralCallbacks::endInstance() {
#ifdef _VERBOSE_
  cout << "End SAX parsing " << endl;
#endif

#ifdef _REF_SOL_
  ifstream valfile("ref_values", ios_base::in);
  ifstream varfile("ref_variables", ios_base::in);
  int v;
  string x;
  reference.resize(declared_var_ids.size());

  std::cout << "VARSIZE = " << solver.variables.size << std::endl;

  solver.reference_solution.resize(solver.variables.size, NOVAL);
#endif

  int i = 0;
  for (auto id : declared_var_ids) {
    // id_map[id] = variables.size;
    Variable X = variable[id];
    variables.add(X);
    var_ids.push_back(id);
    // initial_degree.push_back(X.get_degree());
    initial_degree[i++] += X.get_degree();
  }

// var_set.initialise(0, variables.size, BitSet::empt);

#ifdef _REF_SOL_
  // i=0;
  for (auto id : declared_var_ids) {
    valfile >> v;
    varfile >> x;
    // std::cout << declared_var_ids[id_map[id]] << " <- " << v << std::endl;
    reference[id_map[x]] = v;
    std::cout << x << " <- " << reference[id_map[x]] << "("
              << variables[id_map[x]] << " in "
              << variables[id_map[x]].get_domain() << ")" << std::endl;
    solver.reference_solution[variables[id_map[x]].id()] = v;
    // ++i;
  }
#endif
}


void XCSP3MistralCallbacks::beginVariables() {
#ifdef _VERBOSE_
  cout << " start variables declaration" << endl;
#endif
}


void XCSP3MistralCallbacks::endVariables() {
#ifdef _VERBOSE_
  cout << " end variables declaration" << endl << endl;
#endif

  int count = 0;
  for (auto id : declared_var_ids) {
    id_map[id] = count++;
    initial_degree.push_back(0);
  }
}


void XCSP3MistralCallbacks::beginVariableArray(string id) {
#ifdef _VERBOSE_
  cout << "    array: " << id << endl;
#endif
  //
  // VarArray A;
  // v_array.push_back(A);
}

void XCSP3MistralCallbacks::buildVariableInteger(string id, int minValue, int maxValue) {
#ifdef _VERBOSE_
  cout << "    var b " << id << " : " << minValue << "..." << maxValue << endl;
#endif

  Variable X(minValue, maxValue);
  variable[id] = X;
  // id_map[id] = declared_var_ids.size();
  // variables.add(X);
  declared_var_ids.push_back(id);
  // initial_size.push_back(maxValue-minValue+1);

  solver.add(X);
}


void XCSP3MistralCallbacks::buildVariableInteger(string id, vector<int> &values) {
#ifdef _VERBOSE_
  cout << "    var l " << id << " : ";
  cout << "        ";
  displayList(values);
#endif

  Variable X(values);
  variable[id] = X;
  // variables.add(X);
  declared_var_ids.push_back(id);

  solver.add(X);
}



void XCSP3MistralCallbacks::endVariableArray() {
}


void XCSP3MistralCallbacks::beginConstraints() {
#ifdef _VERBOSE_
  cout << " start constraints declaration" << endl;
#endif
}


void XCSP3MistralCallbacks::endConstraints() {
#ifdef _VERBOSE_
  cout << "\n end constraints declaration" << endl << endl;
#endif
}


void XCSP3MistralCallbacks::beginGroup(string id) {
#ifdef _VERBOSE_
  cout << "   start group of constraint " << id << endl;
#endif
}


void XCSP3MistralCallbacks::endGroup() {
#ifdef _VERBOSE_
  cout << "   end group of constraint" << endl;
#endif
}


void XCSP3MistralCallbacks::beginBlock(string classes) {
#ifdef _VERBOSE_
  cout << "   start block of constraint classes = " << classes << endl;
#endif
}


void XCSP3MistralCallbacks::endBlock() {
#ifdef _VERBOSE_
  cout << "   end group of constraint" << endl;
#endif
}


void XCSP3MistralCallbacks::beginSlide(string id, bool circular) {
#ifdef _VERBOSE_
  cout << "   start slide " << id << endl;
#endif
}


void XCSP3MistralCallbacks::endSlide() {
#ifdef _VERBOSE_
  cout << "   end slide" << endl;
#endif
}


void XCSP3MistralCallbacks::beginObjectives() {
#ifdef _VERBOSE_
  cout << "   start Objective " << endl;
#endif
}


void XCSP3MistralCallbacks::endObjectives() {
#ifdef _VERBOSE_
  cout << "   end Objective " << endl;
#endif
}


void XCSP3MistralCallbacks::buildConstraintExtension(string id, vector<XVariable *> list, vector <vector<int>> &tuples, bool support, bool hasStar) {

#ifdef _VERBOSE_
  cout << "\n    extension constraint : " << id << endl;
  cout << "        " << (support ? "support" : "conflict")
       << " arity:" << list.size() << " nb tuples: " << tuples.size()
       << " star: " << hasStar << endl;
  cout << "        ";
  displayList(list);
#endif

  VarArray scope;
  getVariables(list, scope);

  TableExpression *tab = new TableExpression(scope, support);

  // for ( auto x : scope ) {
  // 	cout << x.get_domain() << endl;
  // }
  // cout << endl;

  vector<int> stared;
  vector<int> fact;
  fact.resize(list.size());
  for (auto pt = begin(tuples); pt != end(tuples); ++pt) {
    // for ( auto& t : tuples ) {

    for (size_t i = 0; i < pt->size(); ++i) {
      // cout << " " << (*(pt))[i];
      if ((*(pt))[i] == 2147483647) {
        stared.push_back(i);
        fact[i] = scope[i].get_min();
      } else {
        fact[i] = (*(pt))[i];
      }
    }
    // cout << endl;

    // int count = 10000;
    if (stared.size() > 0) {

      int numtuples = 1;
      for (auto j : stared) {
        numtuples *= scope[j].get_size();
        if (numtuples > 1000000000) {
          cout << "s UNSUPPORTED" << _ID_(": too many tuples") << "\n";
          exit(1);
        }
      }

      int i = stared.size() - 1;
      bool cont = true;
      do {

        tab->add(&(fact[0]));
        numtuples--;

        int j = stared[i];
        int n = scope[j].next(fact[j]);

        // cout << " ++t[" << j << "]: " << " " << fact[j] << " -> " << n << "
        // (" << i << ")" << endl;

        while (n == fact[j]) {
          if (i == 0) {
            cont = false;
            break;
          }
          for (size_t k = i; k < stared.size(); ++k) {
            fact[stared[k]] = scope[stared[k]].get_min();
          }
          j = stared[--i];
          n = scope[j].next(fact[j]);
        }

        i = stared.size() - 1;

        fact[j] = n;

      } while (cont);

      stared.clear();

      assert(numtuples == 0);

      // std::cout << "#tuples = " << numtuples << std::endl;

    } else {
      tab->add(&((*(pt))[0]));
    }
  }

  last_table = tab->tuples;
  last_tuple_list = &tuples;

  solver.add(Variable(tab));
}

void XCSP3MistralCallbacks::buildConstraintExtension(string id, XVariable *var,
                                                     vector<int> &tuples,
                                                     bool support,
                                                     bool hasStar) {
#ifdef _VERBOSE_
  cout << "\n    extension constraint with one variable: " << id << endl;
  cout << "        " << (support ? "support" : "conflict")
       << " nb tuples: " << tuples.size() << " star: " << hasStar << endl;
  cout << (*var) << endl;
#endif

  ++initial_degree[id_map[var->id]];

  // cout << "here " << var->id << " / " << variable.size() << "\n";

  if (support) {
    auto the_min = *std::min_element(begin(tuples), end(tuples));
    auto the_max = *std::max_element(begin(tuples), end(tuples));

    if (last_domain)
      delete last_domain;
    last_domain = new BitSet(the_min, the_max, BitSet::empt);

    for (auto t : tuples) {
      last_domain->add(t);
    }

    // cout << *last_domain << endl;

    solver.add(Member(variable[var->id], *last_domain));

    // cout << 22 << endl;

  } else {

    Variable X = variable[var->id];

    if (last_domain)
      delete last_domain;
    last_domain = new BitSet(X.get_min(), X.get_max(), BitSet::full);

    for (auto t : tuples) {
      last_domain->remove(t);
    }

    solver.add(Member(X, *last_domain));
  }
}

void XCSP3MistralCallbacks::buildConstraintExtensionAs(string id,
                                                       vector<XVariable *> list,
                                                       bool support,
                                                       bool hasStar) {
#ifdef _VERBOSE_
  cout << "\n    extension constraint similar as previous one on ";
  displayList(list);
  cout << ": " << id << endl;
#endif

  if (hasStar) {
    buildConstraintExtension(id, list, *last_tuple_list, support, hasStar);
  } else {

    VarArray scope;
    getVariables(list, scope);

    if (scope.size > 1) {
      TableExpression *tab = new TableExpression(scope, support);
      for (auto t : *last_table) {
        tab->add(t);
      }

      solver.add(Variable(tab));
    } else {

      ++initial_degree[id_map[list[0]->id]];
      solver.add(Member(scope[0], *last_domain));
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintIntension(string id, string expr) {
#ifdef _VERBOSE_
  cout << "\n    intension constraint : " << id << " : " << expr << endl;
#endif

#ifdef _OLD_TREE
  XCSP3Mistral::Tree tree(expr);
  solver.add(postExpression(tree.root, true));
  tree.dispose();
#else
  Tree tree(expr);
  solver.add(postExpression(tree.root, true));
// delete tree.root;
#endif
}

void XCSP3MistralCallbacks::buildConstraintIntension(string id, Tree *tree) {
#ifdef _VERBOSE_
  cout << "\n    intension constraint with tree : " << id << endl;
#endif

  solver.add(postExpression(tree->root, true));
}

void XCSP3MistralCallbacks::buildConstraintPrimitive(string id, OrderType op,
                                                     XVariable *x, int k,
                                                     XVariable *y) {
#ifdef _VERBOSE_
  cout << "\n   intension constraint " << id << ": " << x->id
       << (k >= 0 ? "+" : "") << k << " op " << y->id << endl;
#endif

  Variable X = variable[x->id];
  Variable Y = variable[y->id];

  ++initial_degree[id_map[x->id]];
  ++initial_degree[id_map[y->id]];

  // LT, GE, GT, IN, EQ, NE)
  if (op == LE) {
    if (k == 0) {
      solver.add(X <= Y);
    } else if (k == 1) {
      solver.add(X < Y);
    } else {
      solver.add((X + k) <= Y);
    }
  } else if (op == LT) {
    if (k == 0) {
      solver.add(X < Y);
    } else if (k == -1) {
      solver.add(X <= Y);
    } else {
      solver.add((X + k) < Y);
    }
  } else if (op == GE) {
    if (k == 0) {
      solver.add(X >= Y);
    } else if (k == -1) {
      solver.add(X > Y);
    } else {
      solver.add((X + k) >= Y);
    }
  } else if (op == GT) {
    if (k == 0) {
      solver.add(X > Y);
    } else if (k == 1) {
      solver.add(X >= Y);
    } else {
      solver.add((X + k) > Y);
    }
  } else if (op == IN) {
    cout << "s UNSUPPORTED" << _ID_(": IN non-interval") << "\n";
    exit(1);
  } else if (op == EQ) {
    if (k == 0)
      solver.add(X == Y);
    else
      solver.add((X + k) == Y);
  } else if (op == NE) {
    if (k == 0)
      solver.add(X != Y);
    else
      solver.add((X + k) != Y);
  }
}

void XCSP3MistralCallbacks::buildConstraintPrimitive(string id, OrderType op,
                                                     XVariable *x, int k) {

#ifdef _VERBOSE_
  cout << "\n  intension constraint " << id << ": " << x->id << " op " << k
       << endl;
#endif

  Variable X = variable[x->id];

  ++initial_degree[id_map[x->id]];

  // LT, GE, GT, IN, EQ, NE)
  if (op == LE) {
    solver.add(X <= k);
  } else if (op == LT) {
    solver.add(X < k);
  } else if (op == GE) {
    solver.add(X >= k);
  } else if (op == GT) {
    solver.add(X > k);
  } else if (op == IN) {
    cout << "s UNSUPPORTED" << _ID_(": IN non-interval") << "\n";
    exit(1);
  } else if (op == EQ) {
    solver.add(X == k);
  } else if (op == NE) {
    solver.add(X != k);
  }
}

/**
 * If  #recognizeSpecialIntensionCases is enabled (this is the case by default)
 * intensional constraint of the form : x in/notin [min max] are recognized
 * If such a intensional constraint is recognized, a callback to this function
 * is done and not to  #buildConstraintIntension
 *
 * @param id the id (name) of the constraint
 * @param x the variable
 * @param in true if x is in this interval
 * @param min the constant
 * @param max the constant
 *
 */
void XCSP3MistralCallbacks::buildConstraintPrimitive(string id, XVariable *x,
                                                     bool in, int min,
                                                     int max) {

#ifdef _VERBOSE_
  cout << "\n  intension constraint " << id << ": " << x->id
       << (in ? " in " : " not in ") << "[" << min << ", " << max << "]\n";
#endif

  Variable X = variable[x->id];

  ++initial_degree[id_map[x->id]];

  if (in) {

    solver.add(X <= max);
    solver.add(X >= min);

  } else {

    solver.add((X < min) || (X > max));
  }
}

void XCSP3MistralCallbacks::buildConstraintPrecedence(string id,
                                                      vector<XVariable *> &list,
                                                      vector<int> values) {
  cout << "Prec: TODO!\n";
}

void XCSP3MistralCallbacks::buildConstraintClause(
    string id, vector<XVariable *> &positive, vector<XVariable *> &negative) {
  cout << "Clause: TODO!\n";
}

void XCSP3MistralCallbacks::buildConstraintRegular(
    string id, vector<XVariable *> &list, string start, vector<string> & final,
    vector<XTransition> &transitions) {
#ifdef _VERBOSE_
  cout << "\n    regular constraint " << id << endl;
  cout << "        ";
  displayList(list);
  cout << "        start: " << start << endl;
  cout << "        final: ";
  displayList(final, ",");
  cout << endl;
  cout << "        transitions: ";
  for (unsigned int i = 0;
       i < (transitions.size() > 4 ? 4 : transitions.size()); i++) {
    cout << "(" << transitions[i].from << "," << transitions[i].val << ","
         << transitions[i].to << ") ";
  }
  if (transitions.size() > 4)
    cout << "...";
  cout << endl;
#endif

  map<string, int> state_map;
  map<string, bool> final_map;

  for (auto t : transitions) {
    if (state_map.find(t.from) == end(state_map)) {
      state_map[t.from] = state_map.size();
      final_map[t.from] = false;
    }
    if (state_map.find(t.to) == end(state_map)) {
      state_map[t.to] = state_map.size();
      final_map[t.to] = false;
    }
  }

  for (auto f : final) {
    final_map[f] = true;
  }

  VarArray scope;
  getVariables(list, scope);

  VarArray state_var(scope.size - 1, 0, state_map.size() - 1);

#ifdef _DEBUG_REG
  std::cout << " states: " << state_var.size << ": "
            << state_var[0].get_domain() << std::endl;
#endif

  Vector<const int *> *trans = new Vector<const int *>;
  Vector<const int *> *first = new Vector<const int *>;
  Vector<const int *> *last = new Vector<const int *>;

  for (const auto &t : transitions) {
    int *tuple = new int[3];
    tuple[0] = state_map[t.from];
    tuple[1] = t.val;
    tuple[2] = state_map[t.to];

#ifdef _DEBUG_REG
    std::cout << "  T " << tuple[0] << " " << tuple[1] << " " << tuple[2]
              << std::endl;
#endif

    trans->add(tuple);
    if (start == t.from) {
      int *pair = new int[2];
      pair[0] = tuple[1];
      pair[1] = tuple[2];
      first->add(pair);

#ifdef _DEBUG_REG
      std::cout << "  F " << tuple[0] << " " << tuple[1] << std::endl;
#endif
    }
    if (final_map[t.to]) {
      int *pair = new int[2];
      pair[0] = tuple[0];
      pair[1] = tuple[1];
      last->add(pair);

#ifdef _DEBUG_REG
      std::cout << "  L " << tuple[0] << " " << tuple[1] << std::endl;
#endif
    }
  }

  VarArray S;
  S.add(scope[0]);
  S.add(state_var[0]);
  solver.add(Table(S, first));

  S.clear();
  S.add(state_var.back());
  S.add(scope.back());
  solver.add(Table(S, last));

  for (size_t i = 1; i < state_var.size; ++i) {
    S.clear();
    S.add(state_var[i - 1]);
    S.add(scope[i]);
    S.add(state_var[i]);
    solver.add(Table(S, trans));
  }

  // exit(1);
}

void XCSP3MistralCallbacks::buildConstraintMDD(
    string id, vector<XVariable *> &list, vector<XTransition> &transitions) {
#ifdef _VERBOSE_
  cout << "\n    mdd constraint" << endl;
  cout << "        ";
  displayList(list);
  cout << "        transitions: ";
  for (unsigned int i = 0;
       i < (transitions.size() > 4 ? 4 : transitions.size()); i++) {
    cout << "(" << transitions[i].from << "," << transitions[i].val << ","
         << transitions[i].to << ") ";
  }
  if (transitions.size() > 4)
    cout << "...";
  cout << endl;

// for(unsigned int i = 0; i < transitions.size(); i++) {
//     cout << "(" << transitions[i].from << "," << transitions[i].val << "," <<
//     transitions[i].to << ") ";
// }
// cout << endl;
#endif

  // for( const auto t : transitions ) {
  // 	std::cout << t.from << "-" << t.val << "-" << t.to << std::endl;
  // }

  // for( unsigned int i = 0; i < transitions.size(); ++i ) {
  // 	std::cout << transitions[i].from ;
  // 	std::cout << "-" << transitions[i].val ;
  // 	std::cout << "-" << transitions[i].to << std::endl;
  // }
  //
  // exit(1);

  VarArray scope;
  getVariables(list, scope);

  string root = begin(transitions)->from;
  // string final = transitions.rbegin()->to;
  map<string, bool> is_final;

  map<string, int> layer_map;
  map<string, int> node_map;
  map<string, string> prev_map;
  auto num_layers{0};

  // get layer of every node
  for (const auto &t : transitions) {
    auto l{0};
    if (t.from != root) {
      l = layer_map[t.from] + 1;
      is_final[t.from] = false;
    }

    if (layer_map.count(t.to)) {
      if (!is_final[t.to] && l != layer_map[t.to]) {
        cout << "old trans " << prev_map[t.to] << "("
             << layer_map[prev_map[t.to]] << ") -> " << t.to << std::endl;
        cout << "new trans " << t.from << "(" << layer_map[t.from] << ") -> "
             << t.to << std::endl;

        cout << "s UNSUPPORTED"
             << _ID_(": MDD non final node on multiple layers") << " " << t.to
             << " " << is_final[t.to] << " "
             << "\n";
        exit(1);
      }
    } else {
      is_final[t.to] = true;
    }

    prev_map[t.to] = t.from;
    layer_map[t.to] = l;
    if (l >= num_layers)
      num_layers = l + 1;
  }

  vector<int> layer_size;
  layer_size.resize(num_layers, 0);

  // compute the size of every layer and compute the map node -> value (from 1
  // to layer_size because value 0 is for the final node)
  for (const auto &l : layer_map) {
    if (is_final[l.first])
      node_map[l.first] = 0;
    else
      node_map[l.first] = ++layer_size[l.second];
  }

  // last layer should contain only the final state
  if (layer_size[num_layers - 1] > 0) {
    cout << "s UNSUPPORTED" << _ID_(": MDD non final node on last layer")
         << "\n";
    exit(1);
  }
  layer_size.pop_back();
  --num_layers;

  // initialise the state variables, one value per node of the layer, plus the
  // value 0 for final node
  VarArray state_var;
  for (auto m : layer_size) {
    state_var.add(Variable(0, m));
    // std::cout << state_var.back().get_domain() << std::endl;
  }

#ifdef _DEBUG_MDD
  std::cout << "MDD cons: " << scope.size << " vars, " << layer_size.size()
            << " layers " << std::endl;
// exit(1);
#endif

  // create the tables for layer transitions
  VarArray S;
  vector<TableExpression *> tables;
  for (size_t i = 0; i < scope.size; ++i) {
    if (i > 0)
      S.add(state_var[i - 1]); // first transition is binary because it always
                               // come from the root
    S.add(scope[i]);
    if ((int)i < num_layers)
      S.add(state_var[i]); // last transition is binary because it always goes
                           // to the final node

    tables.push_back(new TableExpression(S));
    S.clear();
  }

  int tuple[3];
  for (const auto &t : transitions) {
    auto cons{0};

    if (t.from == root) {
      tuple[0] = t.val;
      tuple[1] = node_map[t.to];

      if (tables[cons]->children.size != 2) {
        std::cout << "ERROR ROOT" << std::endl;
        exit(1);
      }

      // std::cout << "add a tuple of size 2 (root)" ;

    } else {
      cons = layer_map[t.from] + 1;

      if (cons == num_layers) {
        if (!is_final[t.to]) {
          cout << "s UNSUPPORTED" << _ID_(": MDD dead-end transition") << "\n";
          exit(1);
        }
        tuple[0] = node_map[t.from];
        tuple[1] = t.val;

        // std::cout << "add a tuple of size 2 (final)" ;

        if (tables[cons]->children.size != 2) {
          std::cout << "ERROR FINAL" << std::endl;
          exit(1);
        }

      } else {
        tuple[0] = node_map[t.from];
        tuple[1] = t.val;
        tuple[2] = node_map[t.to];

        // std::cout << "add a tuple of size 3 (trans) " << t.from << "-" <<
        // t.val << "-" << t.to << " / ";
        // std::cout.flush();
        //
        // std::cout << node_map[t.from] ;
        // std::cout.flush();
        // std::cout << "-" << t.val ;
        // std::cout.flush();
        // std::cout << "-" << node_map[t.to];
        // std::cout.flush();

        if (tables[cons]->children.size != 3) {
          std::cout << "ERROR" << std::endl;
          exit(1);
        }
      }
    }

    // for(int i=0; i<tables[cons]->children.size; ++i) {
    // 	std::cout << " " << tuple[i];
    // }

    // std::cout << " to table " << cons << "(" << tables[cons]->children.size
    // << ")\n";

    tables[cons]->add(tuple);
  }

  for (auto tab : tables) {

    // #ifdef _DEBUG_MDD
    // 		std::cout << "add tab: " << tab->children << std::endl;
    //
    // 		for(int i=0; i<tab->tuples->size; ++i) {
    // 			std::cout << (int*)(tab->tuples+i) ;
    // 			for(int j=0; j<tab->children.size; ++j) {
    // 				std::cout << " " << (*(tab->tuples))[i][j]
    // ;
    // 			}
    // 			std::cout << std::endl;
    // 		}
    //
    // #endif
    //

    solver.add(Variable(tab));

    // if(tab->children.size == 3) {
    // 	exit(1);
    // }
  }
}

void XCSP3MistralCallbacks::buildConstraintAlldifferent(
    string id, vector<XVariable *> &list) {
#ifdef _VERBOSE_
  cout << "\n    allDiff constraint" << id << endl;
  cout << "        ";
  displayList(list);
#endif

  VarArray scope;
  getVariables(list, scope);

  solver.add(AllDiff(scope));
}

void XCSP3MistralCallbacks::buildConstraintAlldifferentExcept(
    string id, vector<XVariable *> &list, vector<int> &except) {
#ifdef _VERBOSE_
  cout << "\n    allDiff constraint with exceptions" << id << endl;
  cout << "        ";
  displayList(list);
  cout << "        Exceptions:";
  displayList(except);
#endif

  if (except.size() > 1) {
    cout << "s UNSUPPORTED" << _ID_(": AllDiff except SET") << "\n";
    exit(1);
  }

  VarArray scope;
  getVariables(list, scope);

  solver.add(AllDiffExcept(scope, *begin(except)));
}

void XCSP3MistralCallbacks::buildConstraintAlldifferentList(
    string id, vector<vector<XVariable *>> &lists) {
#ifdef _VERBOSE_
  cout << "\n    allDiff list constraint" << id << endl;
  for (unsigned int i = 0; i < (lists.size() < 4 ? lists.size() : 3); i++) {
    cout << "        ";
    displayList(lists[i]);
  }
#endif

  for (size_t i = 0; i < lists[0].size(); ++i)
    ++initial_degree[id_map[lists[0][i]->id]];
  for (auto lpt1 = begin(lists); lpt1 != end(lists); ++lpt1) {
    for (auto lpt2 = lpt1 + 1; lpt2 != end(lists); ++lpt2) {
      // assert(lpt1->size() == lpt2->size());
      VarArray differences;
      for (size_t i = 0; i < lists[0].size(); ++i) {
        ++initial_degree[id_map[(*lpt2)[i]->id]];
        differences.add(variable[(*lpt1)[i]->id] != variable[(*lpt2)[i]->id]);
      }
      solver.add(Sum(differences) > 0);
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintAlldifferentMatrix(
    string id, vector<vector<XVariable *>> &matrix) {
#ifdef _VERBOSE_
  cout << "\n    allDiff matrix constraint" << id << endl;
  for (unsigned int i = 0; i < matrix.size(); i++) {
    cout << "        ";
    displayList(matrix[i]);
  }
#endif

  VarArray scope;
  for (unsigned int i = 0; i < matrix.size(); i++) {
    getVariables(matrix[i], scope);
    solver.add(AllDiff(scope));
    scope.clear();
    for (unsigned int j = 0; j < matrix.size(); j++) {
      scope.add(variable[matrix[j][i]->id]);
    }
    solver.add(AllDiff(scope));
    scope.clear();
  }
}

void XCSP3MistralCallbacks::buildConstraintAllEqual(string id,
                                                    vector<XVariable *> &list) {
#ifdef _VERBOSE_
  cout << "\n    allEqual constraint" << id << endl;
  cout << "        ";
  displayList(list);
#endif

  ++initial_degree[id_map[list[0]->id]];
  for (size_t i = 1; i < list.size(); ++i) {
    ++initial_degree[id_map[list[i]->id]];
    solver.add(variable[list[i - 1]->id] == variable[list[i]->id]);
  }
}

void XCSP3MistralCallbacks::buildConstraintNotAllEqual(
    string id, vector<XVariable *> &list) {
#ifdef _VERBOSE_
  cout << "\n    not allEqual constraint" << id << endl;
  cout << "        ";
  displayList(list);
#endif

  VarArray differences;

  for (size_t i = 0; i < list.size(); ++i)
    ++initial_degree[id_map[list[i]->id]];

  for (auto xpt = begin(list); xpt != end(list); ++xpt) {
    for (auto ypt = xpt + 1; ypt != end(list); ++ypt) {
      differences.add(variable[(*xpt)->id] != variable[(*ypt)->id]);
    }
  }

  solver.add(Sum(differences) > 0);
}

void XCSP3MistralCallbacks::buildConstraintOrdered(string id,
                                                   vector<XVariable *> &list,
                                                   OrderType order) {
#ifdef _VERBOSE_
  cout << "\n    ordered constraint" << endl;
  string sep;
  if (order == LT)
    sep = " < ";
  if (order == LE)
    sep = " <= ";
  if (order == GT)
    sep = " > ";
  if (order == GE)
    sep = " >= ";
  cout << "        ";
  displayList(list, sep);
#endif

  VarArray X;
  getVariables(list, X);

  for (size_t i = 1; i < X.size; ++i) {
    if (order == LT)
      solver.add(X[i - 1] < X[i]);
    if (order == LE)
      solver.add(X[i - 1] <= X[i]);
    if (order == GT)
      solver.add(X[i - 1] > X[i]);
    if (order == GE)
      solver.add(X[i - 1] >= X[i]);
  }
}

void XCSP3MistralCallbacks::buildConstraintLex(
    string id, vector<vector<XVariable *>> &lists, OrderType order) {
#ifdef _VERBOSE_
  cout << "\n    lex constraint   nb lists: " << lists.size() << endl;
  string sep;
  if (order == LT)
    sep = " < ";
  if (order == LE)
    sep = " <= ";
  if (order == GT)
    sep = " > ";
  if (order == GE)
    sep = " >= ";
  cout << "        operator: " << sep << endl;
  for (unsigned int i = 0; i < lists.size(); i++) {
    cout << "        list " << i << ": ";
    cout << "        ";
    displayList(lists[i], " ");
  }
#endif

  vector<VarArray> scope;
  for (auto l = begin(lists); l != end(lists); ++l) {
    VarArray X;
    getVariables(*l, X);
    scope.push_back(X);
  }

  for (size_t i = 1; i < lists.size(); ++i) {
    if (order == LT)
      solver.add(scope[i - 1] < scope[i]);
    if (order == LE)
      solver.add(scope[i - 1] <= scope[i]);
    if (order == GT)
      solver.add(scope[i - 1] > scope[i]);
    if (order == GE)
      solver.add(scope[i - 1] >= scope[i]);
  }
}

void XCSP3MistralCallbacks::buildConstraintLexMatrix(
    string id, vector<vector<XVariable *>> &matrix, OrderType order) {
#ifdef _VERBOSE_
  cout << "\n    lex matrix constraint   matrix  " << endl;
  string sep;
  if (order == LT)
    sep = " < ";
  if (order == LE)
    sep = " <= ";
  if (order == GT)
    sep = " > ";
  if (order == GE)
    sep = " >= ";

  // for(unsigned int i = 0; i < (matrix.size() < 4 ? matrix.size() : 3); i++) {
  for (unsigned int i = 0; i < matrix.size(); i++) {
    cout << "        ";
    displayList(matrix[i]);
  }
  cout << "        Order " << sep << endl;
#endif

  vector<VarArray> scope;
  for (auto l = begin(matrix); l != end(matrix); ++l) {
    VarArray X;
    getVariables(*l, X);
    scope.push_back(X);
  }

  for (size_t i = 1; i < scope.size(); ++i) {
    if (order == LT)
      solver.add(scope[i - 1] < scope[i]);
    if (order == LE)
      solver.add(scope[i - 1] <= scope[i]);
    if (order == GT)
      solver.add(scope[i - 1] > scope[i]);
    if (order == GE)
      solver.add(scope[i - 1] >= scope[i]);

    // if(order == LT) cout << scope[i-1] << "\n<\n" <<  scope[i] << endl;
    // if(order == LE) cout << scope[i-1] << "\n<=\n" <<  scope[i] << endl;
    // if(order == GT) cout << scope[i-1] << "\n>\n" <<  scope[i] << endl;
    // if(order == GE) cout << scope[i-1] << "\n>=\n" <<  scope[i] << endl;
  }

  scope.clear();

  for (size_t i = 0; i < begin(matrix)->size(); ++i) {
    VarArray X;
    for (size_t j = 0; j < matrix.size(); ++j) {
      X.add(variable[matrix[j][i]->id]);
    }
    scope.push_back(X);
  }

  for (size_t i = 1; i < scope.size(); ++i) {
    if (order == LT)
      solver.add(scope[i - 1] < scope[i]);
    if (order == LE)
      solver.add(scope[i - 1] <= scope[i]);
    if (order == GT)
      solver.add(scope[i - 1] > scope[i]);
    if (order == GE)
      solver.add(scope[i - 1] >= scope[i]);

    // if(order == LT) cout << scope[i-1] << "\n<\n" <<  scope[i] << endl;
    // if(order == LE) cout << scope[i-1] << "\n<=\n" <<  scope[i] << endl;
    // if(order == GT) cout << scope[i-1] << "\n>\n" <<  scope[i] << endl;
    // if(order == GE) cout << scope[i-1] << "\n>=\n" <<  scope[i] << endl;
  }
}

void XCSP3MistralCallbacks::buildConstraintSum(string id,
                                               vector<XVariable *> &list,
                                               vector<int> &coeffs,
                                               XCondition &cond) {
#ifdef _VERBOSE_
  cout << "\n        sum constraint: ";
  if (list.size() > 8) {
    for (int i = 0; i < 3; i++)
      cout << (coeffs.size() == 0 ? 1 : coeffs[i]) << "*" << *(list[i]) << " ";
    cout << " ... ";
    for (unsigned int i = list.size() - 4; i < list.size(); i++)
      cout << (coeffs.size() == 0 ? 1 : coeffs[i]) << "*" << *(list[i]) << " ";
  } else {
    for (unsigned int i = 0; i < list.size(); i++)
      cout << (coeffs.size() == 0 ? 1 : coeffs[i]) << "*" << *(list[i]) << " ";
  }
  cout << cond << endl;
#endif

  vector<Variable> scope;
  for (auto x : list) {
    ++initial_degree[id_map[x->id]];
    scope.push_back(variable[x->id]);
  }

  if (cond.operandType == VARIABLE) {
    Variable total = variable[cond.var];

    if (cond.op == EQ) {

      // std::cout << "HERE " << total << " in " << total.get_domain() << "\n" ;

      solver.add(Sum(scope, coeffs, total));
    } else if (cond.op == NE) {
      solver.add(Sum(scope, coeffs) != total);
    } else if (cond.op == LE) {
      solver.add(Sum(scope, coeffs) <= total);
    } else if (cond.op == LT) {
      solver.add(Sum(scope, coeffs) < total);
    } else if (cond.op == GE) {
      solver.add(Sum(scope, coeffs) >= total);
    } else if (cond.op == GT) {
      solver.add(Sum(scope, coeffs) > total);
    }

  } else if (cond.operandType == INTERVAL) {
    if (cond.op != IN) {
      cout << "s UNSUPPORTED" << _ID_(": IN non-interval") << "\n";
      exit(1);
    }

    solver.add(Sum(scope, coeffs, cond.min, cond.max));

  } else {

    if (cond.op == EQ) {
      solver.add(Sum(scope, coeffs, cond.val, cond.val));
    } else if (cond.op == NE) {
      Variable cst(cond.val);
      solver.add(Sum(scope, coeffs) != cst);
    } else if (cond.op == LE) {
      solver.add(Sum(scope, coeffs, -INFTY, cond.val));
    } else if (cond.op == LT) {
      solver.add(Sum(scope, coeffs, -INFTY, cond.val - 1));
    } else if (cond.op == GE) {
      solver.add(Sum(scope, coeffs, cond.val, INFTY));
    } else if (cond.op == GT) {
      solver.add(Sum(scope, coeffs, cond.val + 1, INFTY));
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintSum(string id, vector<Tree *> &trees,
                                               XCondition &cond) {

// Not sure if this is correct!
#ifdef _VERBOSE_
  cout << "\n        sum constraint with tree:";

  cout << endl;
#endif

  vector<Variable> scope;
  for (auto x : trees) {
    // initial_degree ???
    //++initial_degree[x->id];
    scope.push_back(postExpression(x->root, true).get_var());
  }

  if (cond.operandType == VARIABLE) {
    Variable total = variable[cond.var];

    if (cond.op == EQ) {
      solver.add(Sum(scope, total));
    } else if (cond.op == NE) {
      solver.add(Sum(scope) != total);
    } else if (cond.op == LE) {
      solver.add(Sum(scope) <= total);
    } else if (cond.op == LT) {
      solver.add(Sum(scope) < total);
    } else if (cond.op == GE) {
      solver.add(Sum(scope) >= total);
    } else if (cond.op == GT) {
      solver.add(Sum(scope) > total);
    }

  } else if (cond.operandType == INTERVAL) {
    if (cond.op != IN) {
      cout << "s UNSUPPORTED" << _ID_(": IN non-interval") << "\n";
      exit(1);
    }

    solver.add(Sum(scope, cond.min, cond.max));

  } else {

    if (cond.op == EQ) {
      solver.add(Sum(scope, cond.val, cond.val));
    } else if (cond.op == NE) {
      Variable cst(cond.val);
      solver.add(Sum(scope) != cst);
    } else if (cond.op == LE) {
      solver.add(Sum(scope, -INFTY, cond.val));
    } else if (cond.op == LT) {
      solver.add(Sum(scope, -INFTY, cond.val - 1));
    } else if (cond.op == GE) {
      solver.add(Sum(scope, cond.val, INFTY));
    } else if (cond.op == GT) {
      solver.add(Sum(scope, cond.val + 1, INFTY));
    }
  }
}
//
// void XCSP3MistralCallbacks::buildConstraintSum(string id, vector<Tree *>
// &trees,
//                                vector<int> &coefs,
//                                XCondition &cond) {
//
//                                }

void XCSP3MistralCallbacks::buildConstraintSum(string id, vector<Tree *> &trees,
                                               vector<int> &coeffs,
                                               XCondition &cond) {

// Not sure if this is correct!
#ifdef _VERBOSE_
  cout << "\n        sum constraint with tree:";

  cout << endl;
#endif

  vector<Variable> scope;
  for (auto x : trees) {
    // initial_degree ???
    //++initial_degree[x->id];
    scope.push_back(postExpression(x->root, true).get_var());
  }

  if (cond.operandType == VARIABLE) {
    Variable total = variable[cond.var];

    if (cond.op == EQ) {
      solver.add(Sum(scope, coeffs, total));
    } else if (cond.op == NE) {
      solver.add(Sum(scope, coeffs) != total);
    } else if (cond.op == LE) {
      solver.add(Sum(scope, coeffs) <= total);
    } else if (cond.op == LT) {
      solver.add(Sum(scope, coeffs) < total);
    } else if (cond.op == GE) {
      solver.add(Sum(scope, coeffs) >= total);
    } else if (cond.op == GT) {
      solver.add(Sum(scope, coeffs) > total);
    }

  } else if (cond.operandType == INTERVAL) {
    if (cond.op != IN) {
      cout << "s UNSUPPORTED" << _ID_(": IN non-interval") << "\n";
      exit(1);
    }

    solver.add(Sum(scope, coeffs, cond.min, cond.max));

  } else {

    if (cond.op == EQ) {
      solver.add(Sum(scope, coeffs, cond.val, cond.val));
    } else if (cond.op == NE) {
      Variable cst(cond.val);
      solver.add(Sum(scope, coeffs) != cst);
    } else if (cond.op == LE) {
      solver.add(Sum(scope, coeffs, -INFTY, cond.val));
    } else if (cond.op == LT) {
      solver.add(Sum(scope, coeffs, -INFTY, cond.val - 1));
    } else if (cond.op == GE) {
      solver.add(Sum(scope, coeffs, cond.val, INFTY));
    } else if (cond.op == GT) {
      solver.add(Sum(scope, coeffs, cond.val + 1, INFTY));
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintSum(string id,
                                               vector<XVariable *> &list,
                                               XCondition &cond) {
#ifdef _VERBOSE_
  cout << "\n        unweighted sum constraint:";
  cout << "        ";
  displayList(list, "+");
  cout << cond << endl;
#endif

  vector<Variable> scope;
  for (auto x : list) {
    ++initial_degree[id_map[x->id]];
    scope.push_back(variable[x->id]);
  }

  if (cond.operandType == VARIABLE) {
    Variable total = variable[cond.var];

    if (cond.op == EQ) {
      solver.add(Sum(scope, total));
    } else if (cond.op == NE) {
      solver.add(Sum(scope) != total);
    } else if (cond.op == LE) {
      solver.add(Sum(scope) <= total);
    } else if (cond.op == LT) {
      solver.add(Sum(scope) < total);
    } else if (cond.op == GE) {
      solver.add(Sum(scope) >= total);
    } else if (cond.op == GT) {
      solver.add(Sum(scope) > total);
    }

  } else if (cond.operandType == INTERVAL) {
    if (cond.op != IN) {
      cout << "s UNSUPPORTED" << _ID_(": IN non-interval") << "\n";
      exit(1);
    }

    solver.add(Sum(scope, cond.min, cond.max));

  } else {

    if (cond.op == EQ) {
      solver.add(Sum(scope, cond.val, cond.val));
    } else if (cond.op == NE) {
      Variable cst(cond.val);
      solver.add(Sum(scope) != cst);
    } else if (cond.op == LE) {
      solver.add(Sum(scope, -INFTY, cond.val));
    } else if (cond.op == LT) {
      solver.add(Sum(scope, -INFTY, cond.val - 1));
    } else if (cond.op == GE) {
      solver.add(Sum(scope, cond.val, INFTY));
    } else if (cond.op == GT) {
      solver.add(Sum(scope, cond.val + 1, INFTY));
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintSum(string id,
                                               vector<XVariable *> &list,
                                               vector<XVariable *> &coeffs,
                                               XCondition &cond) {
#ifdef _VERBOSE_
  cout << "\n        scalar sum constraint:";
  if (list.size() > 8) {
    for (int i = 0; i < 3; i++)
      cout << coeffs[i]->id << "*" << *(list[i]) << " ";
    cout << " ... ";
    for (unsigned int i = list.size() - 4; i < list.size(); i++)
      cout << coeffs[i]->id << "*" << *(list[i]) << " ";
  } else {
    for (unsigned int i = 0; i < list.size(); i++)
      cout << coeffs[i]->id << "*" << *(list[i]) << " ";
  }
  cout << cond << endl;
#endif

  vector<Variable> scope;
  for (unsigned int i = 0; i < list.size(); i++) {
    ++initial_degree[id_map[list[i]->id]];
    ++initial_degree[id_map[coeffs[i]->id]];
    scope.push_back(variable[list[i]->id] * variable[coeffs[i]->id]);
  }

  if (cond.operandType == VARIABLE) {
    Variable total = variable[cond.var];

    if (cond.op == EQ) {
      solver.add(Sum(scope, total));
    } else if (cond.op == NE) {
      solver.add(Sum(scope) != total);
    } else if (cond.op == LE) {
      solver.add(Sum(scope) <= total);
    } else if (cond.op == LT) {
      solver.add(Sum(scope) < total);
    } else if (cond.op == GE) {
      solver.add(Sum(scope) >= total);
    } else if (cond.op == GT) {
      solver.add(Sum(scope) > total);
    }

  } else if (cond.operandType == INTERVAL) {
    if (cond.op != IN) {
      cout << "s UNSUPPORTED" << _ID_(": IN non-interval") << "\n";
      exit(1);
    }

    solver.add(Sum(scope, cond.min, cond.max));

  } else {

    if (cond.op == EQ) {
      solver.add(Sum(scope, cond.val, cond.val));
    } else if (cond.op == NE) {
      Variable cst(cond.val);
      solver.add(Sum(scope) != cst);
    } else if (cond.op == LE) {
      solver.add(Sum(scope, -INFTY, cond.val));
    } else if (cond.op == LT) {
      solver.add(Sum(scope, -INFTY, cond.val - 1));
    } else if (cond.op == GE) {
      solver.add(Sum(scope, cond.val, INFTY));
    } else if (cond.op == GT) {
      solver.add(Sum(scope, cond.val + 1, INFTY));
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintAtMost(string id,
                                                  vector<XVariable *> &list,
                                                  int value, int k) {
#ifdef _VERBOSE_
  cout << "\n    AtMost constraint: val=" << value << " k=" << k << endl;
  cout << "        ";
  displayList(list);
#endif

  VarArray scope;
  getVariables(list, scope);

  VarArray equalities;
  for (auto X : scope) {
    equalities.add(X == value);
  }

  solver.add(Sum(equalities) <= k);
}

void XCSP3MistralCallbacks::buildConstraintAtLeast(string id,
                                                   vector<XVariable *> &list,
                                                   int value, int k) {
#ifdef _VERBOSE_
  cout << "\n    Atleast constraint: val=" << value << " k=" << k << endl;
  cout << "        ";
  displayList(list);
#endif

  XCondition cond;
  cond.operandType = INTEGER;
  cond.op = GE;
  cond.val = k;
  vector<int> values;
  values.push_back(value);
  return buildConstraintCount(id, list, values, cond);
}

void XCSP3MistralCallbacks::buildConstraintExactlyK(string id,
                                                    vector<XVariable *> &list,
                                                    int value, int k) {
#ifdef _VERBOSE_
  cout << "\n    Exactly constraint: val=" << value << " k=" << k << endl;
  cout << "        ";
  displayList(list);
#endif

  XCondition cond;
  cond.operandType = INTEGER;
  cond.op = EQ;
  cond.val = k;
  vector<int> values;
  values.push_back(value);
  return buildConstraintCount(id, list, values, cond);
}

void XCSP3MistralCallbacks::buildConstraintAmong(string id,
                                                 vector<XVariable *> &list,
                                                 vector<int> &values, int k) {
#ifdef _VERBOSE_
  cout << "\n    Among constraint: k=" << k << endl;
  cout << "        ";
  displayList(list);
  cout << "        values:";
  displayList(values);
#endif

  XCondition cond;
  cond.operandType = INTEGER;
  cond.op = EQ;
  cond.val = k;
  return buildConstraintCount(id, list, values, cond);
}

void XCSP3MistralCallbacks::buildConstraintExactlyVariable(
    string id, vector<XVariable *> &list, int value, XVariable *x) {
#ifdef _VERBOSE_
  cout << "\n    Exactly Variable constraint: val=" << value
       << " variable=" << *x << endl;
  cout << "        ";
  displayList(list);
#endif

  XCondition cond;
  cond.operandType = VARIABLE;
  cond.op = EQ;
  cond.var = x->id;
  vector<int> values;
  values.push_back(value);
  return buildConstraintCount(id, list, values, cond);
}

void XCSP3MistralCallbacks::buildConstraintCount(string id,
                                                 vector<XVariable *> &list,
                                                 vector<int> &values,
                                                 XCondition &cond) {
#ifdef _VERBOSE_
  cout << "\n    count constraint" << endl;
  cout << "        ";
  displayList(list);
  cout << "        values: ";
  cout << "        ";
  displayList(values);
  cout << "        condition: " << cond << endl;
#endif

  VarArray scope;
  getVariables(list, scope);

  if (cond.operandType == VARIABLE) {
    Variable nocc = variable[cond.var];

    if (values.size() == 1) {
      int val = *begin(values);
      if (cond.op == EQ) {
        solver.add(Occurrence(scope, val) == nocc);
      } else if (cond.op == NE) {
        solver.add(Occurrence(scope, val) != nocc);
      } else if (cond.op == LE) {
        solver.add(Occurrence(scope, val) <= nocc);
      } else if (cond.op == LT) {
        solver.add(Occurrence(scope, val) < nocc);
      } else if (cond.op == GE) {
        solver.add(Occurrence(scope, val) >= nocc);
      } else if (cond.op == GT) {
        solver.add(Occurrence(scope, val) > nocc);
      }
    } else {
      if (cond.op == EQ) {
        solver.add(Occurrence(scope, values) == nocc);
      } else if (cond.op == NE) {
        solver.add(Occurrence(scope, values) != nocc);
      } else if (cond.op == LE) {
        solver.add(Occurrence(scope, values) <= nocc);
      } else if (cond.op == LT) {
        solver.add(Occurrence(scope, values) < nocc);
      } else if (cond.op == GE) {
        solver.add(Occurrence(scope, values) >= nocc);
      } else if (cond.op == GT) {
        solver.add(Occurrence(scope, values) > nocc);
      }
    }

  } else if (cond.operandType == INTERVAL) {
    if (cond.op != IN) {
      cout << "s UNSUPPORTED" << _ID_(": IN non-interval") << "\n";
      exit(1);
    }

    if (values.size() == 1) {
      solver.add(Occurrence(scope, *begin(values), cond.min, cond.max));
    } else {
      solver.add(Occurrence(scope, values, cond.min, cond.max));
    }

  } else {

    if (values.size() == 1) {
      int val = *begin(values);

      if (cond.op == EQ) {
        solver.add(Occurrence(scope, val, cond.val, cond.val));
      } else if (cond.op == NE) {
        Variable cst(cond.val);
        solver.add(Occurrence(scope, val) != cst);
      } else if (cond.op == LE) {
        solver.add(Occurrence(scope, val, -INFTY, cond.val));
      } else if (cond.op == LT) {
        solver.add(Occurrence(scope, val, -INFTY, cond.val - 1));
      } else if (cond.op == GE) {
        solver.add(Occurrence(scope, val, cond.val, INFTY));
      } else if (cond.op == GT) {
        solver.add(Occurrence(scope, val, cond.val + 1, INFTY));
      }
    } else {
      if (cond.op == EQ) {
        solver.add(Occurrence(scope, values, cond.val, cond.val));
      } else if (cond.op == NE) {
        Variable cst(cond.val);
        solver.add(Occurrence(scope, values) != cst);
      } else if (cond.op == LE) {
        solver.add(Occurrence(scope, values, -INFTY, cond.val));
      } else if (cond.op == LT) {
        solver.add(Occurrence(scope, values, -INFTY, cond.val - 1));
      } else if (cond.op == GE) {
        solver.add(Occurrence(scope, values, cond.val, INFTY));
      } else if (cond.op == GT) {
        solver.add(Occurrence(scope, values, cond.val + 1, INFTY));
      }
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintCount(string id,
                                                 vector<XVariable *> &list,
                                                 vector<XVariable *> &values,
                                                 XCondition &cond) {
#ifdef _VERBOSE_
  cout << "\n    count constraint" << endl;
  cout << "        ";
  displayList(list);
  cout << "        values: ";
  displayList(values);
  cout << "        condition: " << cond << endl;
#endif

  VarArray vars;
  getVariables(list, vars);

  if (values.size() > 1)
    throw runtime_error(
        "s UNSUPPORTED: Count with several values represented as variables");

  VarArray vals;
  getVariables(values, vals);

  VarArray count;
  for (auto x : vars) {
    for (auto v : vals) {
      count.add(x == v);
    }
  }

  if (cond.operandType == VARIABLE) {
    Variable nocc = variable[cond.var];

    if (cond.op == EQ) {
      solver.add(Sum(count) == nocc);
    } else if (cond.op == NE) {
      solver.add(Sum(count) != nocc);
    } else if (cond.op == LE) {
      solver.add(Sum(count) <= nocc);
    } else if (cond.op == LT) {
      solver.add(Sum(count) < nocc);
    } else if (cond.op == GE) {
      solver.add(Sum(count) >= nocc);
    } else if (cond.op == GT) {
      solver.add(Sum(count) > nocc);
    }

  } else if (cond.operandType == INTERVAL) {
    Variable C = Sum(count);
    if (cond.op == NE) {
      solver.add((C <= cond.min) || (C > cond.max));
    } else if (cond.op == IN) {
      solver.add(C >= cond.min);
      solver.add(C <= cond.max);
    } else {
      throw runtime_error(
          "s UNSUPPORTED: Count with conditions other than IN or NE");
    }

  } else {
    if (cond.op == EQ) {
      solver.add(Sum(count) == cond.val);
    } else if (cond.op == NE) {
      Variable cst(cond.val);
      solver.add(Sum(count) != cst);
    } else if (cond.op == LE) {
      solver.add(Sum(count) <= cond.val);
    } else if (cond.op == LT) {
      solver.add(Sum(count) < cond.val);
    } else if (cond.op == GE) {
      solver.add(Sum(count) >= cond.val);
    } else if (cond.op == GT) {
      solver.add(Sum(count) > cond.val);
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintCount(string id,
                                                 vector<Tree *> &trees,
                                                 vector<int> &values,
                                                 XCondition &cond) {
#ifdef _VERBOSE_
  cout << "\n    count constraint" << endl;
  cout << "        ";
  displayList(trees);
  cout << "        values: ";
  cout << "        ";
  displayList(values);
  cout << "        condition: " << cond << endl;
#endif

  VarArray scope;
  for (auto x : trees) {
    scope.push_back(postExpression(x->root, true).get_var());
  }

  // VarArray scope;
  // getVariables(list, scope);

  if (cond.operandType == VARIABLE) {
    Variable nocc = variable[cond.var];

    if (values.size() == 1) {
      int val = *begin(values);
      if (cond.op == EQ) {
        solver.add(Occurrence(scope, val) == nocc);
      } else if (cond.op == NE) {
        solver.add(Occurrence(scope, val) != nocc);
      } else if (cond.op == LE) {
        solver.add(Occurrence(scope, val) <= nocc);
      } else if (cond.op == LT) {
        solver.add(Occurrence(scope, val) < nocc);
      } else if (cond.op == GE) {
        solver.add(Occurrence(scope, val) >= nocc);
      } else if (cond.op == GT) {
        solver.add(Occurrence(scope, val) > nocc);
      }
    } else {
      if (cond.op == EQ) {
        solver.add(Occurrence(scope, values) == nocc);
      } else if (cond.op == NE) {
        solver.add(Occurrence(scope, values) != nocc);
      } else if (cond.op == LE) {
        solver.add(Occurrence(scope, values) <= nocc);
      } else if (cond.op == LT) {
        solver.add(Occurrence(scope, values) < nocc);
      } else if (cond.op == GE) {
        solver.add(Occurrence(scope, values) >= nocc);
      } else if (cond.op == GT) {
        solver.add(Occurrence(scope, values) > nocc);
      }
    }

  } else if (cond.operandType == INTERVAL) {
    if (cond.op != IN) {
      cout << "s UNSUPPORTED" << _ID_(": IN non-interval") << "\n";
      exit(1);
    }

    if (values.size() == 1) {
      solver.add(Occurrence(scope, *begin(values), cond.min, cond.max));
    } else {
      solver.add(Occurrence(scope, values, cond.min, cond.max));
    }

  } else {

    if (values.size() == 1) {
      int val = *begin(values);

      if (cond.op == EQ) {
        solver.add(Occurrence(scope, val, cond.val, cond.val));
      } else if (cond.op == NE) {
        Variable cst(cond.val);
        solver.add(Occurrence(scope, val) != cst);
      } else if (cond.op == LE) {
        solver.add(Occurrence(scope, val, -INFTY, cond.val));
      } else if (cond.op == LT) {
        solver.add(Occurrence(scope, val, -INFTY, cond.val - 1));
      } else if (cond.op == GE) {
        solver.add(Occurrence(scope, val, cond.val, INFTY));
      } else if (cond.op == GT) {
        solver.add(Occurrence(scope, val, cond.val + 1, INFTY));
      }
    } else {
      if (cond.op == EQ) {
        solver.add(Occurrence(scope, values, cond.val, cond.val));
      } else if (cond.op == NE) {
        Variable cst(cond.val);
        solver.add(Occurrence(scope, values) != cst);
      } else if (cond.op == LE) {
        solver.add(Occurrence(scope, values, -INFTY, cond.val));
      } else if (cond.op == LT) {
        solver.add(Occurrence(scope, values, -INFTY, cond.val - 1));
      } else if (cond.op == GE) {
        solver.add(Occurrence(scope, values, cond.val, INFTY));
      } else if (cond.op == GT) {
        solver.add(Occurrence(scope, values, cond.val + 1, INFTY));
      }
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintCount(string id,
                                                 vector<Tree *> &trees,
                                                 vector<XVariable *> &values,
                                                 XCondition &cond) {
#ifdef _VERBOSE_
  cout << "\n    count constraint" << endl;
  cout << "        ";
  displayList(trees);
  cout << "        values: ";
  displayList(values);
  cout << "        condition: " << cond << endl;
#endif

  VarArray vars;
  for (auto x : trees) {
    vars.add(postExpression(x->root, true).get_var());
  }

  if (values.size() > 1)
    throw runtime_error(
        "s UNSUPPORTED: Count with several values represented as variables");

  VarArray vals;
  getVariables(values, vals);

  VarArray count;
  for (auto x : vars) {
    for (auto v : vals) {
      count.add(x == v);
    }
  }

  if (cond.operandType == VARIABLE) {
    Variable nocc = variable[cond.var];

    if (cond.op == EQ) {
      solver.add(Sum(count) == nocc);
    } else if (cond.op == NE) {
      solver.add(Sum(count) != nocc);
    } else if (cond.op == LE) {
      solver.add(Sum(count) <= nocc);
    } else if (cond.op == LT) {
      solver.add(Sum(count) < nocc);
    } else if (cond.op == GE) {
      solver.add(Sum(count) >= nocc);
    } else if (cond.op == GT) {
      solver.add(Sum(count) > nocc);
    }

  } else if (cond.operandType == INTERVAL) {
    Variable C = Sum(count);
    if (cond.op == NE) {
      solver.add((C <= cond.min) || (C > cond.max));
    } else if (cond.op == IN) {
      solver.add(C >= cond.min);
      solver.add(C <= cond.max);
    } else {
      throw runtime_error(
          "s UNSUPPORTED: Count with conditions other than IN or NE");
    }

  } else {
    if (cond.op == EQ) {
      solver.add(Sum(count) == cond.val);
    } else if (cond.op == NE) {
      Variable cst(cond.val);
      solver.add(Sum(count) != cst);
    } else if (cond.op == LE) {
      solver.add(Sum(count) <= cond.val);
    } else if (cond.op == LT) {
      solver.add(Sum(count) < cond.val);
    } else if (cond.op == GE) {
      solver.add(Sum(count) >= cond.val);
    } else if (cond.op == GT) {
      solver.add(Sum(count) > cond.val);
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintNValues(string id,
                                                   vector<XVariable *> &list,
                                                   vector<int> &except,
                                                   XCondition &cond) {
#ifdef _VERBOSE_
  cout << "\n    NValues with exceptions constraint" << endl;
  cout << "        ";
  displayList(list);
  cout << "        exceptions: ";
  displayList(except);
  cout << "        condition:" << cond << endl;
#endif

  VarArray scope;
  getVariables(list, scope);

  int min = INFTY;
  int max = -INFTY;

  for (auto X : scope) {
    if (min > X.get_min()) {
      min = X.get_min();
    }
    if (max < X.get_max()) {
      max = X.get_max();
    }
  }

  int maxe = *max_element(begin(except), end(except));
  int mine = *min_element(begin(except), end(except));

  BitSet exceptions(mine, maxe, BitSet::empt);

  for (auto e : except) {
    exceptions.add(e);
  }

  VarArray used;
  for (int i = min; i <= max; ++i) {
    if (!exceptions.contain(i))
      used.add(Occurrence(scope, i));
  }

  if (cond.operandType == VARIABLE) {
    Variable nval = variable[cond.var];

    if (cond.op == EQ) {
      solver.add(BoolSum(used) == nval);
    } else if (cond.op == NE) {
      solver.add(BoolSum(used) != nval);
    } else if (cond.op == LE) {
      solver.add(BoolSum(used) <= nval);
    } else if (cond.op == LT) {
      solver.add(BoolSum(used) < nval);
    } else if (cond.op == GE) {
      solver.add(BoolSum(used) >= nval);
    } else if (cond.op == GT) {
      solver.add(BoolSum(used) > nval);
    }

  } else if (cond.operandType == INTERVAL) {
    if (cond.op != IN) {
      cout << "s UNSUPPORTED" << _ID_(": IN non-interval") << "\n";
      exit(1);
    }

    solver.add(BoolSum(used, cond.min, cond.max));
  } else {

    if (cond.op == EQ) {
      solver.add(BoolSum(used, cond.val, cond.val));
    } else if (cond.op == NE) {
      Variable cst(cond.val);
      solver.add(BoolSum(used) != cst);
    } else if (cond.op == LE) {
      solver.add(BoolSum(used, -INFTY, cond.val));
    } else if (cond.op == LT) {
      solver.add(BoolSum(used, -INFTY, cond.val - 1));
    } else if (cond.op == GE) {
      solver.add(BoolSum(used, cond.val, INFTY));
    } else if (cond.op == GT) {
      solver.add(BoolSum(used, cond.val + 1, INFTY));
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintNValues(string id,
                                                   vector<XVariable *> &list,
                                                   XCondition &cond) {
#ifdef _VERBOSE_
  cout << "\n    NValues  constraint" << endl;
  cout << "        ";
  displayList(list);
  cout << "        condition:" << cond << endl;
#endif

  VarArray scope;
  getVariables(list, scope);

  int min = INFTY;
  int max = -INFTY;

  for (auto X : scope) {
    if (min > X.get_min()) {
      min = X.get_min();
    }
    if (max < X.get_max()) {
      max = X.get_max();
    }
  }

  VarArray used;
  for (int i = min; i <= max; ++i) {
    used.add(Occurrence(scope, i));
  }

  if (cond.operandType == VARIABLE) {
    Variable nval = variable[cond.var];

    if (cond.op == EQ) {
      solver.add(BoolSum(used) == nval);
    } else if (cond.op == NE) {
      solver.add(BoolSum(used) != nval);
    } else if (cond.op == LE) {
      solver.add(BoolSum(used) <= nval);
    } else if (cond.op == LT) {
      solver.add(BoolSum(used) < nval);
    } else if (cond.op == GE) {
      solver.add(BoolSum(used) >= nval);
    } else if (cond.op == GT) {
      solver.add(BoolSum(used) > nval);
    }

  } else if (cond.operandType == INTERVAL) {
    if (cond.op != IN) {
      cout << "s UNSUPPORTED" << _ID_(": IN non-interval") << "\n";
      exit(1);
    }

    solver.add(BoolSum(used, cond.min, cond.max));
  } else {

    if (cond.op == EQ) {
      solver.add(BoolSum(used, cond.val, cond.val));
    } else if (cond.op == NE) {
      Variable cst(cond.val);
      solver.add(BoolSum(used) != cst);
    } else if (cond.op == LE) {
      solver.add(BoolSum(used, -INFTY, cond.val));
    } else if (cond.op == LT) {
      solver.add(BoolSum(used, -INFTY, cond.val - 1));
    } else if (cond.op == GE) {
      solver.add(BoolSum(used, cond.val, INFTY));
    } else if (cond.op == GT) {
      solver.add(BoolSum(used, cond.val + 1, INFTY));
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintCardinality(
    string id, vector<XVariable *> &list, vector<int> values,
    vector<int> &occurs, bool closed) {
#ifdef _VERBOSE_
  cout << "\n    Cardinality constraint (int values, int occurs)  constraint "
          "closed: "
       << closed << endl;
  cout << "        ";
  displayList(list);
  cout << "        values:";
  displayList(values);
  cout << "        occurs:";
  displayList(occurs);
#endif

  VarArray scope;
  getVariables(list, scope);

  if (closed)
    for (auto X : scope)
      solver.add(Member(X, values));

  if (values.back() - (*begin(values)) + 1 == (int)(values.size())) {
    solver.add(Occurrences(scope, values, occurs, occurs));
  } else {
    for (size_t i = 0; i < values.size(); ++i) {
      solver.add(Occurrence(scope, values[i]) == occurs[i]);
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintCardinality(
    string id, vector<XVariable *> &list, vector<int> values,
    vector<XVariable *> &occurs, bool closed) {
#ifdef _VERBOSE_
  cout << "\n    Cardinality constraint (int values, var occurs)  constraint "
          "closed: "
       << closed << endl;
  cout << "        ";
  displayList(list);
  cout << "        values:";
  displayList(values);
  cout << "        occurs:";
  displayList(occurs);
#endif

  VarArray scope;
  getVariables(list, scope);
  VarArray count;
  getVariables(occurs, count);

  for (size_t i = 0; i < values.size(); ++i) {
    solver.add(Occurrence(scope, values[i]) == count[i]);
  }

  if (closed) {
    solver.add(Sum(count, scope.size, scope.size));
    for (auto X : scope)
      solver.add(Member(X, values));
  }

  // solver.consolidate();
  //
  // std::cout << solver << std::endl;
  //
  // exit(1);
}

void XCSP3MistralCallbacks::buildConstraintCardinality(
    string id, vector<XVariable *> &list, vector<int> values,
    vector<XInterval> &occurs, bool closed) {
#ifdef _VERBOSE_
  cout << "\n    Cardinality constraint (int values, interval occurs)  "
          "constraint closed: "
       << closed << endl;
  cout << "        ";
  displayList(list);
  cout << "        values:";
  displayList(values);
  cout << "        occurs:";
  displayList(occurs);
#endif

  VarArray scope;
  getVariables(list, scope);

  vector<int> lb;
  vector<int> ub;

  int total_lb = 0;
  int total_ub = 0;
  for (auto I : occurs) {
    total_lb += I.min;
    total_ub += I.max;

    lb.push_back(I.min);
    ub.push_back(I.max);
  }

  if (total_lb > (int)(scope.size)) {
    solver.fail();
  }

  if (closed) {
    if (total_ub < (int)(scope.size)) {
      solver.fail();
    }
    for (auto X : scope)
      solver.add(Member(X, values));
  }

  solver.add(Occurrences(scope, values, lb, ub));
}

void XCSP3MistralCallbacks::buildConstraintCardinality(
    string id, vector<XVariable *> &list, vector<XVariable *> values,
    vector<int> &occurs, bool closed) {
#ifdef _VERBOSE_
  cout << "\n    Cardinality constraint (var values, int occurs)  constraint "
          "closed: "
       << closed << endl;
  cout << "        ";
  displayList(list);
  cout << "        values:";
  displayList(values);
  cout << "        occurs:";
  displayList(occurs);
#endif

  VarArray scope;
  getVariables(list, scope);
  VarArray vals;
  getVariables(values, vals);

  solver.add(AllDiff(vals));

  int total = 0;
  for (size_t i = 0; i < values.size(); ++i) {
    total += occurs[i];
    solver.add(Occurrence(scope, vals[i]) == occurs[i]);
  }

  if (total > (int)(scope.size)) {
    solver.fail();
  }

  if (closed || total == (int)(scope.size))
    for (auto X : scope) {
      VarArray M;
      for (auto V : vals)
        M.add(Member(X, V));
      solver.add(Sum(M) > 0);
    }
  else {
    // TODO: set the number of occurrences of values not in values to scope.size
    // - total
  }
}

void XCSP3MistralCallbacks::buildConstraintCardinality(
    string id, vector<XVariable *> &list, vector<XVariable *> values,
    vector<XVariable *> &occurs, bool closed) {
#ifdef _VERBOSE_
  cout << "\n    Cardinality constraint (var values, var occurs)  constraint "
          "closed: "
       << closed << endl;
  cout << "        ";
  displayList(list);
  cout << "        values:";
  displayList(values);
  cout << "        occurs:";
  displayList(occurs);
#endif

  VarArray scope;
  getVariables(list, scope);
  VarArray count;
  getVariables(occurs, count);
  VarArray vals;
  getVariables(values, vals);

  for (size_t i = 0; i < values.size(); ++i) {
    solver.add(Occurrence(scope, vals[i]) == count[i]);
  }

  solver.add(AllDiff(vals));

  if (closed) {
    solver.add(Sum(count, scope.size, scope.size));
    for (auto X : scope) {
      VarArray M;
      for (auto V : vals)
        M.add(Member(X, V));
      solver.add(Sum(M) > 0);
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintCardinality(
    string id, vector<XVariable *> &list, vector<XVariable *> values,
    vector<XInterval> &occurs, bool closed) {
#ifdef _VERBOSE_
  cout << "\n    Cardinality constraint (var values, interval occurs)  "
          "constraint closed: "
       << closed << endl;
  cout << "        ";
  displayList(list);
  cout << "        values:";
  displayList(values);
  cout << "        occurs:";
  displayList(occurs);
#endif

  VarArray scope;
  getVariables(list, scope);
  VarArray vals;
  getVariables(values, vals);

  solver.add(AllDiff(vals));

  vector<int> lb;
  vector<int> ub;

  for (auto I : occurs) {
    lb.push_back(I.min);
    ub.push_back(I.max);
  }

  if (closed)
    for (auto X : scope) {
      VarArray M;
      for (auto V : vals)
        M.add(Member(X, V));
      solver.add(Sum(M) > 0);
    }

  for (size_t i = 0; i < values.size(); ++i) {
    solver.add(Occurrence(scope, vals[i], lb[i], ub[i]));
  }
}

void XCSP3MistralCallbacks::buildConstraintMinimum(string id,
                                                   vector<XVariable *> &list,
                                                   XCondition &cond) {
#ifdef _VERBOSE_
  cout << "\n    minimum  constraint" << endl;
  cout << "        ";
  displayList(list);
  cout << "        condition: " << cond << endl;
#endif

  VarArray scope;
  getVariables(list, scope);

  Variable MinVar = Min(scope);

  if (cond.operandType == VARIABLE) {
    Variable minv = variable[cond.var];

    if (cond.op == EQ) {
      solver.add(MinVar == minv);
    } else if (cond.op == NE) {
      solver.add(MinVar != minv);
    } else if (cond.op == LE) {
      solver.add(MinVar <= minv);
    } else if (cond.op == LT) {
      solver.add(MinVar < minv);
    } else if (cond.op == GE) {
      for (auto X : scope) {
        solver.add(X >= minv);
      }
    } else if (cond.op == GT) {
      for (auto X : scope) {
        solver.add(X > minv);
      }
    }

  } else if (cond.operandType == INTERVAL) {
    if (cond.op != IN) {
      cout << "s UNSUPPORTED" << _ID_(": IN non-interval") << "\n";
      exit(1);
    }

    Variable MV(cond.min, cond.max);

    solver.add(MinVar == MV);
  } else {

    if (cond.op == EQ) {
      solver.add(MinVar == cond.val);
    } else if (cond.op == NE) {
      Variable cst(cond.val);
      solver.add(MinVar != cst);
    } else if (cond.op == LE) {
      solver.add(MinVar <= cond.val);
    } else if (cond.op == LT) {
      solver.add(MinVar < cond.val);
    } else if (cond.op == GE) {
      for (auto X : scope) {
        solver.add(X >= cond.val);
      }
    } else if (cond.op == GT) {
      for (auto X : scope) {
        solver.add(X > cond.val);
      }
    }
  }

  last_var = MinVar;
}

void XCSP3MistralCallbacks::buildConstraintMinimum(
    string id, vector<XVariable *> &list, XVariable *index, int startIndex,
    RankType rank, XCondition &xc) {
#ifdef _VERBOSE_
  cout << "\n    arg_minimum  constraint" << endl;
  cout << "        ";
  displayList(list);
  cout << "        index:" << *index << endl;
  cout << "        Start index : " << startIndex << endl;
  cout << "        condition: " << xc << endl;
#endif

  VarArray scope;
  getVariables(list, scope);
  buildConstraintMinimum(id, list, xc);
  solver.add(scope[variable[index->id]] == last_var);
}

void XCSP3MistralCallbacks::buildConstraintMinimum(string id,
                                                   vector<Tree *> &trees,
                                                   XCondition &cond) {
#ifdef _VERBOSE_
  cout << "\n    minimum  constraint over trees" << endl;
  cout << "        ";
  displayList(trees);
  cout << "        condition: " << cond << endl;
#endif

  cout << "post trees min\n";

  Vector<Variable> scope;
  for (auto x : trees) {
    scope.add(postExpression(x->root, true).get_var());
  }

  cout << "create min var \n";

  // getVariables(list, scope);

  Variable minVar = Min(scope);

  if (cond.operandType == VARIABLE) {
    Variable minv = variable[cond.var];

    if (cond.op == EQ) {
      solver.add(minVar == minv);
    } else if (cond.op == NE) {
      solver.add(minVar != minv);
    } else if (cond.op == GE) {
      solver.add(minVar >= minv);
    } else if (cond.op == GT) {
      solver.add(minVar > minv);
    } else if (cond.op == LE) {
      for (auto X : scope) {
        solver.add(X <= minv);
      }
    } else if (cond.op == LT) {
      for (auto X : scope) {
        solver.add(X < minv);
      }
    }

  } else if (cond.operandType == INTERVAL) {
    if (cond.op != IN) {
      cout << "s UNSUPPORTED" << _ID_(": IN non-interval") << "\n";
      exit(1);
    }

    Variable MV(cond.min, cond.min);

    solver.add(minVar == MV);
  } else {

    if (cond.op == EQ) {
      solver.add(minVar == cond.val);
    } else if (cond.op == NE) {
      Variable cst(cond.val);
      solver.add(minVar != cst);
    } else if (cond.op == GE) {
      solver.add(minVar >= cond.val);
    } else if (cond.op == GT) {
      solver.add(minVar > cond.val);
    } else if (cond.op == LE) {
      for (auto X : scope) {
        solver.add(X <= cond.val);
      }
    } else if (cond.op == LT) {
      for (auto X : scope) {
        solver.add(X < cond.val);
      }
    }
  }

  std::cout << "end min\n";

  last_var = minVar;
}

void XCSP3MistralCallbacks::buildConstraintMaximum(string id,
                                                   vector<XVariable *> &list,
                                                   XCondition &cond) {
#ifdef _VERBOSE_
  cout << "\n    maximum  constraint" << endl;
  cout << "        ";
  displayList(list);
  cout << "        condition: " << cond << endl;
#endif

  VarArray scope;
  getVariables(list, scope);

  Variable MaxVar = Max(scope);

  if (cond.operandType == VARIABLE) {
    Variable maxv = variable[cond.var];

    if (cond.op == EQ) {
      solver.add(MaxVar == maxv);
    } else if (cond.op == NE) {
      solver.add(MaxVar != maxv);
    } else if (cond.op == GE) {
      solver.add(MaxVar >= maxv);
    } else if (cond.op == GT) {
      solver.add(MaxVar > maxv);
    } else if (cond.op == LE) {
      for (auto X : scope) {
        solver.add(X <= maxv);
      }
    } else if (cond.op == LT) {
      for (auto X : scope) {
        solver.add(X < maxv);
      }
    }

  } else if (cond.operandType == INTERVAL) {
    if (cond.op != IN) {
      cout << "s UNSUPPORTED" << _ID_(": IN non-interval") << "\n";
      exit(1);
    }

    Variable MV(cond.max, cond.max);

    solver.add(MaxVar == MV);
  } else {

    if (cond.op == EQ) {
      solver.add(MaxVar == cond.val);
    } else if (cond.op == NE) {
      Variable cst(cond.val);
      solver.add(MaxVar != cst);
    } else if (cond.op == GE) {
      solver.add(MaxVar >= cond.val);
    } else if (cond.op == GT) {
      solver.add(MaxVar > cond.val);
    } else if (cond.op == LE) {
      for (auto X : scope) {
        solver.add(X <= cond.val);
      }
    } else if (cond.op == LT) {
      for (auto X : scope) {
        solver.add(X < cond.val);
      }
    }
  }

  last_var = MaxVar;
}

void XCSP3MistralCallbacks::buildConstraintMaximum(
    string id, vector<XVariable *> &list, XVariable *index, int startIndex,
    RankType rank, XCondition &xc) {
#ifdef _VERBOSE_
  cout << "\n    arg_maximum  constraint" << endl;
  cout << "        ";
  displayList(list);
  cout << "        index:" << *index << endl;
  cout << "        Start index : " << startIndex << endl;
  cout << "        condition: " << xc << endl;
#endif

  VarArray scope;
  getVariables(list, scope);
  buildConstraintMaximum(id, list, xc);
  solver.add(scope[variable[index->id]] == last_var);
}

void XCSP3MistralCallbacks::buildConstraintMaximum(string id,
                                                   vector<Tree *> &trees,
                                                   XCondition &cond) {
#ifdef _VERBOSE_
  cout << "\n    maximum  constraint" << endl;
  cout << "        ";
  displayList(trees);
  cout << "        condition: " << cond << endl;
#endif

  cout << "post trees min\n";

  Vector<Variable> scope;
  for (auto x : trees) {
    scope.add(postExpression(x->root, true).get_var());
  }

  cout << "create min var \n";

  // getVariables(list, scope);

  Variable MaxVar = Max(scope);

  if (cond.operandType == VARIABLE) {
    Variable maxv = variable[cond.var];

    if (cond.op == EQ) {
      solver.add(MaxVar == maxv);
    } else if (cond.op == NE) {
      solver.add(MaxVar != maxv);
    } else if (cond.op == GE) {
      solver.add(MaxVar >= maxv);
    } else if (cond.op == GT) {
      solver.add(MaxVar > maxv);
    } else if (cond.op == LE) {
      for (auto X : scope) {
        solver.add(X <= maxv);
      }
    } else if (cond.op == LT) {
      for (auto X : scope) {
        solver.add(X < maxv);
      }
    }

  } else if (cond.operandType == INTERVAL) {
    if (cond.op != IN) {
      cout << "s UNSUPPORTED" << _ID_(": IN non-interval") << "\n";
      exit(1);
    }

    Variable MV(cond.max, cond.max);

    solver.add(MaxVar == MV);
  } else {

    if (cond.op == EQ) {
      solver.add(MaxVar == cond.val);
    } else if (cond.op == NE) {
      Variable cst(cond.val);
      solver.add(MaxVar != cst);
    } else if (cond.op == GE) {
      solver.add(MaxVar >= cond.val);
    } else if (cond.op == GT) {
      solver.add(MaxVar > cond.val);
    } else if (cond.op == LE) {
      for (auto X : scope) {
        solver.add(X <= cond.val);
      }
    } else if (cond.op == LT) {
      for (auto X : scope) {
        solver.add(X < cond.val);
      }
    }
  }

  std::cout << "end max\n";

  last_var = MaxVar;
}

void XCSP3MistralCallbacks::buildConstraintElement(string id,
                                                   vector<XVariable *> &list,
                                                   int value) {
#ifdef _VERBOSE_
  cout << "\n    element constant constraint" << endl;
  cout << "        ";
  displayList(list);
  cout << "        value: " << value << endl;
#endif

  VarArray equalities;
  for (auto x : list) {
    ++initial_degree[id_map[x->id]];
    equalities.add(variable[x->id] == value);
  }
  solver.add(BoolSum(equalities) > 0);
}

void XCSP3MistralCallbacks::buildConstraintElement(string id,
                                                   vector<XVariable *> &list,
                                                   XVariable *value) {
#ifdef _VERBOSE_
  cout << "\n    element variable constraint" << endl;
  cout << "        ";
  displayList(list);
  cout << "        value: " << *value << endl;
#endif

  ++initial_degree[id_map[value->id]];
  VarArray equalities;
  for (auto x : list) {
    ++initial_degree[id_map[x->id]];
    equalities.add(variable[x->id] == variable[value->id]);
  }
	
	++initial_degree[id_map[value->id]];
	
  solver.add(BoolSum(equalities) > 0);
}

void XCSP3MistralCallbacks::buildConstraintElement(string id,
                                                   vector<XVariable *> &list,
                                                   int startIndex,
                                                   XVariable *index,
                                                   RankType rank, int value) {
#ifdef _VERBOSE_
  cout << "\n    element constant (with index) constraint" << endl;
  cout << "        ";
  displayList(list);
  cout << "        value: " << value << endl;
  cout << "        Start index : " << startIndex << endl;
  cout << "        index : " << *index << endl;
#endif

  if (rank != ANY) {
    cout << "s UNSUPPORTED" << _ID_(": Element rank != Any") << "\n";
    exit(1);
  }

  // cout << variable[index->id].get_domain() << endl;

  VarArray scope;
  getVariables(list, scope);
	++initial_degree[id_map[index->id]];
  solver.add(Element(scope, variable[index->id], startIndex) == value);
	
	
	
}

void XCSP3MistralCallbacks::buildConstraintElement(
    string id, vector<XVariable *> &list, int startIndex, XVariable *index,
    RankType rank, XVariable *value) {
#ifdef _VERBOSE_
  cout << "\n    element variable (with index) constraint" << endl;
  cout << "        ";
  displayList(list);
  cout << "        value: " << *value << endl;
  cout << "        Start index : " << startIndex << endl;
  cout << "        index : " << *index << endl;
#endif

  if (rank != ANY) {
    cout << "s UNSUPPORTED" << _ID_(": Element rank != Any") << "\n";
    exit(1);
  }

  // for( auto x : list ) {
  // 	std::cout << " " << x->id ;
  // }
  // std::cout << std::endl;
  //
  // // for( auto x : list ) {
  // // 	std::cout << " [" << x->domain->minimum() << "," <<
  // x->domain->maximum() << "]" ;
  // // }
  //
  // for( auto x : list ) {
  // 	std::cout << " " << (int*)(x->domain) ;
  // }
  // std::cout << std::endl;

  VarArray scope;
  getVariables(list, scope);

  // std::cout << scope << " " << variable[index->id] << " " <<
  // variable[value->id] << std::endl;

  // for( auto x : scope ) {
  // 	std::cout << " " << x.get_value() ;
  // }
  // std::cout << std::endl;
	
	++initial_degree[id_map[index->id]];
	++initial_degree[id_map[value->id]];

  solver.add(Element(scope, variable[index->id], startIndex) ==
             variable[value->id]);
}

void XCSP3MistralCallbacks::buildConstraintElement(string id, vector<int> &list,
                                                   int startIndex,
                                                   XVariable *index,
                                                   RankType rank,
                                                   XVariable *value) {
#ifdef _VERBOSE_
  cout << "\n    element variable with list of integers (with index) constraint"
       << endl;
  cout << "        ";
  displayList(list);
  cout << "        value: " << *value << endl;
  cout << "        Start index : " << startIndex << endl;
  cout << "        index : " << *index << endl;
#endif

  VarArray scope;
  for (auto v : list) {
    scope.add(Variable(v, v));
  }
	
	++initial_degree[id_map[index->id]];
	++initial_degree[id_map[value->id]];
	
  solver.add(Element(scope, variable[index->id], startIndex) ==
             variable[value->id]);
}

void XCSP3MistralCallbacks::buildConstraintElement(string id,
                                                   vector<XVariable *> &list,
                                                   XVariable *index,
                                                   int startIndex,
                                                   XCondition &cond) {

#ifdef _VERBOSE_
  cout << "\n    element variable (with index and condition) constraint"
       << endl;
  cout << "        ";
  displayList(list);
  cout << "        Start index : " << startIndex << endl;
  cout << "        index : " << *index << endl;
  cout << "        condittion : " << cond << endl;
#endif

  VarArray scope;
  getVariables(list, scope);
	
	++initial_degree[id_map[index->id]];
	// ++initial_degree[id_map[value->id]];

  Variable Elt = Element(scope, variable[index->id], startIndex);

  if (cond.operandType == VARIABLE) {
    Variable rhs = variable[cond.var];

    if (cond.op == EQ) {
      solver.add(Elt == rhs);
    } else if (cond.op == NE) {
      solver.add(Elt != rhs);
    } else if (cond.op == LE) {
      solver.add(Elt <= rhs);
    } else if (cond.op == LT) {
      solver.add(Elt < rhs);
    } else if (cond.op == GE) {
      solver.add(Elt >= rhs);
    } else if (cond.op == GT) {
      solver.add(Elt > rhs);
    }

  } else if (cond.operandType == INTERVAL) {
    if (cond.op != IN) {
      cout << "s UNSUPPORTED" << _ID_(": IN non-interval") << "\n";
      exit(1);
    }

    solver.add(Elt >= cond.min);
    solver.add(Elt <= cond.max);

  } else {

    if (cond.op == EQ) {
      solver.add(Elt >= cond.min);
      solver.add(Elt <= cond.max);
      // solver.add(Sum(scope, coeffs, cond.val, cond.val));
    } else if (cond.op == NE) {
      Variable cst(cond.val);
      solver.add(Elt != cond.val);
    } else if (cond.op == LE) {
      solver.add(Elt <= cond.val);
    } else if (cond.op == LT) {
      solver.add(Elt < cond.val);
    } else if (cond.op == GE) {
      solver.add(Elt >= cond.val);
    } else if (cond.op == GT) {
      solver.add(Elt > cond.val);
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintElement(
    string id, vector<vector<int>> &matrix, int startRowIndex,
    XVariable *rowIndex, int startColIndex, XVariable *colIndex,
    XVariable *value) {

#ifdef _VERBOSE_
  cout << "\n    element variable (with index) constraint over matrix " << endl;
  for (auto &list : matrix) {
    cout << "        ";
    displayList(list);
  }
  cout << "        value: " << *value << endl;
  cout << "        Start col index : " << startColIndex << endl;
  cout << "        Col index : " << *colIndex << endl;
  cout << "        Start row index : " << startRowIndex << endl;
  cout << "        Row index : " << *rowIndex << endl;
#endif

  VarArray rows;
  for (auto &row : matrix) {
    VarArray scope;
    for (auto v : row) {
      scope.add(Variable(v, v));
    }
    rows.add(Element(scope, variable[rowIndex->id], startRowIndex));
  }
	
	++initial_degree[id_map[rowIndex->id]];
	++initial_degree[id_map[colIndex->id]];
	++initial_degree[id_map[value->id]];
	

  solver.add(Element(rows, variable[colIndex->id], startColIndex) ==
             variable[value->id]);
}

void XCSP3MistralCallbacks::buildConstraintElement(
    string id, vector<vector<XVariable *>> &matrix, int startRowIndex,
    XVariable *rowIndex, int startColIndex, XVariable *colIndex, int value) {

#ifdef _VERBOSE_
  cout << "\n    element variable (with index) constraint over matrix " << endl;
  for (auto &list : matrix) {
    cout << "        ";
    displayList(list);
  }
  cout << "        value: " << value << endl;
  cout << "        Start col index : " << startColIndex << endl;
  cout << "        Col index : " << *colIndex << endl;
  cout << "        Start row index : " << startRowIndex << endl;
  cout << "        Row index : " << *rowIndex << endl;
#endif

  VarArray rows;
  for (auto &row : matrix) {
    VarArray scope;
    getVariables(row, scope);
    rows.add(Element(scope, variable[rowIndex->id], startRowIndex));
  }
	
	++initial_degree[id_map[rowIndex->id]];
	++initial_degree[id_map[colIndex->id]];
	// ++initial_degree[id_map[value->id]];

  solver.add(Element(rows, variable[colIndex->id], startColIndex) == value);
}

void XCSP3MistralCallbacks::buildConstraintElement(
    string id, vector<vector<XVariable *>> &matrix, int startRowIndex,
    XVariable *rowIndex, int startColIndex, XVariable *colIndex,
    XVariable *value) {

#ifdef _VERBOSE_
  cout << "\n    element variable (with index) constraint over matrix " << endl;
  for (auto &list : matrix) {
    cout << "        ";
    displayList(list);
  }
  cout << "        value: " << *value << endl;
  cout << "        Start col index : " << startColIndex << endl;
  cout << "        Col index : " << *colIndex << endl;
  cout << "        Start row index : " << startRowIndex << endl;
  cout << "        Row index : " << *rowIndex << endl;
#endif

  VarArray rows;
  for (auto &row : matrix) {
    VarArray scope;
    getVariables(row, scope);
    rows.add(Element(scope, variable[rowIndex->id], startRowIndex));
  }
	
	++initial_degree[id_map[rowIndex->id]];
	++initial_degree[id_map[colIndex->id]];
	++initial_degree[id_map[value->id]];

  solver.add(Element(rows, variable[colIndex->id], startColIndex) ==
             variable[value->id]);
}

void XCSP3MistralCallbacks::buildConstraintChannel(string id,
                                                   vector<XVariable *> &list,
                                                   int startIndex) {
#ifdef _VERBOSE_
  cout << "\n    channel constraint" << endl;
  cout << "        ";
  displayList(list);
  cout << "        Start index : " << startIndex << endl;
#endif

  VarArray scope;
  getVariables(list, scope);
  solver.add(AllDiff(scope));
  for (size_t i = 0; i < scope.size; ++i) {
    for (size_t j = i + 1; j < scope.size; ++j) {
      solver.add((scope[i] == j + startIndex) == (scope[j] == i + startIndex));
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintChannel(string id,
                                                   vector<XVariable *> &list1,
                                                   int startIndex1,
                                                   vector<XVariable *> &list2,
                                                   int startIndex2) {
#ifdef _VERBOSE_
  cout << "\n    channel constraint" << endl;
  cout << "        list1 ";
  displayList(list1);
  cout << "        list2 ";
  displayList(list2);
#endif

  // assert(list1.size() == list2.size());

  VarArray scope[2];
  int sindex[2];
  if (list1.size() == list2.size()) {
    getVariables(list1, scope[0]);
    getVariables(list2, scope[1]);
    sindex[0] = startIndex1;
    sindex[1] = startIndex2;

    solver.add(AllDiff(scope[0]));
    solver.add(AllDiff(scope[1]));
    for (size_t i = 0; i < scope[0].size; ++i) {
      for (size_t j = 0; j < scope[1].size; ++j) {
        solver.add((scope[0][i] == j + sindex[0]) ==
                   (scope[1][j] == i + sindex[1]));
      }
    }
    return;
  } else if (list1.size() < list2.size()) {
    getVariables(list1, scope[0]);
    getVariables(list2, scope[1]);
    sindex[0] = startIndex1;
    sindex[1] = startIndex2;
  } else {
    getVariables(list1, scope[1]);
    getVariables(list2, scope[0]);
    sindex[1] = startIndex1;
    sindex[0] = startIndex2;
  }

  for (size_t i = 0; i < scope[0].size; ++i) {
    for (size_t j = 0; j < scope[1].size; ++j) {
      solver.add((scope[0][i] == j + sindex[0]) <=
                 (scope[1][j] == i + sindex[1]));
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintChannel(string id,
                                                   vector<XVariable *> &list,
                                                   int startIndex,
                                                   XVariable *value) {
#ifdef _VERBOSE_
  cout << "\n    channel constraint" << endl;
  cout << "        ";
  displayList(list);
  cout << "        value: " << *value << endl;
#endif

  VarArray scope;
  getVariables(list, scope);

  // solver.add( Occurrence(scope, 1, 1, 1) );
  solver.add(Sum(scope) == 1); // must be Boolean variables
  solver.add(Element(scope, variable[value->id], startIndex) == 1);
}

void XCSP3MistralCallbacks::buildConstraintStretch(string id,
                                                   vector<XVariable *> &list,
                                                   vector<int> &values,
                                                   vector<XInterval> &widths) {
#ifdef _VERBOSE_
  cout << "\n    stretch constraint" << endl;
  cout << "        ";
  displayList(list);
  cout << "        values :";
  displayList(values);
  cout << "        widths:";
  displayList(widths);
#endif

  VarArray scope;
  getVariables(list, scope);
  std::vector<int> lb;
  std::vector<int> ub;

  for (auto I : widths) {
    lb.push_back(I.min);
    ub.push_back(I.max);
  }

  solver.add(Stretch(scope, values, lb, ub));
}

void XCSP3MistralCallbacks::buildConstraintStretch(
    string id, vector<XVariable *> &list, vector<int> &values,
    vector<XInterval> &widths, vector<vector<int>> &patterns) {
#ifdef _VERBOSE_
  cout << "\n    stretch constraint (with patterns)" << endl;
  cout << "        ";
  displayList(list);
  cout << "        values :";
  displayList(values);
  cout << "        widths:";
  displayList(widths);
  cout << "        patterns";
  for (unsigned int i = 0; i < patterns.size(); i++)
    cout << "(" << patterns[i][0] << "," << patterns[i][1] << ") ";
  cout << endl;
#endif

  VarArray scope;
  getVariables(list, scope);
  std::vector<int> lb;
  std::vector<int> ub;
  std::vector<int> trans;

  for (auto I : widths) {
    lb.push_back(I.min);
    ub.push_back(I.max);
  }

  for (unsigned int i = 0; i < patterns.size(); i++) {
    trans.push_back(patterns[i][0]);
    trans.push_back(patterns[i][1]);
  }

  solver.add(Stretch(scope, values, lb, ub, trans));
}

void XCSP3MistralCallbacks::buildConstraintNoOverlap(
    string id, vector<XVariable *> &origins, vector<int> &lengths,
    bool zeroIgnored) {
#ifdef _VERBOSE_
  cout << "\n    nooverlap constraint" << endl;
  cout << "        origins";
  displayList(origins);
  cout << "        lengths";
  displayList(lengths);
#endif

  VarArray tasks;
  getVariables(origins, tasks);

  for (size_t i = 0; i < tasks.size; ++i) {
    for (size_t j = i + 1; j < tasks.size; ++j) {
      solver.add(
          Free(ReifiedDisjunctive(tasks[i], tasks[j], lengths[i], lengths[j])));
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintNoOverlap(
    string id, vector<XVariable *> &origins, vector<XVariable *> &lengths,
    bool zeroIgnored) {
#ifdef _VERBOSE_
  cout << "\n    nooverlap constraint" << endl;
  cout << "        origins:";
  displayList(origins);
  cout << "        lengths";
  displayList(lengths);
#endif

  VarArray tasks;
  getVariables(origins, tasks);
  VarArray durations;
  getVariables(lengths, durations);

  for (size_t i = 0; i < tasks.size; ++i) {
    for (size_t j = i + 1; j < tasks.size; ++j) {
      solver.add((tasks[i] + durations[i] <= tasks[j]) ||
                 (tasks[j] + durations[j] <= tasks[i]));
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintNoOverlap(
    string id, vector<vector<XVariable *>> &origins,
    vector<vector<int>> &lengths, bool zeroIgnored) {
#ifdef _VERBOSE_
  cout << "\n    kdim (int lengths) nooverlap constraint" << endl;
  cout << "origins: " << endl;
  for (unsigned int i = 0; i < origins.size(); i++) {
    cout << "        ";
    displayList(origins[i]);
  }
  cout << "lengths: " << endl;
  for (unsigned int i = 0; i < origins.size(); i++) {
    cout << "        ";
    displayList(lengths[i]);
  }
#endif

  int num_boxes = origins.size();
  int num_dim = begin(origins)->size();

  for (int i = 0; i < num_boxes; ++i) {
    for (int j = i + 1; j < num_boxes; ++j) {
      VarArray distinct;
      for (int k = 0; k < num_dim; ++k) {
        ++initial_degree[id_map[origins[i][k]->id]];
        ++initial_degree[id_map[origins[j][k]->id]];
        Variable Xik = variable[origins[i][k]->id];
        Variable Xjk = variable[origins[j][k]->id];
        distinct.add((Xik + lengths[i][k] <= Xjk) ||
                     (Xjk + lengths[j][k] <= Xik));
      }
      solver.add(BoolSum(distinct, 1, distinct.size));
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintNoOverlap(
    string id, vector<vector<XVariable *>> &origins,
    vector<vector<XVariable *>> &lengths, bool zeroIgnored) {
#ifdef _VERBOSE_
  cout << "\n    kdim (lenghts vars nooverlap constraint" << endl;
  cout << "origins: " << endl;
  for (unsigned int i = 0; i < origins.size(); i++) {
    cout << "        ";
    displayList(origins[i]);
  }
  cout << "lengths: " << endl;
  for (unsigned int i = 0; i < origins.size(); i++) {
    cout << "        ";
    displayList(lengths[i]);
  }
#endif

  int num_boxes = origins.size();
  int num_dim = begin(origins)->size();

  for (int i = 0; i < num_boxes; ++i) {
    for (int j = i + 1; j < num_boxes; ++j) {
      VarArray distinct;
      for (int k = 0; k < num_dim; ++k) {
        ++initial_degree[id_map[origins[i][k]->id]];
        ++initial_degree[id_map[origins[j][k]->id]];
        ++initial_degree[id_map[lengths[i][k]->id]];
        ++initial_degree[id_map[lengths[j][k]->id]];

        Variable Xik = variable[origins[i][k]->id];
        Variable Xjk = variable[origins[j][k]->id];
        Variable Pik = variable[lengths[i][k]->id];
        Variable Pjk = variable[lengths[j][k]->id];
        distinct.add((Xik + Pik <= Xjk) || (Xjk + Pjk <= Xik));
      }
      solver.add(BoolSum(distinct, 1, distinct.size));
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintInstantiation(
    string id, vector<XVariable *> &list, vector<int> &values) {
#ifdef _VERBOSE_
  cout << "\n    instantiation constraint" << endl;
  cout << "        list:";
  displayList(list);
  cout << "        values:";
  displayList(values);
#endif

  for (size_t i = 0; i < list.size(); ++i) {
    ++initial_degree[id_map[list[i]->id]];
    solver.add(variable[list[i]->id] == values[i]);
  }
}

void XCSP3MistralCallbacks::p_encode_cumulative(
    Mistral::Solver &s, Mistral::Vector<Mistral::Variable> &start,
    Mistral::Vector<Mistral::Variable> &dur,
    Mistral::Vector<Mistral::Variable> &req, Mistral::Variable cap) {

  int orig = INFTY, horizon = 0, ddate, rdate, unit_req = true, unit_dur = true,
      disjunctive = false, same_req = true, same_dur = true,
      var_dur = false; //, var_req=false;
  //

  for (size_t i = 0; i < start.size; ++i) {

    // std::cout << start[i].get_domain() << std::endl;

    ddate = start[i].get_max() + dur[i].get_max();
    if (horizon < ddate) {
      horizon = ddate;
    }
    rdate = start[i].get_min();
    if (orig > rdate) {
      orig = rdate;
    }
    // if(!req[i].is_ground())
    // 	var_req = true;
    if (!dur[i].is_ground())
      var_dur = true;
    if (unit_req and (req[i].get_max() > 1 || req[i].get_min() < 1))
      unit_req = false;
    if (unit_dur and (dur[i].get_max() > 1 || dur[i].get_min() < 1))
      unit_dur = false;
    if (same_dur && (!dur[i].is_ground() ||
                     (i && dur[i - 1].get_min() != dur[i].get_min())))
      same_dur = false;
    if (same_req && (!req[i].is_ground() ||
                     (i && req[i - 1].get_min() != req[i].get_min())))
      same_req = false;

    // if(!dur[i].is_ground() || !req[i].is_ground()) {
    //   simple = false;
    // }
  }

#ifdef _DEBUG_CUMULATIVE
  std::cout << "\nhorizon  = " << horizon << std::endl
            << "capacity = " << cap << std::endl
            << "tasks    = " << start.size << std::endl;
#endif

  // int *order = new int[start.size];
  // for(int i=0; i<start.size; ++i) {
  // 	order[i] = i;
  // }
  std::vector<int> order;
  for (size_t i = 0; i < start.size; ++i) {
    order.push_back(i);
  }

  _demand_.copy(req);

  std::sort(begin(order), end(order), [&](int x, int y) {
    return _demand_[x].get_min() < _demand_[y].get_min();
  });
  _demand_.neutralise();

  if (disjunctive and
      (start.size > 1 &&
       req[order[0]].get_min() + req[order[1]].get_min() > cap.get_max())) {

#ifdef _DEBUG_CUMULATIVE
    std::cout << "min (" << req[order[0]].get_min() << ") and next min ("
              << req[order[1]].get_min()
              << ") demands are incompatible => disjunctive" << std::endl;
#endif

    disjunctive = true;
  }

#ifdef _DEBUG_CUMULATIVE
  for (size_t i = 0; i < start.size; ++i) {
    int oi = order[i];
    std::cout << start[oi] << "           [" << start[oi].get_min() << ".."
              << start[oi].get_max() + dur[oi].get_max() << "]";
    if (dur[oi].is_ground())
      std::cout << " p=" << dur[oi].get_max();
    else
      std::cout << " p=" << dur[oi].get_domain();
    if (req[oi].is_ground())
      std::cout << " r=" << req[oi].get_max();
    else
      std::cout << " r=" << req[oi].get_domain();
    std::cout << std::endl;
  }
// exit(1);
#endif

  // delete [] order ;

  int size_discretization = horizon * start.size;
  int size_flow = start.size * start.size * cap.get_max();

#ifdef _DEBUG_CUMULATIVE
  if (unit_dur && disjunctive) {

    std::cout << "this is an alldiff, but do nothing about it!" << std::endl;

    // s.add(AllDiff(start));
  } else if (unit_dur && unit_req) {

    std::cout << "this is a gcc, but do nothing about it!" << std::endl;

    // int *lb = new int[horizon-orig];
    // int *ub = new int[horizon-orig];
    // for(int i=orig; i<horizon; ++i) {
    // 	lb[i-orig] = 0; //cap.get_min();
    // 	ub[i-orig] = cap.get_max();
    // }
    // s.add(Occurrences(start, orig, horizon-1, lb, ub));
    // delete [] lb;
    // delete [] ub;

  } else if (disjunctive && !var_dur) {

    if (same_dur) {
      std::cout << "this is an interdistance, but do nothing about it!"
                << std::endl;
    } else {
      std::cout << "this is a disjunctive, but do nothing about it!"
                << std::endl;
    }

    // p_disjunctive(s, start, dur);
  }
#endif
// else {

#ifdef _DEBUG_CUMULATIVE
  if (same_dur && same_req) {
    std::cout << "this is a symmetric cumulative" << std::endl;
  } else if (same_dur) {
    std::cout << "this is a cumulative with equal processing times"
              << std::endl;
  } else if (same_req) {
    std::cout << "this is a cumulative with equal demands" << std::endl;
  } else {
    std::cout << "this is a general cumulative" << std::endl;
  }

  std::cout << "size of the time-discretization encoding = "
            << size_discretization << std::endl;
  std::cout << "size of the flow encoding = " << size_flow << std::endl;
#endif

  if (false && (size_flow < size_discretization || var_dur)) {
    p_cumulative_flow(solver, start, dur, req, cap, horizon);
  } else {
    p_cumulative_discretization(solver, start, dur, req, cap, horizon);
  }

  // }

  // exit(1);
}

void XCSP3MistralCallbacks::p_cumulative_flow(Solver &s,
                                              Vector<Variable> &start,
                                              Vector<Variable> &dur,
                                              Vector<Variable> &req,
                                              Variable cap, int horizon) {

#ifdef _DEBUG_CUMULATIVE
  std::cout << "CAPACITY = " << cap.get_domain() << std::endl;
#endif

  int n = start.size;
  Vector<Variable> end;

  for (int i = 0; i < n; ++i) {
    end.add(dur[i].is_ground() ? start[i] + dur[i].get_max()
                               : start[i] + dur[i]);

#ifdef _DEBUG_CUMULATIVE
    std::cout << req[i].get_domain() << std::endl;
#endif
  }

  Variable **flow = new Variable *[n + 1];
  for (int i = 0; i <= n; ++i) {
    flow[i] = new Variable[n + 1];
  }

  for (int i = 0; i < n; ++i) {

    flow[i][n] = Variable(0, req[i].get_max());
    flow[n][i] = Variable(0, req[i].get_max());

    for (int j = i + 1; j < n; ++j) {
      // for each pair of tasks, create a variable for the amount flow going
      // from ti to tj
      int mf = (req[i].get_max() > req[j].get_max() ? req[j].get_max()
                                                    : req[i].get_max());
      flow[i][j] = Variable(0, mf);
      flow[j][i] = Variable(0, mf);

      if (mf > 0) {
        // the flow can go only in one direction
        s.add((flow[i][j] <= 0) || (flow[j][i] <= 0));
      }
    }
  }

  // flow conservation constraints
  Vector<Variable> outflow;
  Vector<Variable> inflow;
  for (int i = 0; i < n; ++i) {

#ifdef _DEBUG_CUMULATIVE
    std::cout << " conservation in node " << start[i].id() << std::endl;
#endif

    for (int j = 0; j <= n; ++j) {
      if (i != j) {

#ifdef _DEBUG_CUMULATIVE
        std::cout << i << "." << j << std::endl;
#endif

        if (flow[i][j].get_max() > 0) {

#ifdef _DEBUG_CUMULATIVE
          std::cout << "  + flow to/from " << (j < n ? start[j].id() : n)
                    << std::endl;
#endif

          outflow.add(flow[i][j]);
          inflow.add(flow[j][i]);
        }
      }
    }

#ifdef _DEBUG_CUMULATIVE
    std::cout << "Must be equal to: ";
    std::cout.flush();
#endif

    if (req[i].is_ground()) {

#ifdef _DEBUG_CUMULATIVE
      std::cout << "[" << req[i].get_min() << ", " << req[i].get_max() << "]"
                << std::endl;
#endif

      Variable sum_exp = Sum(outflow, req[i].get_min(), req[i].get_max());

      s.add(sum_exp);

      s.add(Sum(outflow, req[i].get_min(), req[i].get_max()));
      s.add(Sum(inflow, req[i].get_min(), req[i].get_max()));

    } else {

#ifdef _DEBUG_CUMULATIVE
      std::cout << req[i].get_domain() << std::endl;
#endif

      s.add(Sum(outflow) == req[i]);
      s.add(Sum(inflow) == req[i]);
    }
    outflow.clear();
    inflow.clear();
  }

#ifdef _DEBUG_CUMULATIVE
  std::cout << " post total flow constraints:" << std::endl;
#endif

  for (int j = 0; j < n; ++j) {
    if (flow[n][j].get_max() > 0)
      outflow.add(flow[n][j]);
    if (flow[j][n].get_max() > 0)
      inflow.add(flow[j][n]);
  }

  if (cap.is_ground()) {

#ifdef _DEBUG_CUMULATIVE
    std::cout << "fixed capacity" << std::endl;
#endif

    s.add(Sum(outflow, cap.get_max(), cap.get_max()));
    s.add(Sum(inflow, cap.get_max(), cap.get_max()));
  } else {

#ifdef _DEBUG_CUMULATIVE
    std::cout << "variable capacity" << std::endl;
#endif

    s.add(Sum(outflow) == cap);
    s.add(Sum(inflow) == cap);
  }

  // There can be a flow going from i to j only if i precedes j
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i != j) {
        if (flow[i][j].get_max() > 0) {
          s.add((flow[i][j] > 0) <= (end[i] <= start[j]));
          s.add((flow[j][i] > 0) <= (end[j] <= start[i]));
        }
      }
    }
  }

  s.consolidate();

#ifdef _DEBUG_CUMULATIVE
  std::cout << s << std::endl;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i != j) {
        std::cout << " " << std::setw(3) << flow[i][j].id();
      } else {
        std::cout << "    ";
      }
    }
    std::cout << std::endl;
  }
#endif

  for (int i = 0; i <= n; ++i) {
    delete[] flow[i];
  }
  delete[] flow;
}

void XCSP3MistralCallbacks::p_cumulative_discretization(
    Solver &s, Vector<Variable> &start, Vector<Variable> &dur,
    Vector<Variable> &req, Variable cap, int horizon) {

  Variable in_process;
  Vector<Variable> scope;
  Vector<int> demand;
  int process_time;
  int total;
  int min_req;
  int isBoolean;

  // post a sum constraint for each time point
  for (int t = 0; t <= horizon; ++t) {
    demand.clear();
    scope.clear();
    total = 0;
    isBoolean = true;
    min_req = Mistral::INFTY;
    for (size_t i = 0; i < start.size; ++i) {
      if (start[i].get_min() <= t &&
          start[i].get_max() + dur[i].get_max() >= t) {
        total += req[i].get_max();
        if (min_req > req[i].get_min()) {
          min_req = req[i].get_min();
        }

        // this tasks can be in process at time t
        if (dur[i].is_ground()) {
          process_time = dur[i].get_min();
          // constant duration
          if (process_time == 1)
            in_process = (start[i] == t);
          else
            in_process = (Member(start[i], t - process_time + 1, t));
        } else {
          // TODO: tasks with variable durations
          cout << "s UNSUPPORTED"
               << _ID_(": Discretization of Cumulative with variable durations")
               << "\n";
          exit(1);
        }
        if (req[i].is_ground()) {
          // std::cout << "the demand is fixed: " << req[i].get_min() <<
          // std::endl;
          scope.add(in_process);
          demand.add(req[i].get_min());
        } else {
          // std::cout << "the demand is variable: " << req[i].get_domain() <<
          // std::endl;
          scope.add((in_process * req[i]));
          demand.add(1);
          if (req[i].get_max() > 1) {
            isBoolean = false;
          }
        }
      } else {
        // std::cout << "  -task t" << i << " cannot be in process\n";
      }
    }

    if (scope.size > 1 && total > cap.get_min()) {
      if (cap.is_ground()) {
        if (isBoolean) {
          // check if this is in fact a clause
          if (2 * min_req > cap.get_max()) {
            s.add(BoolSum(scope, -Mistral::INFTY, 1));
          } else {
            s.add(BoolSum(scope, demand, -Mistral::INFTY, cap.get_min()));
          }
        } else
          s.add(Sum(scope, demand, -Mistral::INFTY, cap.get_min()));
      } else {
        scope.add(cap);
        demand.add(-1);
        s.add(Sum(scope, demand, -Mistral::INFTY, 0));
      }
    }
  }
}

void XCSP3MistralCallbacks::buildConstraintCumulative(
    string id, vector<XVariable *> &origins, vector<int> &lengths,
    vector<int> &heights, XCondition &cond) {

#ifdef _VERBOSE_
  cout << "\n    Cumulative constraint (int lengths, int heights) " << cond
       << endl;
  cout << "        ";
  displayList(origins);
  cout << "        lengths:";
  displayList(lengths);
  cout << "        heights:";
  displayList(heights);
#endif

#ifdef _DEBUG_CUMULATIVE
  std::cout << "Decompose cumulative constraint" << std::endl;
#endif

  if (cond.op != LE) {
    cout << "s UNSUPPORTED" << _ID_(": Cumulative with lower bound capacity")
         << "\n";
    exit(1);
  }

  Variable cap;
  VarArray start;
  VarArray dur;
  VarArray req;

  getVariables(origins, start);
  for (auto l : lengths) {
    dur.add(Variable(l, l));
  }
  for (auto h : heights) {
    req.add(Variable(h, h));
  }
  if (cond.operandType == VARIABLE) {
    cap = variable[cond.var];
  } else if (cond.operandType == INTEGER) {
    cap = Variable(cond.val, cond.val);
  } else {
    cout << "s UNSUPPORTED" << _ID_(": Weird Cumulative") << "\n";
    exit(1);
  }

  p_encode_cumulative(solver, start, dur, req, cap);
}

void XCSP3MistralCallbacks::buildConstraintCumulative(
    string id, vector<XVariable *> &origins, vector<int> &lengths,
    vector<XVariable *> &heights, XCondition &cond) {

#ifdef _VERBOSE_
  cout << "\n    Cumulative constraint (int lengths, int heights) " << cond
       << endl;
  cout << "        ";
  displayList(origins);
  cout << "        lengths:";
  displayList(lengths);
  cout << "        heights:";
  displayList(heights);
#endif

#ifdef _DEBUG_CUMULATIVE
  std::cout << "Decompose cumulative constraint" << std::endl;
#endif

  if (cond.op != LE) {
    cout << "s UNSUPPORTED" << _ID_(": Cumulative with lower bound capacity")
         << "\n";
    exit(1);
  }

  Variable cap;
  VarArray start;
  VarArray dur;
  VarArray req;

  getVariables(origins, start);
  // getVariables(lengths, dur);
  getVariables(heights, req);
  for (auto l : lengths) {
    dur.add(Variable(l, l));
  }
  // for (auto h : heights) {
  //   req.add(Variable(h, h));
  // }
  if (cond.operandType == VARIABLE) {
    cap = variable[cond.var];
  } else if (cond.operandType == INTEGER) {
    cap = Variable(cond.val, cond.val);
  } else {
    cout << "s UNSUPPORTED" << _ID_(": Weird Cumulative") << "\n";
    exit(1);
  }

  p_encode_cumulative(solver, start, dur, req, cap);
}

void XCSP3MistralCallbacks::buildConstraintCumulative(
    string id, vector<XVariable *> &origins, vector<XVariable *> &lengths,
    vector<int> &heights, XCondition &cond) {
#ifdef _VERBOSE_
  cout << "\n    Cumulative constraint (int lengths, int heights) " << cond
       << endl;
  cout << "        ";
  displayList(origins);
  cout << "        lengths:";
  displayList(lengths);
  cout << "        heights:";
  displayList(heights);
#endif

#ifdef _DEBUG_CUMULATIVE
  std::cout << "Decompose cumulative constraint" << std::endl;
#endif

  if (cond.op != LE) {
    cout << "s UNSUPPORTED" << _ID_(": Cumulative with lower bound capacity")
         << "\n";
    exit(1);
  }

  Variable cap;
  VarArray start;
  VarArray dur;
  VarArray req;

  getVariables(origins, start);
  getVariables(lengths, dur);
  // getVariables(heights, req);
  // for (auto l : lengths) {
  //   dur.add(Variable(l, l));
  // }
  for (auto h : heights) {
    req.add(Variable(h, h));
  }
  if (cond.operandType == VARIABLE) {
    cap = variable[cond.var];
  } else if (cond.operandType == INTEGER) {
    cap = Variable(cond.val, cond.val);
  } else {
    cout << "s UNSUPPORTED" << _ID_(": Weird Cumulative") << "\n";
    exit(1);
  }

  p_encode_cumulative(solver, start, dur, req, cap);
}

void XCSP3MistralCallbacks::buildConstraintCumulative(
    string id, vector<XVariable *> &origins, vector<XVariable *> &lengths,
    vector<XVariable *> &heights, XCondition &cond) {
#ifdef _VERBOSE_
  cout << "\n    Cumulative constraint (int lengths, int heights) " << cond
       << endl;
  cout << "        ";
  displayList(origins);
  cout << "        lengths:";
  displayList(lengths);
  cout << "        heights:";
  displayList(heights);
#endif

#ifdef _DEBUG_CUMULATIVE
  std::cout << "Decompose cumulative constraint" << std::endl;
#endif

  if (cond.op != LE) {
    cout << "s UNSUPPORTED" << _ID_(": Cumulative with lower bound capacity")
         << "\n";
    exit(1);
  }

  Variable cap;
  VarArray start;
  VarArray dur;
  VarArray req;

  getVariables(origins, start);
  getVariables(lengths, dur);
  getVariables(heights, req);
  // for (auto l : lengths) {
  //   dur.add(Variable(l, l));
  // }
  // for (auto h : heights) {
  //   req.add(Variable(h, h));
  // }
  if (cond.operandType == VARIABLE) {
    cap = variable[cond.var];
  } else if (cond.operandType == INTEGER) {
    cap = Variable(cond.val, cond.val);
  } else {
    cout << "s UNSUPPORTED" << _ID_(": Weird Cumulative") << "\n";
    exit(1);
  }

  p_encode_cumulative(solver, start, dur, req, cap);
}

void XCSP3MistralCallbacks::buildConstraintCircuit(string id,
                                                   vector<XVariable *> &list,
                                                   int startIndex) {

#ifdef _VERBOSE_
  cout << "\n    circuit constraint" << endl;
  cout << "        list:";
  displayList(list);
  cout << "        startIndex:" << startIndex << endl;
#endif
	
	if(startIndex != 0)
  throw runtime_error(
      "s non-null start inde in circuit not supported");

  VarArray next;
  getVariables(list, next);
	
	
	int n{static_cast<int>(list.size()) - 1};
  //
  VarArray timestamp(list.size(), 0, n);

  solver.add(AllDiff(timestamp));
  solver.add(AllDiff(next));
	
	// solver.add(timestamp[startIndex] == 0);
	// for(int i=0; i<n; ++i) {
	// 	// int j{(startIndex+i) % n};
	// 	solver.add((next[i - startIndex] != i) <= (timestamp[next[i - startIndex]] == (timestamp[i - startIndex] + 1)));
	// 	// solver.add((next[0] == i) <= (timestamp[i] == 0));
	// }
	// solver.add(timestamp[next[n] - startIndex] == 0);
	
	for(int i=0; i<n; ++i) {
		// int j{(startIndex+i) % n};
		solver.add((next[i] != i) <= (timestamp[next[i]] == (timestamp[i] + 1)));
		// solver.add((next[0] == i) <= (timestamp[i] == 0));
	}
	solver.add(timestamp[next[n]] == 0);
	
	// solver.add(next[(startIndex+n)] == 0);
	//
	//
	//
	//   for (size_t i = 0; i < timestamp.size; ++i) {
	//     solver.add((timestamp[i] == 0) <= (next[i] != i));
	//     solver.add((next[i] != i) <= ((timestamp[next[i]] == 0) ||
	//                                   (timestamp[next[i]] == timestamp[i] + 1)));
	//   }


	// pos[next[0]] == pos[0] + 1
	//
	// pos[next[n-1]] = 0;
	
	// next[0] 1
	// next[1] 2
	// next[2] 3
	// next[3] 0
	//
	// 	pos[0]
	// 	pos[1]
	// 	pos[2]
	// 	pos[3]
		

}

void XCSP3MistralCallbacks::buildConstraintCircuit(string id,
                                                   vector<XVariable *> &list,
                                                   int startIndex, int size) {
  cout << "s UNSUPPORTED" << _ID_(": Fixed-size Circuit") << "\n";
  exit(1);
}

void XCSP3MistralCallbacks::buildConstraintCircuit(string id,
                                                   vector<XVariable *> &list,
                                                   int startIndex,
                                                   XVariable *size) {
  cout << "s UNSUPPORTED" << _ID_(": Circuit?") << "\n";
  exit(1);
}

void XCSP3MistralCallbacks::buildObjectiveMinimizeExpression(string expr) {
#ifdef _VERBOSE_
  cout << "\n    objective: minimize str expr " << expr << endl;
#endif

#ifdef _OLD_TREE
  XCSP3Mistral::Tree tree(expr);
  goal =
      new Goal(Goal::MINIMIZATION, postExpression(tree.root, true).get_var());
  tree.dispose();
#else
  Tree tree(expr);

  Variable obj = postExpression(tree.root, true).get_var();

  solver.add(Free(obj));

  // cout << "domain: " << obj.get_domain() << endl;

  goal = new Goal(Goal::MINIMIZATION, obj);

#endif
}

void XCSP3MistralCallbacks::buildObjectiveMaximizeExpression(string expr) {
#ifdef _VERBOSE_
  cout << "\n    objective: maximize str expr" << expr << endl;
#endif

#ifdef _OLD_TREE
  XCSP3Mistral::Tree tree(expr);
  goal =
      new Goal(Goal::MAXIMIZATION, postExpression(tree.root, true).get_var());
  tree.dispose();
#else
  Tree tree(expr);

  Variable obj = postExpression(tree.root, true).get_var();

  solver.add(Free(obj));

  goal = new Goal(Goal::MAXIMIZATION, obj);
// delete tree.root;
#endif
}

/**
 * The callback function related to an objective minimize a sum/product
 * See http://xcsp.org/specifications/objectives
 *
 * Example:
 * <objectives>
 *   <minimize type="sum">
 *     <list> ne(s[2],2) ne(s[3],3) ... </list>
 *     <coeffs> 2 4 ... </coeffs>
 *   </minimize>
 * <objectives>
 *
 * @param type SUM, PRODUCT...
 * @param trees the expressions
 * @param coefs the vector of coefficients
 */
void XCSP3MistralCallbacks::buildObjectiveMinimize(ExpressionObjective type,
                                                   vector<Tree *> &trees,
                                                   vector<int> &coefs) {

#ifdef _VERBOSE_
  cout << "\n    objective: minimize sum of tree expr" << endl;
#endif

  Vector<Variable> scope;
  for (auto x : trees) {
    scope.add(postExpression(x->root, true).get_var());
  }

  Variable objective;
  switch (type) {
  case SUM_O:
    objective = Sum(scope, coefs).get_var();
    break;
  case MAXIMUM_O:
    for (auto c : coefs) {
      if (c != 1)
        throw runtime_error(
            "Max operand in objective with coefs is not supported");
    }
    objective = Max(scope).get_var();
    break;
  case MINIMUM_O:
    for (auto c : coefs) {
      if (c != 1)
        throw runtime_error(
            "Min operand in objective with coefs is not supported");
    }
    objective = Min(scope).get_var();
    break;
  case EXPRESSION_O:
    throw runtime_error("Expression operand in objective not supported");
  case NVALUES_O:
    throw runtime_error("NValue operand in objective not supported");
  case LEX_O:
    throw runtime_error("Lex operand in objective not supported");
  case PRODUCT_O:
    throw runtime_error("Product operand in objective not supported");
  }

  solver.add(Free(objective));

  goal = new Goal(Goal::MINIMIZATION, objective.get_var());
}

void XCSP3MistralCallbacks::buildObjectiveMaximize(ExpressionObjective type,
                                                   vector<Tree *> &trees,
                                                   vector<int> &coefs) {

#ifdef _VERBOSE_
  cout << "\n    objective: maximize sum of tree expr" << endl;
#endif

  Vector<Variable> scope;
  for (auto x : trees) {
    scope.add(postExpression(x->root, true).get_var());
  }

  Variable objective;
  switch (type) {
  case SUM_O:
    objective = Sum(scope, coefs).get_var();
    break;
  case MAXIMUM_O:
    for (auto c : coefs) {
      if (c != 1)
        throw runtime_error(
            "Max operand in objective with coefs is not supported");
    }
    objective = Max(scope).get_var();
    break;
  case MINIMUM_O:
    for (auto c : coefs) {
      if (c != 1)
        throw runtime_error(
            "Min operand in objective with coefs is not supported");
    }
    objective = Min(scope).get_var();
    break;
  case EXPRESSION_O:
    throw runtime_error("Expression operand in objective not supported");
  case NVALUES_O:
    throw runtime_error("NValue operand in objective not supported");
  case LEX_O:
    throw runtime_error("Lex operand in objective not supported");
  case PRODUCT_O:
    throw runtime_error("Product operand in objective not supported");
  }

  solver.add(Free(objective));

  goal = new Goal(Goal::MAXIMIZATION, objective.get_var());
}

void XCSP3MistralCallbacks::buildObjectiveMinimizeVariable(XVariable *x) {
#ifdef _VERBOSE_
  cout << "\n    objective: minimize variable " << x->id << endl;
#endif

  goal = new Goal(Goal::MINIMIZATION, variable[x->id]);
}

void XCSP3MistralCallbacks::buildObjectiveMaximizeVariable(XVariable *x) {
#ifdef _VERBOSE_
  cout << "\n    objective: maximize variable " << x->id << endl;
#endif

  goal = new Goal(Goal::MAXIMIZATION, variable[x->id].get_var());
}

void XCSP3MistralCallbacks::buildObjectiveMinimize(ExpressionObjective type,
                                                   vector<XVariable *> &list,
                                                   vector<int> &coefs) {
#ifdef _VERBOSE_
  cout << "\n    objective: minimize "
       << (type == SUM_O ? "sum" : "unknown constraint") << endl;
#endif

  VarArray scope;
  getVariables(list, scope);
  Variable objective;

  if (type == EXPRESSION_O) {
    cout << "s UNSUPPORTED" << _ID_(": EXPRESSION_O") << "\n";
    exit(1);
  } else if (type == SUM_O) {

    objective = Sum(scope, coefs).get_var();

  } else if (type == PRODUCT_O) {
    cout << "s UNSUPPORTED" << _ID_(": PRODUCT_O") << "\n";
    exit(1);
  } else if (type == MINIMUM_O) {

    objective = Min(scope).get_var();

  } else if (type == MAXIMUM_O) {

    objective = Max(scope).get_var();

  } else if (type == NVALUES_O) {

    int min = INFTY;
    int max = -INFTY;

    for (auto X : scope) {
      if (min > X.get_min()) {
        min = X.get_min();
      }
      if (max < X.get_max()) {
        max = X.get_max();
      }
    }

    VarArray used;
    for (int i = min; i <= max; ++i) {
      used.add(Occurrence(scope, i));
    }

    objective = BoolSum(used).get_var();

  } else if (type == LEX_O) {
    cout << "s UNSUPPORTED" << _ID_(": LEX_O") << "\n";
    exit(1);
  }

  solver.add(Free(objective));

#ifdef _VERBOSE_
  std::cout << objective << " in " << objective.get_domain() << std::endl;
#endif

  goal = new Goal(Goal::MINIMIZATION, objective.get_var());

  // XCSP3CoreCallbacks::buildObjectiveMinimize(type, list, coefs);
}

void XCSP3MistralCallbacks::buildObjectiveMaximize(ExpressionObjective type,
                                                   vector<XVariable *> &list,
                                                   vector<int> &coefs) {
#ifdef _VERBOSE_
  cout << "\n    objective: maximize "
       << (type == SUM_O ? "sum" : "unknown constraint") << endl;
#endif

  VarArray scope;
  getVariables(list, scope);
  Variable objective;

  if (type == EXPRESSION_O) {
    cout << "s UNSUPPORTED" << _ID_(": EXPRESSION_O") << "\n";
    exit(1);
  } else if (type == SUM_O) {

    objective = Sum(scope, coefs);

  } else if (type == PRODUCT_O) {
    cout << "s UNSUPPORTED" << _ID_(": PRODUCT_O") << "\n";
    exit(1);
  } else if (type == MINIMUM_O) {
    cout << "s UNSUPPORTED" << _ID_(": MINIMUM_O") << "\n";
    exit(1);
  } else if (type == MAXIMUM_O) {
    cout << "s UNSUPPORTED" << _ID_(": MAXIMUM_O") << "\n";
    exit(1);
  } else if (type == NVALUES_O) {
    cout << "s UNSUPPORTED" << _ID_(": NVALUES_O") << "\n";
    exit(1);
  } else if (type == LEX_O) {
    cout << "s UNSUPPORTED" << _ID_(": LEX_O") << "\n";
    exit(1);
  }

  solver.add(Free(objective));

#ifdef _VERBOSE_
  std::cout << objective << " in " << objective.get_domain() << std::endl;
#endif

  goal = new Goal(Goal::MAXIMIZATION, objective.get_var());

  // XCSP3CoreCallbacks::buildObjectiveMaximize(type, list, coefs);
}

void XCSP3MistralCallbacks::buildObjectiveMinimize(ExpressionObjective type,
                                                   vector<XVariable *> &list) {
#ifdef _VERBOSE_
  cout << "\n    objective: minimize "
       << (type == SUM_O
               ? "sum"
               : (type == MINIMUM_O
                      ? "min"
                      : (type == MAXIMUM_O ? "max" : "unknown constraint")))
       << endl;
#endif

  VarArray scope;
  getVariables(list, scope);
  Variable objective;

  if (type == EXPRESSION_O) {
    cout << "s UNSUPPORTED" << _ID_(": EXPRESSION_O") << "\n";
    exit(1);
  } else if (type == SUM_O) {

    objective = Sum(scope);

  } else if (type == PRODUCT_O) {
    cout << "s UNSUPPORTED" << _ID_(": PRODUCT_O") << "\n";
    exit(1);
  } else if (type == MINIMUM_O) {

    objective = Min(scope);

  } else if (type == MAXIMUM_O) {

    objective = Max(scope);

  } else if (type == NVALUES_O) {
    cout << "s UNSUPPORTED" << _ID_(": NVALUES_O") << "\n";
    exit(1);
  } else if (type == LEX_O) {
    cout << "s UNSUPPORTED" << _ID_(": LEX_O") << "\n";
    exit(1);
  }

  solver.add(Free(objective));

#ifdef _VERBOSE_
  std::cout << objective << " in " << objective.get_domain() << std::endl;
#endif

  goal = new Goal(Goal::MINIMIZATION, objective.get_var());

  // XCSP3CoreCallbacks::buildObjectiveMinimize(type, list);
}

void XCSP3MistralCallbacks::buildObjectiveMaximize(ExpressionObjective type,
                                                   vector<XVariable *> &list) {
#ifdef _VERBOSE_
  cout << "\n    objective: maximize "
       << (type == SUM_O
               ? "sum"
               : (type == MINIMUM_O
                      ? "min"
                      : (type == MAXIMUM_O ? "max" : "unknown constraint")))
       << endl;
#endif

  VarArray scope;
  getVariables(list, scope);
  Variable objective;

  if (type == EXPRESSION_O) {
    cout << "s UNSUPPORTED" << _ID_(": EXPRESSION_O") << "\n";
    exit(1);
  } else if (type == SUM_O) {

    objective = Sum(scope);

  } else if (type == PRODUCT_O) {
    cout << "s UNSUPPORTED" << _ID_(": PRODUCT_O") << "\n";
    exit(1);
  } else if (type == MINIMUM_O) {

    objective = Min(scope);

  } else if (type == MAXIMUM_O) {

    objective = Max(scope);

  } else if (type == NVALUES_O) {
    cout << "s UNSUPPORTED" << _ID_(": NVALUES_O") << "\n";
    exit(1);
  } else if (type == LEX_O) {
    cout << "s UNSUPPORTED" << _ID_(": LEX_O") << "\n";
    exit(1);
  }

// solver.add(Free(objective));

#ifdef _VERBOSE_
  std::cout << objective << std::endl;
#endif

  goal = new Goal(Goal::MINIMIZATION, objective.get_var());

  // XCSP3CoreCallbacks::buildObjectiveMaximize(type, list);
}

Variable XCSP3MistralCallbacks::postExpression(XCSP3Mistral::Node *n,
                                               int level) {
  Variable rv;

  if (n->type == XCSP3Mistral::NT_VARIABLE) {
    assert(level > 0);

    XCSP3Mistral::NodeVariable *nv = (XCSP3Mistral::NodeVariable *)n;
    rv = variable[nv->var];
    ++initial_degree[id_map[nv->var]];
  }

  if (n->type == XCSP3Mistral::NT_CONSTANT) {
    assert(level > 0);

    XCSP3Mistral::NodeConstant *nc = (XCSP3Mistral::NodeConstant *)n;
    rv = Variable(nc->val, nc->val);
  }

  XCSP3Mistral::NodeOperator *fn = (XCSP3Mistral::NodeOperator *)n;

  if (fn->type == XCSP3Mistral::NT_EQ) {

    if (fn->args.size() == 2) {
      Variable x1 = postExpression(fn->args[0]);
      Variable x2 = postExpression(fn->args[1]);
      rv = (x1 == x2);
    } else {
      VarArray scope;
      for (auto expr : fn->args) {
        scope.add(postExpression(expr));
      }
      VarArray equalities;
      for (size_t i = 1; i < scope.size; ++i) {
        equalities.add(scope[i - 1] == scope[i]);
      }

      rv = (BoolSum(equalities) == (scope.size - 1));
    }
  }

  if (fn->type == XCSP3Mistral::NT_NE) {

    Variable x1 = postExpression(fn->args[0]);
    Variable x2 = postExpression(fn->args[1]);
    rv = (x1 != x2);
  }

  if (fn->type == XCSP3Mistral::NT_GE) {

    Variable x1 = postExpression(fn->args[0]);
    Variable x2 = postExpression(fn->args[1]);
    rv = (x1 >= x2);
  }

  if (fn->type == XCSP3Mistral::NT_GT) {

    Variable x1 = postExpression(fn->args[0]);
    Variable x2 = postExpression(fn->args[1]);
    rv = (x1 > x2);
  }

  if (fn->type == XCSP3Mistral::NT_LE) {

    Variable x1 = postExpression(fn->args[0]);
    Variable x2 = postExpression(fn->args[1]);
    rv = (x1 <= x2);
  }

  if (fn->type == XCSP3Mistral::NT_LT) {

    Variable x1 = postExpression(fn->args[0]);
    Variable x2 = postExpression(fn->args[1]);
    rv = (x1 < x2);
  }

  if (fn->type == XCSP3Mistral::NT_IMP) { // IMP(X,Y) = NOT X OR Y

    Variable x1 = postExpression(fn->args[0]);
    Variable x2 = postExpression(fn->args[1]);
    rv = (x1 <= x2);
  }

  if (fn->type == XCSP3Mistral::NT_OR) {

    if (fn->args.size() == 2) {
      Variable x1 = postExpression(fn->args[0]);
      Variable x2 = postExpression(fn->args[1]);
      rv = (x1 || x2);
    } else {
      VarArray scope;
      for (auto x : fn->args) {
        scope.add(postExpression(x));
      }
      rv = (Sum(scope) > 0);
    }
    //
    // assert(fn->args.size() == 2);
    // 	Variable x1 = postExpression(fn->args[0]);
    //     		Variable x2 = postExpression(fn->args[1]);
    // 	rv = ( x1 || x2 );
  }

  if (fn->type == XCSP3Mistral::NT_AND) {

    if (fn->args.size() == 2) {
      Variable x1 = postExpression(fn->args[0]);
      Variable x2 = postExpression(fn->args[1]);
      rv = (x1 && x2);
    } else {
      VarArray scope;
      for (auto x : fn->args) {
        scope.add(postExpression(x));
      }
      rv = (Sum(scope) == scope.size);
    }
  }

  if (fn->type == XCSP3Mistral::NT_NOT) {

    Variable x1 = postExpression(fn->args[0]);
    rv = (!x1);
  }

  if (fn->type == XCSP3Mistral::NT_IFF) {

    Variable x1 = postExpression(fn->args[0]);
    Variable x2 = postExpression(fn->args[1]);
    rv = (x1 == x2);
  }

  if (fn->type == XCSP3Mistral::NT_XOR) {

    if (fn->args.size() == 2) {
      Variable x1 = postExpression(fn->args[0]);
      Variable x2 = postExpression(fn->args[1]);
      rv = (x1 != x2);
    } else {
      VarArray scope;
      for (auto x : fn->args) {
        scope.add(postExpression(x));
      }
      rv = Parity(scope, 1);
    }
  }

  // function stuff
  if (fn->type == XCSP3Mistral::NT_NEG) {

    Variable x1 = postExpression(fn->args[0]);
    rv = -x1;
  }

  if (fn->type == XCSP3Mistral::NT_ABS) {

    Variable x1 = postExpression(fn->args[0]);
    rv = Abs(x1);
  }

  if (fn->type == XCSP3Mistral::NT_SUB) {

    Variable x1 = postExpression(fn->args[0]);
    Variable x2 = postExpression(fn->args[1]);
    rv = (x1 - x2);
  }

  if (fn->type == XCSP3Mistral::NT_DIST) { // Simulate DIST(X,Y) = ABS(SUB(X,Y))

    Variable x1 = postExpression(fn->args[0]);
    Variable x2 = postExpression(fn->args[1]);
    rv = Abs(x1 - x2);
  }

  if (fn->type == XCSP3Mistral::NT_ADD) {

    if (fn->args.size() == 2) {
      Variable x1 = postExpression(fn->args[0]);
      Variable x2 = postExpression(fn->args[1]);
      rv = (x1 + x2);
    } else {
      VarArray scope;
      Vector<int> weights;

      map<int, int> id_map;
      bool weighted = false;
      for (auto x : fn->args) {
        Variable y = postExpression(x);

        if (y.id() >= 0 && id_map.count(y.id())) {
          weights[id_map[y.id()]]++;
          weighted = true;
        } else {
          id_map[y.id()] = weights.size;
          weights.add(1);
          scope.add(y);
        }
      }

      if (weighted) {
        rv = Sum(scope, weights);
      } else {
        rv = Sum(scope);
      }
    }
  }

  if (fn->type == XCSP3Mistral::NT_MULT) {

    if (fn->args.size() == 2) {
      Variable x1 = postExpression(fn->args[0]);
      Variable x2 = postExpression(fn->args[1]);
      rv = (x1 * x2);
    } else {
      rv = postExpression(fn->args[0]);
      for (size_t i = 1; i < fn->args.size(); ++i) {
        rv = (rv * postExpression(fn->args[i]));
      }
    }
  }

  if (fn->type == XCSP3Mistral::NT_DIV) {

    assert(fn->args.size() == 2);

    Variable x1 = postExpression(fn->args[0]);
    Variable x2 = postExpression(fn->args[1]);
    rv = (x1 / x2);
  }

  if (fn->type == XCSP3Mistral::NT_MOD) {

    assert(fn->args.size() == 2);

    Variable x1 = postExpression(fn->args[0]);
    Variable x2 = postExpression(fn->args[1]);
    rv = (x1 % x2);
  }

  if (fn->type == XCSP3Mistral::NT_MIN) {
    assert(level > 0);

    if (fn->args.size() == 2) {
      Variable x1 = postExpression(fn->args[0]);
      Variable x2 = postExpression(fn->args[1]);
      rv = Min(x1, x2);
    } else {
      VarArray scope;
      for (auto x : fn->args) {
        scope.add(postExpression(x));
      }
      rv = Min(scope);
    }
  }
  if (fn->type == XCSP3Mistral::NT_MAX) {
    assert(level > 0);

    if (fn->args.size() == 2) {
      Variable x1 = postExpression(fn->args[0]);
      Variable x2 = postExpression(fn->args[1]);
      rv = Max(x1, x2);
    } else {
      VarArray scope;
      for (auto x : fn->args) {
        scope.add(postExpression(x));
      }
      rv = Max(scope);
    }
  }

  if (fn->type == XCSP3Mistral::NT_IF) {
    assert(level > 0);

    Variable x = postExpression(fn->args[0]);
    VarArray A;
    A.add(postExpression(fn->args[2]));
    A.add(postExpression(fn->args[1]));
    rv = Element(A, x);
  }

  if (fn->type == XCSP3Mistral::NT_IN) {

    Variable x = postExpression(fn->args[0]);
    Variable y = postExpression(fn->args[1]);

    if (y.is_void())
      rv = Member(x, util);
    else
      rv = (x == y);
  }

  if (fn->type == XCSP3Mistral::NT_SET) {
    assert(level > 0);

    vector<int> &elements(util);
    elements.clear();

    for (auto x : fn->args) {
      if (x->type == XCSP3Mistral::NT_CONSTANT) {
        elements.push_back(((XCSP3Mistral::NodeConstant *)x)->val);
      } else {
        break;
      }
    }

    if (elements.size() != fn->args.size()) {
      VarArray scope;
      for (auto x : fn->args) {
        scope.add(postExpression(x));
      }

      Variable any(0, scope.size - 1);
      rv = Element(scope, any);
    } else {
      rv = Variable();
    }
  }

  return rv;
}

Variable XCSP3MistralCallbacks::postExpression(Node *n, int level) {
  Variable rv;

  if (n->type == OVAR) {
    assert(level > 0);

    NodeVariable *nv = (NodeVariable *)n;
    rv = variable[nv->var];
    ++initial_degree[id_map[nv->var]];

#ifdef _VERBOSE_
    for (int l = 0; l < level; ++l)
      std::cout << "  ";
    std::cout << rv << " in " << rv.get_domain() << std::endl;
#endif

  }

  else if (n->type == ODECIMAL) {
    assert(level > 0);

    NodeConstant *nc = (NodeConstant *)n;
    rv = Variable(nc->val, nc->val);

#ifdef _VERBOSE_
    for (int l = 0; l < level; ++l)
      std::cout << "  ";
    std::cout << rv << " in " << rv.get_domain() << std::endl;
#endif

  }

  else {

    NodeOperator *fn = (NodeOperator *)n;

    if (fn->type == OEQ) {

      if (fn->parameters.size() == 2) {
        Variable x1 = postExpression(fn->parameters[0], level + 1);
        Variable x2 = postExpression(fn->parameters[1], level + 1);
        rv = (x1 == x2);
      } else {
        VarArray scope;
        for (auto expr : fn->parameters) {
          scope.add(postExpression(expr, level + 1));
        }
        VarArray equalities;
        for (size_t i = 1; i < scope.size; ++i) {
          equalities.add(scope[i - 1] == scope[i]);
        }

        rv = (BoolSum(equalities) == (scope.size - 1));
      }

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == ONE) {

      Variable x1 = postExpression(fn->parameters[0], level + 1);
      Variable x2 = postExpression(fn->parameters[1], level + 1);
      rv = (x1 != x2);

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OGE) {

      Variable x1 = postExpression(fn->parameters[0], level + 1);
      Variable x2 = postExpression(fn->parameters[1], level + 1);
      rv = (x1 >= x2);

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OGT) {

      Variable x1 = postExpression(fn->parameters[0], level + 1);
      Variable x2 = postExpression(fn->parameters[1], level + 1);
      rv = (x1 > x2);

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OLE) {

      Variable x1 = postExpression(fn->parameters[0], level + 1);
      Variable x2 = postExpression(fn->parameters[1], level + 1);
      rv = (x1 <= x2);

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OLT) {

      Variable x1 = postExpression(fn->parameters[0], level + 1);
      Variable x2 = postExpression(fn->parameters[1], level + 1);
      rv = (x1 < x2);

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OIMP) { // IMP(X,Y) = NOT X OR Y

      Variable x1 = postExpression(fn->parameters[0], level + 1);
      Variable x2 = postExpression(fn->parameters[1], level + 1);
      rv = (x1 <= x2);

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OOR) {

      if (fn->parameters.size() == 2) {
        Variable x1 = postExpression(fn->parameters[0], level + 1);
        Variable x2 = postExpression(fn->parameters[1], level + 1);
        rv = (x1 || x2);
      } else {
        VarArray scope;
        for (auto x : fn->parameters) {
          scope.add(postExpression(x, level + 1));
        }
        rv = (Sum(scope) > 0);
      }
//
// assert(fn->parameters.size() == 2);
// 	Variable x1 = postExpression(fn->parameters[0]);
//     		Variable x2 = postExpression(fn->parameters[1]);
// 	rv = ( x1 || x2 );

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OAND) {

      if (fn->parameters.size() == 2) {
        Variable x1 = postExpression(fn->parameters[0], level + 1);
        Variable x2 = postExpression(fn->parameters[1], level + 1);
        rv = (x1 && x2);
      } else {
        VarArray scope;
        for (auto x : fn->parameters) {
          scope.add(postExpression(x, level + 1));
        }
        rv = (Sum(scope) == scope.size);
      }

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == ONOT) {

      Variable x1 = postExpression(fn->parameters[0], level + 1);
      rv = (!x1);

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OIFF) {

      Variable x1 = postExpression(fn->parameters[0], level + 1);
      Variable x2 = postExpression(fn->parameters[1], level + 1);
      rv = (x1 == x2);

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OXOR) {

      if (fn->parameters.size() == 2) {
        Variable x1 = postExpression(fn->parameters[0], level + 1);
        Variable x2 = postExpression(fn->parameters[1], level + 1);
        rv = (x1 != x2);
      } else {
        VarArray scope;
        for (auto x : fn->parameters) {
          scope.add(postExpression(x, level + 1));
        }
        rv = Parity(scope, 1);
      }

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    // function stuff
    else if (fn->type == ONEG) {

      Variable x1 = postExpression(fn->parameters[0], level + 1);
      rv = -x1;

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OABS) {

      Variable x1 = postExpression(fn->parameters[0], level + 1);
      rv = Abs(x1);

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OSUB) {

      Variable x1 = postExpression(fn->parameters[0], level + 1);
      Variable x2 = postExpression(fn->parameters[1], level + 1);
      rv = (x1 - x2);

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == ODIST) { // Simulate DIST(X,Y) = ABS(SUB(X,Y))

      Variable x1 = postExpression(fn->parameters[0], level + 1);
      Variable x2 = postExpression(fn->parameters[1], level + 1);
      rv = Abs(x1 - x2);

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OADD) {

      if (fn->parameters.size() == 2) {
        Variable x1 = postExpression(fn->parameters[0], level + 1);
        Variable x2 = postExpression(fn->parameters[1], level + 1);
        rv = (x1 + x2);
      } else {
        VarArray scope;
        Vector<int> weights;

        map<int, int> id_map;
        bool weighted = false;
        for (auto x : fn->parameters) {
          Variable y = postExpression(x, level + 1);

          if (y.id() >= 0 && id_map.count(y.id())) {
            weights[id_map[y.id()]]++;
            weighted = true;
          } else {
            id_map[y.id()] = weights.size;
            weights.add(1);
            scope.add(y);
          }
        }

        if (weighted) {
          rv = Sum(scope, weights);
        } else {
          rv = Sum(scope);
        }
      }

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif
    }

    else if (fn->type == OMUL) {

      if (fn->parameters.size() == 2) {
        Variable x1 = postExpression(fn->parameters[0], level + 1);
        Variable x2 = postExpression(fn->parameters[1], level + 1);
        rv = (x1 * x2);
      } else {
        rv = postExpression(fn->parameters[0], level + 1);
        for (size_t i = 1; i < fn->parameters.size(); ++i) {
          rv = (rv * postExpression(fn->parameters[i], level + 1));
        }
      }

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == ODIV) {

      assert(fn->parameters.size() == 2);

      Variable x1 = postExpression(fn->parameters[0], level + 1);
      Variable x2 = postExpression(fn->parameters[1], level + 1);
      rv = (x1 / x2);

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OMOD) {

      assert(fn->parameters.size() == 2);

      Variable x1 = postExpression(fn->parameters[0], level + 1);
      Variable x2 = postExpression(fn->parameters[1], level + 1);
      rv = (x1 % x2);

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OMIN) {
      assert(level > 0);

      if (fn->parameters.size() == 2) {
        Variable x1 = postExpression(fn->parameters[0], level + 1);
        Variable x2 = postExpression(fn->parameters[1], level + 1);
        rv = Min(x1, x2);
      } else {
        VarArray scope;
        for (auto x : fn->parameters) {
          scope.add(postExpression(x, level + 1));
        }
        rv = Min(scope);
      }

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    } else if (fn->type == OMAX) {
      assert(level > 0);

      if (fn->parameters.size() == 2) {
        Variable x1 = postExpression(fn->parameters[0], level + 1);
        Variable x2 = postExpression(fn->parameters[1], level + 1);
        rv = Max(x1, x2);
      } else {
        VarArray scope;
        for (auto x : fn->parameters) {
          scope.add(postExpression(x, level + 1));
        }
        rv = Max(scope);
      }

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OIF) {
      assert(level > 0);

      Variable x = postExpression(fn->parameters[0], level + 1);
      VarArray A;
      A.add(postExpression(fn->parameters[2], level + 1));
      A.add(postExpression(fn->parameters[1], level + 1));

      rv = Element(A, x);

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OIN) {

      Variable x = postExpression(fn->parameters[0], level + 1);
      Variable y = postExpression(fn->parameters[1], level + 1);

      if (y.is_void())
        rv = Member(x, util);
      else
        rv = (x == y);

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else if (fn->type == OSET) {
      assert(level > 0);

      vector<int> &elements(util);
      elements.clear();

      for (auto x : fn->parameters) {
        if (x->type == ODECIMAL) {
          elements.push_back(((NodeConstant *)x)->val);
        } else {
          break;
        }
      }

      if (elements.size() != fn->parameters.size()) {
        VarArray scope;
        for (auto x : fn->parameters) {
          scope.add(postExpression(x, level + 1));
        }

        Variable any(0, scope.size - 1);
        rv = Element(scope, any);
      } else {
        rv = Variable();
      }

#ifdef _VERBOSE_
      for (int l = 0; l < level; ++l)
        std::cout << "  ";
      std::cout << rv << std::endl;
#endif

    }

    else {
      std::cout << "unknown operator: " << operatorToString(fn->type) << "\n";
      exit(1);
    }
  }

  return rv;
}

#endif //COSOCO_XCSP3PRINTCALLBACKS_H
