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

#include "Tree.hpp"

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

namespace XCSP3Core {

    class XCSP3MistralCallbacks : public XCSP3CoreCallbacks {
    public:
		
			Mistral::Solver& solver;
			
			map<string, Mistral::VarArray> array;
			map<string, Mistral::Variable> variable;
			
			vector<Mistral::VarArray> v_array;
			vector<Mistral::Variable> v_variable;
			
			
			
			
        XCSP3MistralCallbacks(Mistral::Solver& s);

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


        virtual void buildVariableInteger(string id, int minValue, int maxValue) override;

        virtual void buildVariableInteger(string id, vector<int> &values) override;

        virtual void buildConstraintExtension(string id, vector<XVariable *> list, vector <vector<int>> &tuples, bool support, bool hasStar) override;

        virtual void buildConstraintExtension(string id, XVariable *variable, vector<int> &tuples, bool support, bool hasStar) override;

        virtual void buildConstraintExtensionAs(string id, vector<XVariable *> list, bool support, bool hasStar) override;

        virtual void buildConstraintIntension(string id, string expr) override;

        virtual void buildConstraintPrimitive(string id, OrderType op, XVariable *x, int k, XVariable *y) override;


        virtual void buildConstraintRegular(string id, vector<XVariable *> &list, string st, vector <string> &final, vector <XTransition> &transitions) override;

        virtual void buildConstraintMDD(string id, vector<XVariable *> &list, vector <XTransition> &transitions) override;

        virtual void buildConstraintAlldifferent(string id, vector<XVariable *> &list) override;

        virtual void buildConstraintAlldifferentExcept(string id, vector<XVariable *> &list, vector<int> &except) override;

        virtual void buildConstraintAlldifferentList(string id, vector <vector<XVariable *>> &lists) override;

        virtual void buildConstraintAlldifferentMatrix(string id, vector <vector<XVariable *>> &matrix) override;

        virtual void buildConstraintAllEqual(string id, vector<XVariable *> &list) override;

        virtual void buildConstraintNotAllEqual(string id, vector<XVariable *> &list) override;

        virtual void buildConstraintOrdered(string id, vector<XVariable *> &list, OrderType order) override;

        virtual void buildConstraintLex(string id, vector <vector<XVariable *>> &lists, OrderType order) override;

        virtual void buildConstraintLexMatrix(string id, vector <vector<XVariable *>> &matrix, OrderType order) override;

        virtual void buildConstraintSum(string id, vector<XVariable *> &list, vector<int> &coeffs, XCondition &cond) override;

        virtual void buildConstraintSum(string id, vector<XVariable *> &list, XCondition &cond) override;

        virtual void buildConstraintSum(string id, vector<XVariable *> &list, vector<XVariable *> &coeffs, XCondition &cond) override;

        virtual void buildConstraintAtMost(string id, vector<XVariable *> &list, int value, int k) override;

        virtual void buildConstraintAtLeast(string id, vector<XVariable *> &list, int value, int k) override;

        virtual void buildConstraintExactlyK(string id, vector<XVariable *> &list, int value, int k) override;

        virtual void buildConstraintAmong(string id, vector<XVariable *> &list, vector<int> &values, int k) override;

        virtual void buildConstraintExactlyVariable(string id, vector<XVariable *> &list, int value, XVariable *x) override;

        virtual void buildConstraintCount(string id, vector<XVariable *> &list, vector<int> &values, XCondition &xc) override;

        virtual void buildConstraintCount(string id, vector<XVariable *> &list, vector<XVariable *> &values, XCondition &xc) override;

        virtual void buildConstraintNValues(string id, vector<XVariable *> &list, vector<int> &except, XCondition &xc) override;

        virtual void buildConstraintNValues(string id, vector<XVariable *> &list, XCondition &xc) override;

        virtual void buildConstraintCardinality(string id, vector<XVariable *> &list, vector<int> values, vector<int> &occurs, bool closed) override;

        virtual void buildConstraintCardinality(string id, vector<XVariable *> &list, vector<int> values, vector<XVariable *> &occurs,
                                                bool closed) override;

        virtual void buildConstraintCardinality(string id, vector<XVariable *> &list, vector<int> values, vector <XInterval> &occurs,
                                                bool closed) override;

        virtual void buildConstraintCardinality(string id, vector<XVariable *> &list, vector<XVariable *> values, vector<int> &occurs,
                                                bool closed) override;

        virtual void buildConstraintCardinality(string id, vector<XVariable *> &list, vector<XVariable *> values, vector<XVariable *> &occurs,
                                                bool closed) override;

        virtual void buildConstraintCardinality(string id, vector<XVariable *> &list, vector<XVariable *> values, vector <XInterval> &occurs,
                                                bool closed) override;

        virtual void buildConstraintMinimum(string id, vector<XVariable *> &list, XCondition &xc) override;

        virtual void buildConstraintMinimum(string id, vector<XVariable *> &list, XVariable *index, int startIndex, RankType rank,
                                            XCondition &xc) override;

        virtual void buildConstraintMaximum(string id, vector<XVariable *> &list, XCondition &xc) override;

        virtual void buildConstraintMaximum(string id, vector<XVariable *> &list, XVariable *index, int startIndex, RankType rank,
                                            XCondition &xc) override;

        virtual void buildConstraintElement(string id, vector<XVariable *> &list, int value) override;

        virtual void buildConstraintElement(string id, vector<XVariable *> &list, XVariable *value) override;

        virtual void buildConstraintElement(string id, vector<XVariable *> &list, int startIndex, XVariable *index, RankType rank, int value) override;

        virtual void buildConstraintElement(string id, vector<XVariable *> &list, int startIndex, XVariable *index, RankType rank, XVariable *value) override;

        virtual void buildConstraintChannel(string id, vector<XVariable *> &list, int startIndex) override;

        virtual void buildConstraintChannel(string id, vector<XVariable *> &list1, int startIndex1, vector<XVariable *> &list2,
                                            int startIndex2) override;

        virtual void buildConstraintChannel(string id, vector<XVariable *> &list, int startIndex, XVariable *value) override;

        virtual void buildConstraintStretch(string id, vector<XVariable *> &list, vector<int> &values, vector <XInterval> &widths) override;

        virtual void buildConstraintStretch(string id, vector<XVariable *> &list, vector<int> &values, vector <XInterval> &widths, vector <vector<int>> &patterns) override;

        virtual void buildConstraintNoOverlap(string id, vector<XVariable *> &origins, vector<int> &lengths, bool zeroIgnored) override;

        virtual void buildConstraintNoOverlap(string id, vector<XVariable *> &origins, vector<XVariable *> &lengths, bool zeroIgnored) override;

        virtual void buildConstraintNoOverlap(string id, vector <vector<XVariable *>> &origins, vector <vector<int>> &lengths, bool zeroIgnored) override;

        virtual void buildConstraintNoOverlap(string id, vector <vector<XVariable *>> &origins, vector <vector<XVariable *>> &lengths, bool zeroIgnored) override;

        virtual void buildConstraintInstantiation(string id, vector<XVariable *> &list, vector<int> &values) override;

        virtual void buildObjectiveMinimizeExpression(string expr) override;

        virtual void buildObjectiveMaximizeExpression(string expr) override;


        virtual void buildObjectiveMinimizeVariable(XVariable *x) override;


        virtual void buildObjectiveMaximizeVariable(XVariable *x) override;


        virtual void buildObjectiveMinimize(ExpressionObjective type, vector<XVariable *> &list, vector<int> &coefs) override;


        virtual void buildObjectiveMaximize(ExpressionObjective type, vector<XVariable *> &list, vector<int> &coefs) override;


        virtual void buildObjectiveMinimize(ExpressionObjective type, vector<XVariable *> &list) override;


        virtual void buildObjectiveMaximize(ExpressionObjective type, vector<XVariable *> &list) override;


				void getVariables(vector < XVariable* > &list, Mistral::Vector<Mistral::Variable>& scope);

				Mistral::Variable postExpression(Node *n, bool isRoot = false);

    };


}

using namespace XCSP3Core;
using namespace Mistral;


XCSP3MistralCallbacks::XCSP3MistralCallbacks(Solver& s) : XCSP3CoreCallbacks(), solver(s) {

}


void XCSP3MistralCallbacks::getVariables(vector < XVariable* > &list, Vector<Variable>& scope) {	
    for(int i = 0; i < list.size(); i++)
				scope.add(variable[list[i]->id]);
}


template<class T>
void displayList(vector <T> &list, string separator = " ") {
    if(list.size() > 8) {
        for(int i = 0; i < 3; i++)
            cout << list[i] << separator;
        cout << " ... ";
        for(int i = list.size() - 4; i < list.size(); i++)
            cout << list[i] << separator;
        cout << endl;
        return;
    }
    for(int i = 0; i < list.size(); i++)
        cout << list[i] << separator;
    cout << endl;
}


void displayList(vector < XVariable * > &list, string
separator = " "
) {
if(list.

size()

> 8) {
for(
int i = 0;
i < 3; i++)
cout << list[i]->id <<
separator;
cout << " ... ";
for(
unsigned int i = list.size() - 4;
i<list.

size();

i++)
cout << list[i]->id <<
separator;
cout <<
endl;
return;
}
for(
unsigned int i = 0;
i<list.

size();

i++)
cout << list[i]->id <<
separator;
cout <<
endl;
}


void XCSP3MistralCallbacks::beginInstance(InstanceType type) {
    cout << "Start Instance - type=" << type << endl;
}


void XCSP3MistralCallbacks::endInstance() {
    cout << "End SAX parsing " << endl;
}


void XCSP3MistralCallbacks::beginVariables() {
    cout << " start variables declaration" << endl;
}


void XCSP3MistralCallbacks::endVariables() {
    cout << " end variables declaration" << endl << endl;
}


void XCSP3MistralCallbacks::beginVariableArray(string id) {
    cout << "    array: " << id << endl;
		
		VarArray A;
		v_array.push_back(A);
}

void XCSP3MistralCallbacks::buildVariableInteger(string id, int minValue, int maxValue) {
    cout << "    var " << id << " : " << minValue << "..." << maxValue << endl;
		
		Variable X(minValue, maxValue);	
		variable[id] = X;
}


void XCSP3MistralCallbacks::buildVariableInteger(string id, vector<int> &values) {
    cout << "    var " << id << " : ";
    cout << "        ";
    displayList(values);

		Variable X(values);	
		variable[id] = X;
}



void XCSP3MistralCallbacks::endVariableArray() {
}


void XCSP3MistralCallbacks::beginConstraints() {
    cout << " start constraints declaration" << endl;
}


void XCSP3MistralCallbacks::endConstraints() {
    cout << "\n end constraints declaration" << endl << endl;
}


void XCSP3MistralCallbacks::beginGroup(string id) {
    cout << "   start group of constraint " << id << endl;
}


void XCSP3MistralCallbacks::endGroup() {
    cout << "   end group of constraint" << endl;
}


void XCSP3MistralCallbacks::beginBlock(string classes) {
    cout << "   start block of constraint classes = " << classes << endl;
}


void XCSP3MistralCallbacks::endBlock() {
    cout << "   end group of constraint" << endl;
}


void XCSP3MistralCallbacks::beginSlide(string id, bool circular) {
    cout << "   start slide " << id << endl;
}


void XCSP3MistralCallbacks::endSlide() {
    cout << "   end slide" << endl;
}


void XCSP3MistralCallbacks::beginObjectives() {
    cout << "   start Objective " << endl;
}


void XCSP3MistralCallbacks::endObjectives() {
    cout << "   end Objective " << endl;
}


void XCSP3MistralCallbacks::buildConstraintExtension(string id, vector<XVariable *> list, vector <vector<int>> &tuples, bool support, bool hasStar) {
    cout << "\n    extension constraint : " << id << endl;
    cout << "        " << (support ? "support" : "conflict") << " arity:" << list.size() << " nb tuples: " << tuples.size() << " star: " << hasStar << endl;
    cout << "        ";
    displayList(list);
		
		VarArray scope;
		getVariables(list, scope);
		
		TableExpression* tab = new TableExpression(scope);
		
		for( auto pt=begin(tuples); pt!=end(tuples); ++pt ) {
				tab->add(&((*pt)[0]));
		}
		
		solver.add( Variable(tab) );
}


void XCSP3MistralCallbacks::buildConstraintExtension(string id, XVariable *variable, vector<int> &tuples, bool support, bool hasStar) {
    cout << "\n    extension constraint with one variable: " << id << endl;
    cout << "        " << (support ? "support" : "conflict") << " nb tuples: " << tuples.size() << " star: " << hasStar << endl;
    cout << (*variable) << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintExtensionAs(string id, vector<XVariable *> list, bool support, bool hasStar) {
    cout << "\n    extension constraint similar as previous one: " << id << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintIntension(string id, string expr) {
    cout << "\n    intension constraint : " << id << " : " << expr << endl;
		
    Tree tree(expr);
    solver.add(postExpression(tree.root, true));
    tree.dispose();
		
}


void XCSP3MistralCallbacks::buildConstraintPrimitive(string id, OrderType op, XVariable *x, int k, XVariable *y) {
    cout << "\n   intension constraint " << id << ": " << x->id << (k >= 0 ? "+" : "") << k << " op " << y->id << endl;

		Variable X = variable[x->id];
		Variable Y = variable[y->id];
		
		 // LT, GE, GT, IN, EQ, NE)
		if(op == LE)
			solver.add( (X + k) <= Y );
		else if(op == LT)
			solver.add( (X + k) < Y );
		else if(op == GE)
			solver.add( (X + k) >= Y );
		else if(op == GT)
			solver.add( (X + k) > Y );
		else if(op == IN) {
			cout << "NOT IMPLEMENTED!" << endl;
			exit(1);	
		}
		else if(op == EQ)
			solver.add( (X + k) == Y );
		else if(op == NE)
			solver.add( (X + k) != Y );
		
}


void XCSP3MistralCallbacks::buildConstraintRegular(string id, vector<XVariable *> &list, string start, vector <string> &final, vector <XTransition> &transitions) {
    cout << "\n    regular constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        start: " << start << endl;
    cout << "        final: ";
    displayList(final, ",");
    cout << endl;
    cout << "        transitions: ";
    for(unsigned int i = 0; i < (transitions.size() > 4 ? 4 : transitions.size()); i++) {
        cout << "(" << transitions[i].from << "," << transitions[i].val << "," << transitions[i].to << ") ";
    }
    if(transitions.size() > 4) cout << "...";
    cout << endl;
		
		
		map<string, int > state_map;
		map<string, bool> final_map;
		
		for( auto t : transitions ) {
			if(state_map.find(t.from) == end(state_map)) {
					state_map[t.from] = state_map.size();
					final_map[t.from] = false;
			}
			if(state_map.find(t.to) == end(state_map)) {
					state_map[t.to] = state_map.size();
					final_map[t.to] = false;
			}
		}
		
		for( auto f : final ) {
			final_map[f] = true;
		}
		
		VarArray scope;
		getVariables(list, scope);

		VarArray state_var(scope.size-1, 0, state_map.size()-1);
		
		Vector<const int*>* trans = new Vector<const int*>;
		Vector<const int*>* first = new Vector<const int*>;
		Vector<const int*>* last  = new Vector<const int*>;
		
		for( const auto& t : transitions ) {
				int *tuple = new int[3];
				tuple[0] = state_map[t.from];
				tuple[1] = t.val;
				tuple[2] = state_map[t.to];
				trans->add(tuple);
				if(start == t.from) {
					int *pair = new int[2];
					pair[0] = tuple[1];
					pair[1] = tuple[2];
					first->add(pair);
				}
				if(final_map[t.to]) {
					int *pair = new int[2];
					pair[0] = tuple[0];
					pair[1] = tuple[1];
					last->add(pair);
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
		
		for(int i=1; i<state_var.size; ++i) {
			S.clear();
			S.add(state_var[i-1]);
			S.add(scope[i]);
			S.add(state_var[i]);
			solver.add(Table(S, trans));
		}
}


void XCSP3MistralCallbacks::buildConstraintMDD(string id, vector<XVariable *> &list, vector <XTransition> &transitions) {
    cout << "\n    mdd constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        transitions: ";
    for(unsigned int i = 0; i < (transitions.size() > 4 ? 4 : transitions.size()); i++) {
        cout << "(" << transitions[i].from << "," << transitions[i].val << "," << transitions[i].to << ") ";
    }
    if(transitions.size() > 4) cout << "...";
    cout << endl;
		
    for(unsigned int i = 0; i < transitions.size(); i++) {
        cout << "(" << transitions[i].from << "," << transitions[i].val << "," << transitions[i].to << ") ";
    }
    cout << endl;
		
		
		VarArray scope;
		getVariables(list, scope);
		
		string root = begin(transitions)->from;
		string final = transitions.rbegin()->to;
		
		map<string, int> layer_map;
		map<string, int> node_map;
		auto num_layers{0};
		
		for( const auto& t : transitions ) {
			if( t.to == final) break;
			auto l{0};
			if( t.from != root ) {
				l = layer_map[t.from]+1;
			}
			layer_map[t.to] = l;
			if(l>=num_layers) num_layers = l+1;
		}
		
		vector<int> layer_size;
		layer_size.resize(num_layers, 0);
		
		
		for( const auto& l : layer_map ) {
			node_map[l.first] = layer_size[l.second]++;
		}
		
		VarArray state_var;
		for( auto m : layer_size ) {
			state_var.add(Variable(0, m-1));
		}
		
		
		VarArray S;
		vector<TableExpression*> tables;
		for(int i=0; i<scope.size; ++i) {
			if(i>0) S.add(state_var[i-1]);
			S.add(scope[i]);
			if(i<num_layers) S.add(state_var[i]);
			
			tables.push_back(new TableExpression(S));
			S.clear();
		}
		
		
		int tuple[3];
		for( const auto& t : transitions ) {
			auto cons{0};
		
			if( t.from == root ) {
				tuple[0] = t.val;
				tuple[1] = node_map[t.to];
			} else {
				cons = layer_map[t.from]+1;
				
				if( t.to == final ) {
					if(cons != num_layers) {
						cout << "MDD TRANSITIONS TO FINAL FROM INTERNAL NOT IMPLEMENTED!" << endl;
						exit(1);
					}
					tuple[0] = node_map[t.from];
					tuple[1] = t.val;
				} else {
					tuple[0] = node_map[t.from];
					tuple[1] = t.val;
					tuple[2] = node_map[t.to];
				}
			}

			tables[cons]->add(tuple);
		
		}
	
		for( auto tab : tables ) {
			solver.add( Variable(tab) );
		}

}


void XCSP3MistralCallbacks::buildConstraintAlldifferent(string id, vector<XVariable *> &list) {
    cout << "\n    allDiff constraint" << id << endl;
    cout << "        ";
    displayList(list);
		
		VarArray scope;
		getVariables(list, scope);
	
		solver.add( AllDiff(scope) );
}


void XCSP3MistralCallbacks::buildConstraintAlldifferentExcept(string id, vector<XVariable *> &list, vector<int> &except) {
    cout << "\n    allDiff constraint with exceptions" << id << endl;
    cout << "        ";
    displayList(list);
    cout << "        Exceptions:";
    displayList(except);
		
		assert(except.size()==1);
		
		VarArray scope;
		getVariables(list, scope);
	
		solver.add( AllDiffExcept(scope, *begin(except)) );
}


void XCSP3MistralCallbacks::buildConstraintAlldifferentList(string id, vector <vector<XVariable *>> &lists) {
    cout << "\n    allDiff list constraint" << id << endl;
    for(unsigned int i = 0; i < (lists.size() < 4 ? lists.size() : 3); i++) {
        cout << "        ";
        displayList(lists[i]);

    }
		
		
		for( auto lpt1=begin(lists); lpt1!=end(lists); ++lpt1) {
			for( auto lpt2=lpt1+1; lpt2!=end(lists); ++lpt2) {
				assert(lpt1->size() == lpt2->size());
				VarArray differences;
				for(int i=0; i<lpt1->size(); ++i) {
					differences.add(variable[(*lpt1)[i]->id] != variable[(*lpt2)[i]->id]);
				}
				solver.add( Sum(differences)>0 );
			}
		}

}


void XCSP3MistralCallbacks::buildConstraintAlldifferentMatrix(string id, vector <vector<XVariable *>> &matrix) {
    cout << "\n    allDiff matrix constraint" << id << endl;
    for(unsigned int i = 0; i < matrix.size(); i++) {
        cout << "        ";
        displayList(matrix[i]);
    }
		
		VarArray scope;
    for(unsigned int i = 0; i < matrix.size(); i++) {
			getVariables(matrix[i], scope);
			solver.add(AllDiff(scope));
			scope.clear();
			for(unsigned int j = 0; j < matrix.size(); j++) {
				scope.add(variable[matrix[j][i]->id]);
			}
			solver.add(AllDiff(scope));
			scope.clear();
    }
		
}


void XCSP3MistralCallbacks::buildConstraintAllEqual(string id, vector<XVariable *> &list) {
    cout << "\n    allEqual constraint" << id << endl;
    cout << "        ";
    displayList(list);
		
		for(int i=1; i<list.size(); ++i) {
			solver.add(variable[list[i-1]->id] == variable[list[i]->id]);
		}
}


void XCSP3MistralCallbacks::buildConstraintNotAllEqual(string id, vector<XVariable *> &list) {
    cout << "\n    not allEqual constraint" << id << endl;
    cout << "        ";
    displayList(list);
		
		VarArray differences;
		
		for( auto xpt=begin(list); xpt!=end(list); ++xpt ) {
			for( auto ypt=xpt+1; ypt!=end(list); ++ypt) {
				differences.add( variable[(*xpt)->id] != variable[(*ypt)->id] );
			}
		}
		
		solver.add( Sum(differences)>0 );
}


void XCSP3MistralCallbacks::buildConstraintOrdered(string id, vector<XVariable *> &list, OrderType order) {
    cout << "\n    ordered constraint" << endl;
    string sep;
    if(order == LT) sep = " < ";
    if(order == LE) sep = " <= ";
    if(order == GT) sep = " > ";
    if(order == GE) sep = " >= ";
    cout << "        ";
    displayList(list, sep);
		
		VarArray X;
		getVariables(list, X);
	
		for(int i=1; i<X.size; ++i) {
				if(order == LT) solver.add( X[i-1] <  X[i] ); 
				if(order == LE) solver.add( X[i-1] <= X[i] ); 
				if(order == GT) solver.add( X[i-1] >  X[i] ); 
				if(order == GE) solver.add( X[i-1] >= X[i] ); 
		}
}


void XCSP3MistralCallbacks::buildConstraintLex(string id, vector <vector<XVariable *>> &lists, OrderType order) {
    cout << "\n    lex constraint   nb lists: " << lists.size() << endl;
    string sep;
    if(order == LT) sep = " < ";
    if(order == LE) sep = " <= ";
    if(order == GT) sep = " > ";
    if(order == GE) sep = " >= ";
    cout << "        operator: " << sep << endl;
    for(unsigned int i = 0; i < lists.size(); i++) {
        cout << "        list " << i << ": ";
        cout << "        ";
        displayList(lists[i], " ");
    }
		
		vector<VarArray> scope;
		for( auto l=begin(lists); l!=end(lists); ++l) {
			VarArray X;
			getVariables(*l, X);
			scope.push_back(X);
		}
		
		for(int i=1; i<lists.size(); ++i) {
			if(order == LT) solver.add( scope[i-1] <  scope[i] ); 
			if(order == LE) solver.add( scope[i-1] <= scope[i] ); 
			if(order == GT) solver.add( scope[i-1] >  scope[i] ); 
			if(order == GE) solver.add( scope[i-1] >= scope[i] ); 
		}
		
		
}


void XCSP3MistralCallbacks::buildConstraintLexMatrix(string id, vector <vector<XVariable *>> &matrix, OrderType order) {
    cout << "\n    lex matrix constraint   matrix  " << endl;
    string sep;
    if(order == LT) sep = " < ";
    if(order == LE) sep = " <= ";
    if(order == GT) sep = " > ";
    if(order == GE) sep = " >= ";

    for(unsigned int i = 0; i < (matrix.size() < 4 ? matrix.size() : 3); i++) {
        cout << "        ";
        displayList(matrix[i]);
    }
    cout << "        Order " << sep << endl;
		
		
		vector<VarArray> scope;
		for( auto l=begin(matrix); l!=end(matrix); ++l) {
			VarArray X;
			getVariables(*l, X);
			scope.push_back(X);
		}
		
		for(int i=1; i<scope.size(); ++i) {
			if(order == LT) solver.add( scope[i-1] <  scope[i] ); 
			if(order == LE) solver.add( scope[i-1] <= scope[i] ); 
			if(order == GT) solver.add( scope[i-1] >  scope[i] ); 
			if(order == GE) solver.add( scope[i-1] >= scope[i] ); 
		}
		
		scope.clear();
		
		for( int i=0; i<matrix.size(); ++i) {
			VarArray X;
			for( int j=0; j<matrix.size(); ++j) {	
				X.add(variable[matrix[j][i]->id]);
			}
			scope.push_back(X);
		}
		
		for(int i=1; i<scope.size(); ++i) {
			if(order == LT) solver.add( scope[i-1] <  scope[i] ); 
			if(order == LE) solver.add( scope[i-1] <= scope[i] ); 
			if(order == GT) solver.add( scope[i-1] >  scope[i] ); 
			if(order == GE) solver.add( scope[i-1] >= scope[i] ); 
		}
}


void XCSP3MistralCallbacks::buildConstraintSum(string id, vector<XVariable *> &list, vector<int> &coeffs, XCondition &cond) {
    cout << "\n        sum constraint:";
    if(list.size() > 8) {
        for(int i = 0; i < 3; i++)
            cout << (coeffs.size() == 0 ? 1 : coeffs[i]) << "*" << *(list[i]) << " ";
        cout << " ... ";
        for(unsigned int i = list.size() - 4; i < list.size(); i++)
            cout << (coeffs.size() == 0 ? 1 : coeffs[i]) << "*" << *(list[i]) << " ";
    } else {
        for(unsigned int i = 0; i < list.size(); i++)
            cout << (coeffs.size() == 0 ? 1 : coeffs[i]) << "*" << *(list[i]) << " ";
    }
    cout << cond << endl;
		
		vector<Variable> scope;
		for( auto x : list ) {
			scope.push_back( variable[x->id] );
		}
		
		if(cond.operandType == VARIABLE) {
			Variable total = variable[cond.var];
			
			if(cond.op == EQ) {
				solver.add( Sum(scope, coeffs, total) );
			} else if(cond.op == NE) {
				solver.add( Sum(scope, coeffs) != total);
			} else if(cond.op == LE) {
				solver.add( Sum(scope, coeffs) <= total);
			} else if(cond.op == LT) {
				solver.add( Sum(scope, coeffs) < total);
			} else if(cond.op == GE) {
				solver.add( Sum(scope, coeffs) >= total);
			} else if(cond.op == GT) {
				solver.add( Sum(scope, coeffs) >  total);
			}  
			
		} else if(cond.operandType == INTERVAL) {
			assert(cond.op == IN);
			
			solver.add( Sum(scope, coeffs, cond.min, cond.max) );
			
		} else {
			assert(cond.operandType == INTEGER);
				
			if(cond.op == EQ) {
				solver.add( Sum(scope, coeffs, cond.val, cond.val) );
			} else if(cond.op == NE) {
				Variable cst(cond.val);
				solver.add( Sum(scope, coeffs) != cst);
			} else if(cond.op == LE) {
				solver.add( Sum(scope, coeffs, -INFTY, cond.val) );
			} else if(cond.op == LT) {
				solver.add( Sum(scope, coeffs, -INFTY, cond.val-1) );
			} else if(cond.op == GE) {
				solver.add( Sum(scope, coeffs, cond.val, INFTY) );
			} else if(cond.op == GT) {
				solver.add( Sum(scope, coeffs, cond.val+1, INFTY) );
			}
		}
		
}


void XCSP3MistralCallbacks::buildConstraintSum(string id, vector<XVariable *> &list, XCondition &cond) {
    cout << "\n        unweighted sum constraint:";
    cout << "        ";
    displayList(list, "+");
    cout << cond << endl;

		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintSum(string id, vector<XVariable *> &list, vector<XVariable *> &coeffs, XCondition &cond) {
    cout << "\n        scalar sum constraint:";
    if(list.size() > 8) {
        for(int i = 0; i < 3; i++)
            cout << coeffs[i]->id << "*" << *(list[i]) << " ";
        cout << " ... ";
        for(unsigned int i = list.size() - 4; i < list.size(); i++)
            cout << coeffs[i]->id << "*" << *(list[i]) << " ";
    } else {
        for(unsigned int i = 0; i < list.size(); i++)
            cout << coeffs[i]->id << "*" << *(list[i]) << " ";
    }
    cout << cond << endl;

		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintAtMost(string id, vector<XVariable *> &list, int value, int k) {
    cout << "\n    AtMost constraint: val=" << value << " k=" << k << endl;
    cout << "        ";
    displayList(list);
		
		VarArray scope;
		getVariables(list, scope);
		
		VarArray equalities;
		for( auto X : scope ) {
			equalities.add( X == value );
		}
		
		solver.add( Sum(equalities) <= k );
	
}


void XCSP3MistralCallbacks::buildConstraintAtLeast(string id, vector<XVariable *> &list, int value, int k) {
    cout << "\n    Atleast constraint: val=" << value << " k=" << k << endl;
    cout << "        ";
    displayList(list);

		XCondition cond;
		cond.operandType = INTEGER;
		cond.op = GE;
		cond.val = k;
		vector<int> values;
		values.push_back(value);
		return buildConstraintCount(id, list, values, cond);
}


void XCSP3MistralCallbacks::buildConstraintExactlyK(string id, vector<XVariable *> &list, int value, int k) {
    cout << "\n    Exactly constraint: val=" << value << " k=" << k << endl;
    cout << "        ";
    displayList(list);

		XCondition cond;
		cond.operandType = INTEGER;
		cond.op = EQ;
		cond.val = k;
		vector<int> values;
		values.push_back(value);
		return buildConstraintCount(id, list, values, cond);
}


void XCSP3MistralCallbacks::buildConstraintAmong(string id, vector<XVariable *> &list, vector<int> &values, int k) {
    cout << "\n    Among constraint: k=" << k << endl;
    cout << "        ";
    displayList(list);
    cout << "        values:";
    displayList(values);
		
		XCondition cond;
		cond.operandType = INTEGER;
		cond.op = EQ;
		cond.val = k;
		return buildConstraintCount(id, list, values, cond);
}


void XCSP3MistralCallbacks::buildConstraintExactlyVariable(string id, vector<XVariable *> &list, int value, XVariable *x) {
    cout << "\n    Exactly Variable constraint: val=" << value << " variable=" << *x << endl;
    cout << "        ";
    displayList(list);
		
		XCondition cond;
		cond.operandType = VARIABLE;
		cond.op = EQ;
		cond.var = x->id;
		vector<int> values;
		values.push_back(value);
		return buildConstraintCount(id, list, values, cond);
}


void XCSP3MistralCallbacks::buildConstraintCount(string id, vector<XVariable *> &list, vector<int> &values, XCondition &cond) {
    cout << "\n    count constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        values: ";
    cout << "        ";
    displayList(values);
    cout << "        condition: " << cond << endl;
		
		VarArray scope;
		getVariables(list, scope);
		
		if(cond.operandType == VARIABLE) {
			Variable nocc = variable[cond.var];
			
			if(values.size() == 1) {
				int val = *begin(values);
				if(cond.op == EQ) {
					solver.add( Occurrence(scope, val) == nocc );
				} else if(cond.op == NE) {
					solver.add( Occurrence(scope, val) != nocc );
				} else if(cond.op == LE) {
					solver.add( Occurrence(scope, val) <= nocc );
				} else if(cond.op == LT) {
					solver.add( Occurrence(scope, val) <  nocc );
				} else if(cond.op == GE) {
					solver.add( Occurrence(scope, val) >= nocc );
				} else if(cond.op == GT) {
					solver.add( Occurrence(scope, val) >  nocc );
				}  
			} else {
				if(cond.op == EQ) {
					solver.add( Occurrence(scope, values) == nocc );
				} else if(cond.op == NE) {
					solver.add( Occurrence(scope, values) != nocc );
				} else if(cond.op == LE) {
					solver.add( Occurrence(scope, values) <= nocc );
				} else if(cond.op == LT) {
					solver.add( Occurrence(scope, values) <  nocc );
				} else if(cond.op == GE) {
					solver.add( Occurrence(scope, values) >= nocc );
				} else if(cond.op == GT) {
					solver.add( Occurrence(scope, values) >  nocc );
				}  
			}
			
		} else if(cond.operandType == INTERVAL) {
			assert(cond.op == IN);
			
			if(values.size() == 1) {
				solver.add( Occurrence(scope, *begin(values), cond.min, cond.max) );
			} else {
				solver.add( Occurrence(scope, values, cond.min, cond.max) );
			}
			
		} else {
			assert(cond.operandType == INTEGER);
				
			if(values.size() == 1) {
				int val = *begin(values);
				
				if(cond.op == EQ) {
					solver.add( Occurrence(scope, val, cond.val, cond.val) );
				} else if(cond.op == NE) {
					Variable cst(cond.val);
					solver.add( Occurrence(scope, val) != cst);
				} else if(cond.op == LE) {
					solver.add( Occurrence(scope, val, -INFTY, cond.val) );
				} else if(cond.op == LT) {
					solver.add( Occurrence(scope, val, -INFTY, cond.val-1) );
				} else if(cond.op == GE) {
					solver.add( Occurrence(scope, val, cond.val, INFTY) );
				} else if(cond.op == GT) {
					solver.add( Occurrence(scope, val, cond.val+1, INFTY) );
				}
			} else {
				if(cond.op == EQ) {
					solver.add( Occurrence(scope, values, cond.val, cond.val) );
				} else if(cond.op == NE) {
					Variable cst(cond.val);
					solver.add( Occurrence(scope, values) != cst);
				} else if(cond.op == LE) {
					solver.add( Occurrence(scope, values, -INFTY, cond.val) );
				} else if(cond.op == LT) {
					solver.add( Occurrence(scope, values, -INFTY, cond.val-1) );
				} else if(cond.op == GE) {
					solver.add( Occurrence(scope, values, cond.val, INFTY) );
				} else if(cond.op == GT) {
					solver.add( Occurrence(scope, values, cond.val+1, INFTY) );
				}
			}
				
		}
		
}


void XCSP3MistralCallbacks::buildConstraintCount(string id, vector<XVariable *> &list, vector<XVariable *> &values, XCondition &xc) {
    cout << "\n    count constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        values: ";
    displayList(values);
    cout << "        condition: " << xc << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintNValues(string id, vector<XVariable *> &list, vector<int> &except, XCondition &xc) {
    cout << "\n    NValues with exceptions constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        exceptions: ";
    displayList(except);
    cout << "        condition:" << xc << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintNValues(string id, vector<XVariable *> &list, XCondition &xc) {
    cout << "\n    NValues  constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        condition:" << xc << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintCardinality(string id, vector<XVariable *> &list, vector<int> values, vector<int> &occurs, bool closed) {
    cout << "\n    Cardinality constraint (int values, int occurs)  constraint closed: " << closed << endl;
    cout << "        ";
    displayList(list);
    cout << "        values:";
    displayList(values);
    cout << "        occurs:";
    displayList(occurs);
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintCardinality(string id, vector<XVariable *> &list, vector<int> values, vector<XVariable *> &occurs,
                                                     bool closed) {
    cout << "\n    Cardinality constraint (int values, var occurs)  constraint closed: " << closed << endl;
    cout << "        ";
    displayList(list);
    cout << "        values:";
    displayList(values);
    cout << "        occurs:";
    displayList(occurs);
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintCardinality(string id, vector<XVariable *> &list, vector<int> values, vector <XInterval> &occurs,
                                                     bool closed) {
    cout << "\n    Cardinality constraint (int values, interval occurs)  constraint closed: " << closed << endl;
    cout << "        ";
    displayList(list);
    cout << "        values:";
    displayList(values);
    cout << "        occurs:";
    displayList(occurs);
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintCardinality(string id, vector<XVariable *> &list, vector<XVariable *> values, vector<int> &occurs,
                                                     bool closed) {
    cout << "\n    Cardinality constraint (var values, int occurs)  constraint closed: " << closed << endl;
    cout << "        ";
    displayList(list);
    cout << "        values:";
    displayList(values);
    cout << "        occurs:";
    displayList(occurs);
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintCardinality(string id, vector<XVariable *> &list, vector<XVariable *> values, vector<XVariable *> &occurs,
                                                     bool closed) {
    cout << "\n    Cardinality constraint (var values, var occurs)  constraint closed: " << closed << endl;
    cout << "        ";
    displayList(list);
    cout << "        values:";
    displayList(values);
    cout << "        occurs:";
    displayList(occurs);
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintCardinality(string id, vector<XVariable *> &list, vector<XVariable *> values, vector <XInterval> &occurs,
                                                     bool closed) {
    cout << "\n    Cardinality constraint (var values, interval occurs)  constraint closed: " << closed << endl;
    cout << "        ";
    displayList(list);
    cout << "        values:";
    displayList(values);
    cout << "        occurs:";
    displayList(occurs);
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintMinimum(string id, vector<XVariable *> &list, XCondition &xc) {
    cout << "\n    minimum  constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        condition: " << xc << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintMinimum(string id, vector<XVariable *> &list, XVariable *index, int startIndex, RankType rank,
                                                 XCondition &xc) {
    cout << "\n    arg_minimum  constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        index:" << *index << endl;
    cout << "        Start index : " << startIndex << endl;
    cout << "        condition: " << xc << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintMaximum(string id, vector<XVariable *> &list, XCondition &xc) {
    cout << "\n    maximum  constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        condition: " << xc << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintMaximum(string id, vector<XVariable *> &list, XVariable *index, int startIndex, RankType rank, XCondition &xc) {
    cout << "\n    arg_maximum  constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        index:" << *index << endl;
    cout << "        Start index : " << startIndex << endl;
    cout << "        condition: " << xc << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintElement(string id, vector<XVariable *> &list, int value) {
    cout << "\n    element constant constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        value: " << value << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintElement(string id, vector<XVariable *> &list, XVariable *value) {
    cout << "\n    element variable constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        value: " << *value << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintElement(string id, vector<XVariable *> &list, int startIndex, XVariable *index, RankType rank, int value) {
    cout << "\n    element constant (with index) constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        value: " << value << endl;
    cout << "        Start index : " << startIndex << endl;
    cout << "        index : " << *index << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintElement(string id, vector<XVariable *> &list, int startIndex, XVariable *index, RankType rank, XVariable *value) {
    cout << "\n    element variable (with index) constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        value: " << *value << endl;
    cout << "        Start index : " << startIndex << endl;
    cout << "        index : " << *index << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintChannel(string id, vector<XVariable *> &list, int startIndex) {
    cout << "\n    channel constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        Start index : " << startIndex << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintChannel(string id, vector<XVariable *> &list1, int startIndex1, vector<XVariable *> &list2,
                                                 int startIndex2) {
    cout << "\n    channel constraint" << endl;
    cout << "        list1 ";
    displayList(list1);
    cout << "        list2 ";
    displayList(list2);
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintChannel(string id, vector<XVariable *> &list, int startIndex, XVariable *value) {
    cout << "\n    channel constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        value: " << *value << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintStretch(string id, vector<XVariable *> &list, vector<int> &values, vector <XInterval> &widths) {
    cout << "\n    stretch constraint" << endl;
    cout << "        ";
    displayList(list);
    cout << "        values :";
    displayList(values);
    cout << "        widths:";
    displayList(widths);
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintStretch(string id, vector<XVariable *> &list, vector<int> &values, vector <XInterval> &widths, vector <vector<int>> &patterns) {
    cout << "\n    stretch constraint (with patterns)" << endl;
    cout << "        ";
    displayList(list);
    cout << "        values :";
    displayList(values);
    cout << "        widths:";
    displayList(widths);
    cout << "        patterns";
    for(unsigned int i = 0; i < patterns.size(); i++)
        cout << "(" << patterns[i][0] << "," << patterns[i][1] << ") ";
    cout << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintNoOverlap(string id, vector<XVariable *> &origins, vector<int> &lengths, bool zeroIgnored) {
    cout << "\n    nooverlap constraint" << endl;
    cout << "        origins";
    displayList(origins);
    cout << "        lengths";
    displayList(lengths);
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintNoOverlap(string id, vector<XVariable *> &origins, vector<XVariable *> &lengths, bool zeroIgnored) {
    cout << "\n    nooverlap constraint" << endl;
    cout << "        origins:";
    displayList(origins);
    cout << "        lengths";
    displayList(lengths);
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintNoOverlap(string id, vector <vector<XVariable *>> &origins, vector <vector<int>> &lengths, bool zeroIgnored) {
    cout << "\n    kdim (int lengths) nooverlap constraint" << endl;
    cout << "origins: " << endl;
    for(unsigned int i = 0; i < origins.size(); i++) {
        cout << "        ";
        displayList(origins[i]);
    }
    cout << "lengths: " << endl;
    for(unsigned int i = 0; i < origins.size(); i++) {
        cout << "        ";
        displayList(lengths[i]);
    }

		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintNoOverlap(string id, vector <vector<XVariable *>> &origins, vector <vector<XVariable *>> &lengths, bool zeroIgnored) {
    cout << "\n    kdim (lenghts vars nooverlap constraint" << endl;
    cout << "origins: " << endl;
    for(unsigned int i = 0; i < origins.size(); i++) {
        cout << "        ";
        displayList(origins[i]);
    }
    cout << "lengths: " << endl;
    for(unsigned int i = 0; i < origins.size(); i++) {
        cout << "        ";
        displayList(lengths[i]);
    }
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildConstraintInstantiation(string id, vector<XVariable *> &list, vector<int> &values) {
    cout << "\n    instantiation constraint" << endl;
    cout << "        list:";
    displayList(list);
    cout << "        values:";
    displayList(values);

		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildObjectiveMinimizeExpression(string expr) {
    cout << "\n    objective: minimize" << expr << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildObjectiveMaximizeExpression(string expr) {
    cout << "\n    objective: maximize" << expr << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildObjectiveMinimizeVariable(XVariable *x) {
    cout << "\n    objective: minimize variable " << x << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildObjectiveMaximizeVariable(XVariable *x) {
    cout << "\n    objective: maximize variable " << x << endl;
		
		cout << "NOT IMPLEMENTED!" << endl;
		exit(1);
}


void XCSP3MistralCallbacks::buildObjectiveMinimize(ExpressionObjective type, vector<XVariable *> &list, vector<int> &coefs) {
    XCSP3CoreCallbacks::buildObjectiveMinimize(type, list, coefs);
}


void XCSP3MistralCallbacks::buildObjectiveMaximize(ExpressionObjective type, vector<XVariable *> &list, vector<int> &coefs) {
    XCSP3CoreCallbacks::buildObjectiveMaximize(type, list, coefs);
}


void XCSP3MistralCallbacks::buildObjectiveMinimize(ExpressionObjective type, vector<XVariable *> &list) {
    XCSP3CoreCallbacks::buildObjectiveMinimize(type, list);
}


void XCSP3MistralCallbacks::buildObjectiveMaximize(ExpressionObjective type, vector<XVariable *> &list) {
    XCSP3CoreCallbacks::buildObjectiveMaximize(type, list);
}

Variable XCSP3MistralCallbacks::postExpression(Node *n, bool root) {
	Variable rv;
	
    if(n->type == NT_VARIABLE) {
        assert(!root);
				
        NodeVariable *nv = (NodeVariable *) n;
        return variable[nv->var];
				
    }

    if(n->type == NT_CONSTANT) {
        assert(!root);
				
        NodeConstant *nc = (NodeConstant *) n;
        return Variable(nc->val,nc->val);

    }

    NodeOperator *fn = (NodeOperator *) n;

    if(fn->type == NT_EQ) {
			
        Variable x1 = postExpression(fn->args[0]);
        Variable x2 = postExpression(fn->args[1]);
				rv = (x1 == x2);
   
    }

    if(fn->type == NT_NE) {
			
        Variable x1 = postExpression(fn->args[0]);
        Variable x2 = postExpression(fn->args[1]);
        rv = (x1 != x2);
				
    }

    if(fn->type == NT_GE) {
        
				Variable x1 = postExpression(fn->args[0]);
        Variable x2 = postExpression(fn->args[1]);
				rv = ( x1 >= x2 );

    }

    if(fn->type == NT_GT) {
        
				Variable x1 = postExpression(fn->args[0]);
        Variable x2 = postExpression(fn->args[1]);
				rv = ( x1 > x2 );

    }


    if(fn->type == NT_LE) {
			
        Variable x1 = postExpression(fn->args[0]);
        Variable x2 = postExpression(fn->args[1]);
 				rv = ( x1 <= x2 );
				
    }

    if(fn->type == NT_LT) {
        
				Variable x1 = postExpression(fn->args[0]);
        Variable x2 = postExpression(fn->args[1]);
				rv = ( x1 < x2 );
    
		}

    if(fn->type == NT_IMP) { // IMP(X,Y) = NOT X OR Y
			
				Variable x1 = postExpression(fn->args[0]);
      	Variable x2 = postExpression(fn->args[1]);
				rv = ( x1 <= x2 );
			
    }

    if(fn->type == NT_OR) {
			
				Variable x1 = postExpression(fn->args[0]);
    		Variable x2 = postExpression(fn->args[1]);
				rv = ( x1 || x2 );

    }

    if(fn->type == NT_AND) {
			
				Variable x1 = postExpression(fn->args[0]);
  			Variable x2 = postExpression(fn->args[1]);
				rv = ( x1 && x2 );
    
		}

    if(fn->type == NT_NOT) {
			
				Variable x1 = postExpression(fn->args[0]);
				rv = ( !x1 );
			
    }
		
    if(fn->type == NT_IFF) {
			
				Variable x1 = postExpression(fn->args[0]);
				Variable x2 = postExpression(fn->args[1]);
				rv = ( x1 == x2 );
			
    }

    if(fn->type == NT_XOR) {
			
				Variable x1 = postExpression(fn->args[0]);
				Variable x2 = postExpression(fn->args[1]);
				rv = ( x1 != x2 );
			
    }

    // function stuff
    if(fn->type == NT_NEG) {
			
				Variable x1 = postExpression(fn->args[0]);
				rv = ( -x1 );
				
    }

    if(fn->type == NT_ABS) {

				Variable x1 = postExpression(fn->args[0]);
				rv = ( Abs(x1) );

    }

    if(fn->type == NT_SUB) {
			
				Variable x1 = postExpression(fn->args[0]);
				Variable x2 = postExpression(fn->args[1]);
				rv = ( x1 - x2 );
			
    }

    if(fn->type == NT_DIST) { //Simulate DIST(X,Y) = ABS(SUB(X,Y))
			
				Variable x1 = postExpression(fn->args[0]);
				Variable x2 = postExpression(fn->args[1]);
				rv = ( Abs(x1 - x2) );
			
    }

    if(fn->type == NT_ADD) {
			
				Variable x1 = postExpression(fn->args[0]);
				Variable x2 = postExpression(fn->args[1]);
				rv = ( x1 + x2 );
			
    }
    if(fn->type == NT_MULT) {
			
				Variable x1 = postExpression(fn->args[0]);
				Variable x2 = postExpression(fn->args[1]);
				rv = ( x1 * x2 );
			
    }

    if(fn->type == NT_MIN) {
				assert(!root);
			
				Variable x1 = postExpression(fn->args[0]);
				Variable x2 = postExpression(fn->args[1]);
				rv = Min( x1, x2 );
			
    }
    if(fn->type == NT_MAX) {
        assert(!root);
			
				Variable x1 = postExpression(fn->args[0]);
				Variable x2 = postExpression(fn->args[1]);
				rv = Min( x1, x2 );
			
    }

    if(fn->type == NT_IF) {
	      assert(!root);
		
				Variable x = postExpression(fn->args[0]);
				VarArray A;
				A.add(postExpression(fn->args[2]));
				A.add(postExpression(fn->args[1]));
				rv = Element(A,x);

    }

    return rv;
}


#endif //COSOCO_XCSP3PRINTCALLBACKS_H
