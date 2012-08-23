/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Guido Tack <tack@gecode.org>
 *
 *  Copyright:
 *     Guido Tack, 2007
 *
 *  Last modified:
 *     $Date: 2010-05-11 12:33:38 +0200 (Tue, 11 May 2010) $ by $Author: tack $
 *     $Revision: 10940 $
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

#include "flatzinc.hpp"
#include "registry.hpp"
#include <iomanip>

#include <vector>
#include <assert.h>
#include <string>
#include <set>



#include <mistral_variable.hpp>

//#define _DEBUG_FLATZINC true
//#define _VERBOSE_PARSER 1



using namespace std;

namespace FlatZinc {
inline
set<int> setrange(int min, int max)
{
	set<int> rv;
	for(int i = min; i <= max; ++i)
		rv.insert(i);
	return rv;
}

set<int> vs2is(IntVarSpec* vs) {
	if (vs->assigned) {
		return setrange(vs->i,vs->i);
	}
	if (vs->domain()) {
		AST::SetLit* sl = vs->domain.some();
		if (sl->interval) {
			return setrange(sl->min, sl->max);
		} else {
			set<int> rv;
			for (int i=sl->s.size(); i--;)
				rv.insert(sl->s[i]);
			return rv;
		}
	}
	return setrange(-1000, 1000);
}


Variable vs2var(IntVarSpec* vs) {
	Variable x;
	if (vs->assigned) {
		Variable y(vs->i);
		x = y;
	} else if (vs->domain()) {
		AST::SetLit* sl = vs->domain.some();
		if (sl->interval) {
			Variable y(sl->min, sl->max);
			x = y;
		} else {
			set<int> __tmp_set;
			int size = sl->s.size();
			for (int i=0; i<size;i++)
				__tmp_set.insert(sl->s[i]);
			//	cout << "myset contains:";
			//	for (set<int>::iterator  it=__tmp_set.begin(); it!=__tmp_set.end(); it++)
			//	cout << " " << *it;
			//	cout << endl;
			Vector<int> values;
			for (set<int>::iterator it=__tmp_set.begin(); it!=__tmp_set.end(); it++)
				values.add(*it);

			Variable y(values);
			x = y;
		}
	} else {
		Variable y(-INFTY/1024, +INFTY/1024);
		x = y;
	}

	//std::cout << x << " in " << x.get_domain() << std::endl;

	return x;
}

int vs2bsl(BoolVarSpec* bs) {
	if (bs->assigned) {
		return bs->i;
	}
	if (bs->domain()) {
		AST::SetLit* sl = bs->domain.some();
		assert(sl->interval);
		return std::min(1, std::max(0, sl->min));
	}
	return 0;
}

int vs2bsh(BoolVarSpec* bs) {
	if (bs->assigned) {
		return bs->i;
	}
	if (bs->domain()) {
		AST::SetLit* sl = bs->domain.some();
		assert(sl->interval);
		return std::max(0, std::min(1, sl->max));
	}
	return 1;
}

int ann2ivarsel(AST::Node* ann) {
	if (/*AST::Atom* s =*/ dynamic_cast<AST::Atom*>(ann)) {
		// if (s->id == "input_order")
		//   return TieBreakVarBranch<IntVarBranch>(INT_VAR_NONE);
	}
	std::cerr << "Warning, ignored search annotation: ";
	ann->print(std::cerr);
	std::cerr << std::endl;
	return 0;
}

int ann2ivalsel(AST::Node* ann) {
	if (/*AST::Atom* s =*/ dynamic_cast<AST::Atom*>(ann)) {
		// if (s->id == "indomain_min")
		//   return INT_VAL_MIN;
	}
	std::cerr << "Warning, ignored search annotation: ";
	ann->print(std::cerr);
	std::cerr << std::endl;
	return 0;
}

int ann2asnivalsel(AST::Node* ann) {
	if (/*AST::Atom* s =*/ dynamic_cast<AST::Atom*>(ann)) {
		// if (s->id == "indomain_min")
		//   return INT_ASSIGN_MIN;
	}
	std::cerr << "Warning, ignored search annotation: ";
	ann->print(std::cerr);
	std::cerr << std::endl;
	return 0;
}


FlatZincModel::FlatZincModel(Solver &s)
: solver(s), heuristic(NULL), policy(NULL), use_rewriting(false),
  intVarCount(-1), boolVarCount(-1), setVarCount(-1), _optVar(-1),
  _solveAnnotations(NULL),
  findall(false)
{}

void
FlatZincModel::init(int intVars, int boolVars, int setVars) {
	intVarCount = 0;
	iv = IntVarArray(intVars);
	iv_introduced = std::vector<bool>(intVars);
	iv_boolalias = std::vector<int>(intVars);
	boolVarCount = 0;
	bv = BoolVarArray(boolVars);
	bv_introduced = std::vector<bool>(boolVars);
	setVarCount = 0;
	sv = SetVarArray(setVars);
	sv_introduced = std::vector<bool>(setVars);
}

void
FlatZincModel::newIntVar(IntVarSpec* vs) {

#ifdef _VERBOSE_PARSER
	if((intVarCount % _VERBOSE_PARSER) == 0) {
		std::cout << "i";
		std::cout.flush();
	}
#endif

	if (vs->alias) {
		iv[intVarCount++] = iv[vs->i];
	} else {
		Variable x = vs2var(vs);
		iv[intVarCount++] = x;
	}
	iv_introduced[intVarCount-1] = vs->introduced;
	iv_boolalias[intVarCount-1] = -1;
}

void
FlatZincModel::newSetVar(SetVarSpec* vs) {

#ifdef _VERBOSE_PARSER
	if((setVarCount % _VERBOSE_PARSER) == 0) {
		std::cout << "s";
		std::cout.flush();
	}
#endif

	if (vs->alias) {
		sv[intVarCount++] = sv[vs->i];
	} else if( vs->assigned) {
		assert(vs->upperBound());
		AST::SetLit* vsv = vs->upperBound.some();
		if (vsv->interval) {
			Variable x = SetVariable(vsv->min, vsv->max);
			sv[setVarCount++] = x;
			//for(int i = vsv->min; i <= vsv->max; ++i)
			//x.include(solver, i, NO_REASON);
		} else {
			if( vsv->s.empty() ) {
				Variable x = SetVariable(0, 0);
				//x.exclude(solver, 0, NO_REASON);
				//sv[setVarCount++] = x;
			} else {
				int umin = vsv->s[0], umax = vsv->s[0];
				for(size_t i = 1; i != vsv->s.size(); ++i) {
					umin = std::min(umin, vsv->s[i]);
					umax = std::max(umax, vsv->s[i]);
				}
				Variable x = SetVariable(umin, umax);
				sv[setVarCount++] = x;
				// for(size_t i = 0; i != vsv->s.size(); ++i)
				//   x.include(solver, vsv->s[i], NO_REASON);
				// for(int i = x.umin(solver), iend = x.umax(solver); i <= iend; ++i)
				//   if( !x.includes(solver, i) )
				//     x.exclude(solver, i, NO_REASON);
			}
		}
	} else if( vs->upperBound() ) {
		AST::SetLit* vsv = vs->upperBound.some();
		Variable x = SetVariable(vsv->min, vsv->max);
		//setvar x = solver.newSetVar(vsv->min, vsv->max);
		sv[setVarCount++] = x;
		// if( !vsv->interval ) {
		//   int prev = vsv->min;
		//   for(size_t i = 0; i != vsv->s.size(); ++i) {
		//     if( vsv->s[i] > prev+1 ) {
		//       for(int q = prev+1; q != vsv->s[i]; ++q)
		//         x.exclude(solver, q, NO_REASON);
		//     }
		//     prev = vsv->s[i];
		//   }
		// } // otherwise everything is unset and we are done here
	} else {
		// completely free
		//setvar x = solver.newSetVar(-1000, 1000);
		Variable x = SetVariable(-1000, 1000);
		sv[setVarCount++] = x;
	}
	sv_introduced[setVarCount-1] = vs->introduced;
}

void
FlatZincModel::aliasBool2Int(int iv, int bv) {
	iv_boolalias[iv] = bv;
}
int
FlatZincModel::aliasBool2Int(int iv) {
	return iv_boolalias[iv];
}

void
FlatZincModel::newBoolVar(BoolVarSpec* vs) {

#ifdef _VERBOSE_PARSER
	if((boolVarCount % _VERBOSE_PARSER) == 0) {
		std::cout << "b";
		std::cout.flush();
	}
#endif

	if (vs->alias) {
		bv[boolVarCount++] = bv[vs->i];
	} else {
		bv[boolVarCount++] = Variable(vs2bsl(vs), vs2bsh(vs));
	}
	bv_introduced[boolVarCount-1] = vs->introduced;
}


#ifdef _VERIFICATION
SolutionValue
FlatZincModel::node2SolutionValue(AST::Node * ai ) {
	SolutionValue v;
	int k;
	if (ai->isBool())
	{
		oss.str("");
		//		cout<<"isBool ";
		oss << (ai->getBool() ? "true" : "false");
		v.set_type(1);
		v.set_val(oss.str());
	}
	else
		if (ai->isString())
		{
			oss.str("");
			v.set_type(2);
			//	cout<<" Test if node is a isString node";
			std::string s = ai->getString();
			for (unsigned int i=0; i<s.size(); i++)
			{
				if (s[i] == '\\' && i<s.size()-1)
				{
					switch (s[i+1])
					{
					case 'n': oss << "\n"; break;
					case '\\': oss << "\\"; break;
					case 't': oss << "\t"; break;
					default: oss << "\\" << s[i+1];
					}
					i++;
				}
				else
				{
					oss << s[i];
				}
			}

			v.set_val(oss.str());
		}
		else
			if (ai->isInt(k))
			{
				oss.str("");
				oss << k;
				v.set_type(3);
				v.set_val(oss.str());
			}
			else
				if (ai->isAtom())
				{
					oss.str("");
					//cout<<" Test if node is isAtom node";
					v.set_type(4);
					v.set_val("ATOM is not yet supported");
				}
				else
					if (ai->isSet())
					{
						//		cout<<" Test if node is isSet  node";
						AST::SetLit* setLit = new AST::SetLit();
						setLit->interval =ai->getSet()->interval;
						setLit->min =ai->getSet()->min;
						setLit->max =ai->getSet()->max;
						setLit->s =ai->getSet()->s;
						v.set_ai(setLit);
						v.set_type(5);
					}
					else
						if(	ai->isBoolVar() )
						{
							//		cout<<"isBoolVar ";
							v.set_type(6);
							Variable *_var =  & (bv[ai->getBoolVar()]) ;
							v.set_var( _var);
						}
						else
							if (ai->isArray())
							{
								std::vector<SolutionValue> __tmp;
								//	cout<<" Test if node is isArray node" << endl ;
								for (unsigned int i =0; i< ai->getArray()->a.size(); i++)
								{
									AST::Node *n = ai->getArray()->a[i];
									__tmp.push_back(node2SolutionValue(n));
									v.set_a(__tmp);

								}
								v.set_type(9);
							}
							else
								if (ai->isSetVar())
								{
									//		cout<<" Test if node is a isSetVar node";
									Variable *_var =  & (sv[ai->getSetVar()]);
									v.set_var( _var);
									v.set_type(8);
								}
								else
									if (ai->isIntVar())
									{
										//	cout<<"isIntVar ";
										Variable *_var = &( iv[ai->getIntVar()]) ;
										//std::cout << _var.get_domain() << std::endl;
										v.set_var(_var);
										v.set_type(7);
									}
	return v;
}
#endif

void
FlatZincModel::postConstraint(const ConExpr& ce, AST::Node* ann) {
	try {
#ifdef _VERBOSE_PARSER
		if((solver.constraints.size % _VERBOSE_PARSER) == 0) {
			std::cout << "c";
			std::cout.flush();
		}
#endif
		registry().post(solver, *this, ce, ann);

#ifdef _VERIFICATION
		std::vector<SolutionValue > __vars;
		//	std::cout <<" \n posting " << ce.id <<std::endl ;
		for (unsigned int i=0; i<ce.args->a.size(); i++)
		{
			AST::Node * ai = ce.args->a[i];
			__vars.push_back(node2SolutionValue(ai));
		}
		pair<std::string, std::vector<SolutionValue > > pair (ce.id , __vars);
		verif_constraints.push_back(pair);
#endif

	} catch (AST::TypeError& e) {
		throw FlatZinc::Error("Type error", e.what());
	}

}





void flattenAnnotations(AST::Array* ann, std::vector<AST::Node*>& out) {
	for (unsigned int i=0; i<ann->a.size(); i++) {
		if (ann->a[i]->isCall("seq_search")) {
			AST::Call* c = ann->a[i]->getCall();
			if (c->args->isArray())
				flattenAnnotations(c->args->getArray(), out);
			else
				out.push_back(c->args);
		} else {
			out.push_back(ann->a[i]);
		}
	}
}


bool FlatZincModel::getAnnotations( AST::Call*  c , Vector<Variable> &__vars, std::string & _varHeuristic , std::string & _valHeuristic)
{

	if (c->isCall("int_search") || c->isCall("bool_search"))
	{
	if (c->args->isArray())
	{
	//	cout<< "c array size : "<< c->args->getArray()->a.size() ;
		if (c->args->getArray()->a[0]->isArray())
		{
			AST::Array * __varsArray = c->args->getArray()->a[0]->getArray();
			//cout<< " c number of branching variables : "<< __varsArray->a.size() << endl;

			for (unsigned int j=0; j< __varsArray->a.size(); j++ )
			{
				if (__varsArray->a[j]->isIntVar())
					__vars.push_back(iv[__varsArray->a[j]->getIntVar()]);
				else if (__varsArray->a[j]->isSetVar())
					__vars.push_back( sv[__varsArray->a[j]->getSetVar()]);
				else if (__varsArray->a[j]->isBoolVar())
					__vars.push_back( bv[__varsArray->a[j]->getBoolVar()]);
				else
					return false;
			}
		}

		else return false;

		//test if atom ?
		std::ostringstream __tmpStream;
		c->args->getArray()->a[1]->print(__tmpStream);
		_varHeuristic = __tmpStream.str();
		__tmpStream.str("");
		c->args->getArray()->a[2]->print(__tmpStream);
		_valHeuristic = __tmpStream.str();

	}
	else
		return false;
	}
	else
		return false;



	return true;
}


void
FlatZincModel::createBranchers(AST::Node* ann, bool ignoreUnknown,
		std::ostream& err) {
	if (ann) {
		err << "Warning, ignored search annotation: ";
		ann->print(err);
		err << std::endl;
	}
}

AST::Array*
FlatZincModel::solveAnnotations(void) const {
	return _solveAnnotations;
}

void
FlatZincModel::solve(AST::Array* ann) {
	_method = SATISFACTION;
	_solveAnnotations = ann;
}

// void
// FlatZincModel::enumerate(AST::Array* ann) {
// 	_method = ENUMERATION;
// 	_solveAnnotations = ann;
// 	// Branch on optimization variable to ensure that it is given a value.
// 	/*
// 	AST::Array* args = new AST::Array(4);
// 	args->a[0] = new AST::Array(new AST::IntVar(_optVar));
// 	args->a[1] = new AST::Atom("input_order");
// 	args->a[2] = new AST::Atom("indomain_min");
// 	args->a[3] = new AST::Atom("complete");
// 	AST::Call* c = new AST::Call("int_search", args);
// 	if (!ann)
// 		ann = new AST::Array(c);
// 	else
// 		ann->a.push_back(c);
// 	*/

// }

void
FlatZincModel::minimize(int var, AST::Array* ann) {
	_method = MINIMIZATION;
	_optVar = var;
	_solveAnnotations = ann;
	// Branch on optimization variable to ensure that it is given a value.
	/*
	AST::Array* args = new AST::Array(4);
	args->a[0] = new AST::Array(new AST::IntVar(_optVar));
	args->a[1] = new AST::Atom("input_order");
	args->a[2] = new AST::Atom("indomain_min");
	args->a[3] = new AST::Atom("complete");
	AST::Call* c = new AST::Call("int_search", args);
	if (!ann)
		ann = new AST::Array(c);
	else
		ann->a.push_back(c);
	*/

}

void
FlatZincModel::maximize(int var, AST::Array* ann) {
	_method = MAXIMIZATION;
	_optVar = var;
	_solveAnnotations = ann;
	// Branch on optimization variable to ensure that it is given a value.
	/*
	AST::Array* args = new AST::Array(4);
	args->a[0] = new AST::Array(new AST::IntVar(_optVar));
	args->a[1] = new AST::Atom("input_order");
	args->a[2] = new AST::Atom("indomain_min");
	args->a[3] = new AST::Atom("complete");
	AST::Call* c = new AST::Call("int_search", args);
	if (!ann)
		ann = new AST::Array(c);
	else
		ann->a.push_back(c);
	 */
}

FlatZincModel::~FlatZincModel(void) {
	delete _solveAnnotations;

        // for(unsigned int i=0; i<iv.size; ++i) {
        
        //   int domain_type = iv[i].domain_type;
        //   if     (domain_type ==  BITSET_VAR) delete iv[i].bitset_domain;
        //   else if(domain_type ==    LIST_VAR) delete iv[i].list_domain;
        //   else if(domain_type ==   RANGE_VAR) delete iv[i].range_domain;
        //   else if(domain_type == VIRTUAL_VAR) delete iv[i].virtual_domain;
        //   else if(domain_type ==  EXPRESSION) delete iv[i].expression;
        //   else if(domain_type !=   CONST_VAR) delete iv[i].variable;
 
        // }
}



void
FlatZincModel::set_parameters(SolverParameters& p) {

}


void
FlatZincModel::set_strategy(string var_o, string val_o, string r_pol) {
	heuristic = solver.heuristic_factory(var_o, val_o);
	policy = solver.restart_factory(r_pol);
}

void
FlatZincModel::set_rewriting(const bool on) {
	use_rewriting = on;
}

void 
FlatZincModel::set_enumeration(const bool on) {
  enumerate = on;
}

  void
  FlatZincModel::run(std::ostream& out, const Printer& p) {
    using std::setw;
    using std::setfill;
    

#ifdef _DEBUG_FLATZINC
    std::cout << " c run!" << std::endl;
#endif

    if(use_rewriting) {
#ifdef _DEBUG_FLATZINC
      std::cout << "before rewriting:\n" << solver << std::endl;
#endif

      solver.rewrite() ;

#ifdef _DEBUG_FLATZINC
      std::cout << "after rewriting:\n" << solver << std::endl;
#endif
    }

    solver.consolidate();

#ifdef _DEBUG_FLATZINC
    std::cout << "c mistral representation:\n " << solver << std::endl;
#endif
        
    Outcome result = UNKNOWN;


#ifdef _FLATZINC_OUTPUT
    cout << "%";
#endif

    cout << " c search annotations are not yet supported." << endl;

    /*search annotations :
     * __search_strategies: the number of search strategies.
     * __search_strategies is equal to 0 when no specific search strategy has been specified. With bool_search or int_search, it's egual to 1 and >1 with seq_search.
     int __search_strategies = 0;
     Vector<Variable> banching_variables;
     std::string var_heuristic;
     std::string val_heuristic;

     if (_solveAnnotations!= NULL)
     {
     __search_strategies = 1;
     cout << " c annotations size : " << _solveAnnotations->a.size() << endl;

     if (_solveAnnotations->a[0]->isCall("seq_search")) {
     cout << " c SEQ_SEARCH " << endl ;

     AST::Call* c = _solveAnnotations->a[0]->getCall();
     if (c->args->isArray())
     {
     __search_strategies =c->args->getArray()->a.size();
     cout << " c Total number of search strategies:" << __search_strategies  << endl;

     for (unsigned j = 0; j< __search_strategies; j++)
     {
     Vector<Variable> __banching_variables;
     std::string __var_heuristic;
     std::string __val_heuristic;

     if (getAnnotations(c->args->getArray()->a[j]->getCall(),__banching_variables, __var_heuristic ,__val_heuristic))
     {
     cout << " c Search strategy number "<< j+1 << endl;

     cout << " c var heuristic : "<< __var_heuristic << endl;
     cout << " c val heuristic : "<< __val_heuristic << endl;
     cout << " c _variables size : "<< __banching_variables.size  << endl;
     //	cout << " c branching variables : \n "<< __banching_variables  << endl;

     }
     else
     cout << " c Something wrong with search annotations. The solver will use the default search strategy." << endl;
     }
     cout << " c seq_search is not yet supported." << endl;
     }
     else
     cout << " c Something wrong with search annotations. The solver will use the default search strategy." << endl;
     }
     else
     {
     cout << " c Total number of search strategies:" << __search_strategies  << endl;
     if (getAnnotations(_solveAnnotations->getArray()->a[0]->getCall(), banching_variables, var_heuristic ,val_heuristic ))
     {
     cout << " c var heuristic : "<< var_heuristic  << endl;
     cout << " c val heuristic : "<< val_heuristic  << endl;
     cout << " c number of banching variables : "<< banching_variables.size  << endl;
     //cout << " c branching variables : \n "<< banching_variables  << endl;

     }
     else
     cout << " c Something wrong with search annotations. The solver will use the default search strategy." << endl;
     }
     }
     else
     cout << " c No specific annotation. The solver will use the default search strategy." << endl;
    */

    switch (_method) {
    case MINIMIZATION: {

      //#ifdef _DEBUG_FLATZINC
#ifdef _FLATZINC_OUTPUT
      cout << "%";
#endif
      std::cout << " c Minimize " << iv[_optVar].get_var() << std::endl;
      //#endif

      Goal *goal = new Goal(Goal::MINIMIZATION, iv[_optVar].get_var());
      //	if (__search_strategies == 1)
      //		result = solver.depth_first_search(banching_variables, heuristic, policy, goal);
      //	else
      result = solver.depth_first_search(solver.variables, heuristic, policy, goal);

      //result = solver.minimize(iv[_optVar]);
      break;
    }
    case MAXIMIZATION: {

      //#ifdef _DEBUG_FLATZINC
#ifdef _FLATZINC_OUTPUT
      cout << "%";
#endif
      std::cout << " c Maximize " << iv[_optVar].get_var() << std::endl;
      //#endif

      Goal *goal = new Goal(Goal::MAXIMIZATION, iv[_optVar].get_var());
      //	if (__search_strategies == 1)
      //		result = solver.depth_first_search(banching_variables, heuristic, policy, goal);
      //	else
      result = solver.depth_first_search(solver.variables, heuristic, policy, goal);
      //result = solver.maximize(iv[_optVar]);
      break;
    }
    case SATISFACTION: {

      //#ifdef _DEBUG_FLATZINC
#ifdef _FLATZINC_OUTPUT
      cout << "%";
#endif
      std::cout << " c Solve " << std::endl;
      //#endif

      //Goal *goal = new Goal(Goal::SATISFACTION);


      // for(int i=0; i<solver.variables.size; ++i) {
      //   solver.monitor_list << solver.variables[i] ;
      //   solver.monitor_list << " " ;
      // }


      // int matrix[100];


      // //std::cout << iv << std::endl;
      // for(int i=0; i<100; ++i) {
      //   Variable x = iv[i].get_var();
      //   if(x.domain_type == CONST_VAR)
      //     matrix[i] = x.get_value();
      //   else
      //     matrix[i] = -(x.get_id());
      // }

#ifdef _MONITOR
      for(int i=0; i<10; ++i) {
        for(int j=0; j<10; ++j) {
          solver.monitor_list << " ";
          Variable x = iv[i*10+j].get_var();
          if(x.domain_type == CONST_VAR) {
            //solver.monitor_list << x.get_value();
            solver.monitor_list << " ";
          } else {
            solver.monitor_list << solver.variables[x.id()];
          }
        }
        solver.monitor_list << "\n";
      }
#endif
      
      Goal *goal;
      if(enumerate) 
        goal = new Goal(Goal::ENUMERATION);
      else 
        goal = new Goal(Goal::SATISFACTION);

      //if (__search_strategies == 1)
      //		result = solver.depth_first_search(banching_variables, heuristic, policy, goal);
      //	else
      result = solver.depth_first_search(solver.variables, heuristic, policy, goal);


      // //std::cout << iv << std::endl;
      // for(int i=0; i<10; ++i) {
      //   for(int j=0; j<10; ++j) {
      //     Variable x = iv[i*10+j].get_var();
      //     if(x.is_ground())
      //       std::cout << setw(2) << x.get_value() << "     ";
      //     else
      //       std::cout << setw(2) << x.get_solution_int_value() << " (" << x.id() << ") ";

      //       //std::cout << x << " in " << x.get_domain() << " ";
      //   }
      //   std::cout << std::endl;
      // }
      // //result = solver.solve();
      break;
    }
    }


#ifdef _DEBUG_VERIFICATION
    if ((result == SAT) || (result == OPT) || ((result == LIMITOUT) && (solver.statistics.num_solutions >0)))
      {
        cout<< " \n \n \n solver variables :  " << solver.variables.size ;
        int __size = solver.variables.size;
        for ( int i = 0; i < __size ; i++)
          {
            cout<< "  " << endl;
            std::cout << solver.variables[i].get_solution_int_value() << " (" << solver.variables[i].id() << ") ";
            std::cout << solver.variables[i] << " in " << solver.variables[i].get_domain() << " ";
          }
      }
#endif


    //std::cout << solver.statistics << std::endl;


    // switch(result) {
    // case UNKNOWN: {
    //   out << "c" << setw(5) << setfill('=') << '='
    //       << "UNKNOWN" << setw(5) << '=' << "\n";
    //   break;
    // }
    // case SAT: {
    //   print(out, p);
    //   out << "c" << setw(5) << setfill('=') << '='
    //       << "SAT" << setw(5) << '=' << "\n";
    //   break;
    // }
    // case UNSAT: {
    //   out << "c" << setw(5) << setfill('=') << '='
    //       << "UNSAT" << setw(5) << '=' << "\n";
    //   break;
    // }
    // case OPT: {
    //   print(out, p);
    //   out << "c" << setw(5) << setfill('=') << '='
    //       << "OPTIMAL" << setw(5) << '=' << "\n";
    //   break;
    // }
    // }


    // if(solver.statistics.num_solutions) {
    //   for(unsigned int i=0; i<iv.size; ++i) {
    //     out << iv[i].get_var() << " = " << iv[i].get_solution_str_value() << " ";
    //   }
    //   out << endl;
    // }

  }


FlatZincModel::Meth
FlatZincModel::method(void) const {
	return _method;
}

int
FlatZincModel::optVar(void) const {
	return _optVar;
}

void
FlatZincModel::print_final(std::ostream& out, const Printer& p) const {

#ifdef _FLATZINC_OUTPUT
	 // std::cout << "% c +" << std::setw(90) << std::setfill('=')
	 //  			      //=============================================================================
	 //  				      << "+" << std::endl << std::setfill(' ')
	 //  				      << "% c |      INSTANCE STATS       |                    SEARCH STATS                 | OBJECTIVE |" << std::endl
	 //  				      << "% c |   vars |    vals |   cons |    nodes | filterings | propagations | cpu time |           |" << std::endl;


	 // std::cout << "% c +" << std::setw(90) << std::setfill('=')
	 //    //"=============================================================================
	 //     << "+" << std::endl << std::setfill(' ')
	 //     << std::left << std::setw(46) << "% s  ";

	 //  switch(solver.statistics.outcome) {
	 //  case SAT:
	 //          std::cout << std::right << std::setw(47) << "SATISFIABLE" ;
	 //    break;
	 //  case OPT:
	 //          std::cout << std::right << std::setw(47) << "OPTIMAL" ;
	 //    break;
	 //  case UNSAT:
	 //          std::cout << std::right << std::setw(47) << "UNSATISFIABLE" ;
	 //    break;
	 //  case UNKNOWN:
	 //          std::cout << std::right << std::setw(47) << "UNKNOWN" ;
	 //    break;
	 //  case LIMITOUT:
	 //    if(solver.statistics.num_solutions > 0)
	 //    	 std::cout << std::right << std::setw(47) << "SUBOPTIMAL" ;
	 //    else
	 //    	 std::cout << std::right << std::setw(47) << "LIMITOUT" ;
	 //    //break;
	 //  }
	 //  std::cout << std::endl
	 //     << std::left << std::setw(46) << "% v  0" << std::endl
	 //     << std::left << std::setw(46) << "% d  OBJECTIVE"
	 //     << std::right << std::setw(46) << solver.statistics.objective_value  << std::endl
	 //     << std::left << std::setw(46) << "% d  TIME"
	 //     << std::right << std::setw(46) << (solver.statistics.end_time - solver.statistics.start_time)  << std::endl
	 //     << std::left << std::setw(46) << "% d  MEMORY"
	 //     << std::right << std::setw(46) << (Mistral::mem_used() / 1048576.0) << std::endl
	 //     << std::left << std::setw(46) << "% d  NODES"
	 //     << std::right << std::setw(46) << solver.statistics.num_nodes  << std::endl
	 //     << std::left << std::setw(46) << "% d  RESTARTS"
	 //     << std::right << std::setw(46) << solver.statistics.num_restarts << std::endl
	 //     << std::left << std::setw(46) << "% d  FAILURES"
	 //     << std::right << std::setw(46) << solver.statistics.num_failures << std::endl
	 //     << std::left << std::setw(46) << "% d  BACKTRACKS"
	 //     << std::right << std::setw(46) << solver.statistics.num_backtracks << std::endl
	 //     << std::left << std::setw(46) << "% d  PROPAGATIONS"
	 //     << std::right << std::setw(46) << solver.statistics.num_propagations << std::endl
	 //     << std::left << std::setw(46) << "% d  FILTERINGS"
	 //     << std::right << std::setw(46) << solver.statistics.num_filterings << std::endl
	 //     << "% c +" << std::setw(90) << std::setfill('=') << "+" << std::endl << std::setfill(' ');
	 //  //<< " c +=============================================================================+" << std::endl;

          //p.print(out, solver, iv, bv, sv);
	Mistral::Outcome outcome = solver.statistics.outcome;
			//

        //std::cout << "%% " << outcome2str(outcome) << std::endl;;
        
        if (outcome == OPT)
          out<<"==========";
        else if (solver.statistics.num_solutions && solver.objective && solver.objective->is_optimization())
          out<<"=====UNBOUNDED=====";
	else if (outcome == UNSAT)
          out<<"=====UNSATISFIABLE=====";
        else if (outcome != SAT)
          out<<"=====UNKNOWN=====";

	/* Two missing details:
	 1-		We need to print "==========" if we already explored all the search space and we found at least one solution.
	 2- 	Also we need to check whether the objective of an optimization problem is unbounded. In this case, we should print "=====UNBOUNDED=====";
	 */

	out<< std::endl;
#endif


}


void
FlatZincModel::print_solution(std::ostream& out, const Printer& p) const {

#ifdef _FLATZINC_OUTPUT

  if(solver.statistics.num_solutions) {
    p.print(out, solver, iv, bv, sv);
  }
  out << "----------" << std::endl;

#endif


}



void
Printer::init(AST::Array* output) {
	_output = output;
}

void
Printer::printElem(std::ostream& out,
		Solver& solver,
		AST::Node* ai,
		const IntVarArray& iv,
		const BoolVarArray& bv,
		const SetVarArray& sv
) const {
	int k;
	if (ai->isInt(k)) {
		out << k;
	} else if (ai->isIntVar()) {
/*		int lb = iv[ai->getIntVar()].get_solution_min();
		int ub = iv[ai->getIntVar()].get_solution_max();
		if( lb == ub )
			out << lb;
		else
			out << lb << ".." << ub;
*/
          out << iv[ai->getIntVar()].get_solution_int_value();
	} else if (ai->isBoolVar()) {
		/*
		int lb = bv[ai->getBoolVar()].get_solution_min();
		int ub = bv[ai->getBoolVar()].get_solution_max();
		if (lb == 1) {
			out << "true";
		} else if (ub == 0) {
			out << "false";
		} else {
			out << "false..true";
		}
		*/
		int lb = bv[ai->getBoolVar()].get_solution_int_value();
		if (lb == 1)
		{
			out << "true";
		}
		else
			out << "false";
	} else if( ai->isSetVar()) {
		SetExpression *x = (SetExpression*)(sv[ai->getSetVar()].expression);
		set<int> lb;
		set<int> ub;
		for(unsigned int i=0; i<x->elts_ub.size; ++i) {
			if(x->get_index_var(i).get_solution_min()) {
				lb.insert(x->elts_ub[i]);
				ub.insert(x->elts_ub[i]);
			} else if(x->children[i].get_solution_max())
				ub.insert(x->elts_ub[i]);
		}
		out << "{";
		for( set<int>::const_iterator i = ub.begin(); i != ub.end(); ++i) {
			if( i != ub.begin() ) out << ", ";
			out << *i;
		}
		out << "}";
	} else if (ai->isBool()) {
		out << (ai->getBool() ? "true" : "false");
	} else if (ai->isSet()) {
		AST::SetLit* s = ai->getSet();
		if (s->interval) {
			out << s->min << ".." << s->max;
		} else {
			out << "{";
			for (unsigned int i=0; i<s->s.size(); i++) {
				out << s->s[i] << (i < s->s.size()-1 ? ", " : "}");
			}
		}
	} else if (ai->isString()) {
		std::string s = ai->getString();
		for (unsigned int i=0; i<s.size(); i++) {
			if (s[i] == '\\' && i<s.size()-1) {
				switch (s[i+1]) {
				case 'n': out << "\n"; break;
				case '\\': out << "\\"; break;
				case 't': out << "\t"; break;
				default: out << "\\" << s[i+1];
				}
				i++;
			} else {
				out << s[i];
			}
		}
	}
}

void
Printer::print(std::ostream& out,
		Solver& solver,
		const IntVarArray& iv,
		const BoolVarArray& bv,
		const SetVarArray& sv) const {
	if (_output == NULL)
		return;
	for (unsigned int i=0; i< _output->a.size(); i++) {
		AST::Node* ai = _output->a[i];
		if (ai->isArray()) {
			AST::Array* aia = ai->getArray();
			int size = aia->a.size();
			out << "[";
			for (int j=0; j<size; j++) {
				printElem(out,solver, aia->a[j],iv,bv,sv);
				if (j<size-1)
					out << ", ";
			}
			out << "]";
		} else {
			printElem(out,solver,ai,iv,bv,sv);
		}
	}
}


Printer::~Printer(void) {
	delete _output;
}

}


FlatZinc::SolutionPrinter::SolutionPrinter(Printer *p, FlatZincModel *fm, Mistral::Solver *s) 
  : p_(p), fm_(fm), solver_(s) {
  //solver_->add((SolutionListener*)this);
}

void FlatZinc::SolutionPrinter::notify_solution() {
  fm_->print_solution(std::cout, *p_);
};

// STATISTICS: flatzinc-any
