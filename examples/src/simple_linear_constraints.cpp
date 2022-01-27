#include <mistral_solver.hpp>
#include <mistral_variable.hpp>
#include <mistral_search.hpp>


using namespace std;
using namespace Mistral;

class FilteringStatisticalParity {

  private:

    VarArray scope;
  Solver s;
  Outcome result;
  //
  int nb_protected, nb_unprotected, nb_protected_negative, nb_unprotected_negative;

  //L is a lower bound for the accuracy
  int L;

  //U is a lower bound for the accuracy
  int U;
  //fairness_tolerence is a float in [0,1] that represents the unfairness tolerated.
  //If this value is 0 then the model is 100% fair.
  float fairness_tolerence;

  public:
    /*	ub_sp_plus is the number of examples protected with a positive class
    	ub_sp_minus is the number of examples protected with a negative class
    	ub_su_plus is the number of examples unprotected with a positive class
    	ub_su_minus is the number of examples unprotected with a negative class
     */
    FilteringStatisticalParity(
      int nb_sp_plus,
      int nb_sp_minus,
      int nb_su_plus,
      int nb_su_minus,
      int L,
      int U,
      float fairness_tolerence):
    nb_protected(nb_sp_plus + nb_sp_minus),
    nb_unprotected(nb_su_plus + nb_su_minus),
    nb_protected_negative(nb_sp_minus),
    nb_unprotected_negative(nb_su_minus),
    L(L),
    U(U),
    fairness_tolerence(fairness_tolerence) {

      //scope[0] is sp_plus
      scope.add(Variable(0, nb_sp_plus));
      //scope[1] is sp_minus
      scope.add(Variable(0, nb_sp_minus));
      //scope[2] is su_plus
      scope.add(Variable(0, nb_su_plus));
      //scope[3] is su_minus
      scope.add(Variable(0, nb_su_minus));

      result = UNSAT;

      int constant;

      std::vector < int > accuracy_coefficient;
      // accuracy_coefficient represents the coefficients used in the accuracy.
      accuracy_coefficient.push_back(1);
      accuracy_coefficient.push_back(-1);
      accuracy_coefficient.push_back(1);
      accuracy_coefficient.push_back(-1);

      //Accuracy constraints
      constant = U - nb_sp_minus - nb_su_minus;
      //std::cout <<  " \n c in the CP model, the upper Bound for accuracy constraint: "  << constant	<< std::endl;
      s.add(Sum(scope, accuracy_coefficient) <= constant);

      constant = L - nb_sp_minus - nb_su_minus;

      //std::cout <<  " c in the CP model, the lower Bound for accuracy constraint: "  << constant	<< std::endl;
      s.add(Sum(scope, accuracy_coefficient) >= constant);

      //Fairness constraints
      std::vector < int > fairness_coefficient;

      fairness_coefficient.push_back(nb_unprotected);
      fairness_coefficient.push_back(nb_unprotected);

      fairness_coefficient.push_back(-nb_protected);
      fairness_coefficient.push_back(-nb_protected);

      int fairness_bound = (int)((float) nb_protected * (float) nb_unprotected * fairness_tolerence);

      std::cout << " c fairness_bound : " << fairness_bound << std::endl;
      //Variable bound_fairness (-fairness_bound, fairness_bound);

      s.add(Sum(scope, fairness_coefficient) <= fairness_bound);
      s.add(Sum(scope, fairness_coefficient) >= (-fairness_bound));

      //An alternative way:
      //fairness_coefficient.push_back(-1);
      //scope.add(bound_fairness);
      //s.add( Sum(scope, fairness_coefficient, 0, 0) );
      //s.add( Sum(scope, fairness_coefficient) == bound_fairness) ;
    }

  void run(int verbosity) {
    s.parameters.verbosity = verbosity;
    //std::cout <<  s << std::endl;
    //s.rewrite();
    s.consolidate();
    result = s.solve();
  }

  void print_statistics() {
    s.statistics.print_full(std::cout);
  }

  void print_and_verify_solution() {

    if (result) {
      std::cout << " \n \n c Solution Found! " << std::endl;

      //for ( unsigned int i= 0 ; i< scope.size ; ++i)
      //	std::cout <<  " c Solution value of scope[" << i <<  "]  is " << scope[i].get_solution_int_value() << std::endl;

      int sp_plus = scope[0].get_solution_int_value();
      int sp_minus = scope[1].get_solution_int_value();
      int su_plus = scope[2].get_solution_int_value();
      int su_minus = scope[3].get_solution_int_value();

      //int fairness_sum= scope[4].get_solution_int_value();

      std::cout << " c nb protected positive that are predicted positive: " << sp_plus << std::endl;
      std::cout << " c nb protected negative that are predicted positive: " << sp_minus << std::endl;
      std::cout << " c nb unprotected positive that are predicted positive: " << su_plus << std::endl;
      std::cout << " c nb unprotected negative that are predicted positive:  " << su_minus << std::endl;
      //std::cout <<  " c nb fairness_sum  " <<  fairness_sum << std::endl;

      int accuracy = sp_plus + su_plus - sp_minus - sp_minus + nb_protected_negative + nb_unprotected_negative;

      float fairness = ((float)(sp_plus + sp_minus) / (float) nb_protected) - ((float)(su_plus + su_minus) / (float) nb_unprotected);

      std::cout << " c Accuracy ( Number of examamples well classified) is " << accuracy << " out of " << nb_protected + nb_unprotected << std::endl;

      std::cout << " c Fairness (as float) is " << fairness << std::endl;

      assert(accuracy >= L);
      assert(accuracy <= U);
      assert(fairness <= fairness_tolerence);
      assert(fairness >= -fairness_tolerence);
      //int discrete_fairness= (nb_unprotected * (sp_plus+sp_minus)) - (nb_protected * (su_plus +su_minus) )  ;
      //std::cout <<  " c Fairness (discrete constraint) is "  << discrete_fairness << std::endl;
      std::cout << " c Solution Verified" << std::endl;

    } else
      std::cout << " c No Solution! " << std::endl;
  }

};

int main(int argc, char * argv[]) {

  int nb_sp_plus = 25;
  int nb_sp_minus = 35;
  int nb_su_plus = 20;
  int nb_su_minus = 40;
  int L = 61;
  int U = 105;
  float tolerence = 0.1;

  int total = nb_sp_plus + nb_sp_minus + nb_su_plus + nb_su_minus;
  assert(total > U);
  assert(total < (L * 2));

  std::cout << " c Instance data :  " << std::endl;
  std::cout << " c total number of examples " << total << std::endl;
  std::cout << " c total number of protected examples " << nb_sp_plus + nb_sp_minus << std::endl;
  std::cout << " c total number of unprotected examples " << nb_su_plus + nb_su_minus << std::endl;

  std::cout << " c total number of protected examples with positive class: " << nb_sp_plus << std::endl;
  std::cout << " c total number of protected examples with negative class: " << nb_sp_minus << std::endl;
  std::cout << " c total number of unprotected examples with positive class: " << nb_su_plus << std::endl;
  std::cout << " c total number of unprotected examples with negative class: " << nb_su_minus << std::endl;
  std::cout << " c Lower bound for accuracy: " << L << std::endl;
  std::cout << " c Upper bound for accuracy: " << U << std::endl;
  std::cout << " c Unfairness tolerance (should be in [0,1]): " << tolerence << std::endl;

  FilteringStatisticalParity check_bounds(nb_sp_plus, nb_sp_minus, nb_su_plus, nb_su_minus, L, U, tolerence);

  std::cout << " c run the solver .. and wait for magic  " << std::endl;
  check_bounds.run(0);
  check_bounds.print_statistics();
  check_bounds.print_and_verify_solution();
}
