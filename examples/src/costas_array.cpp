
#include <mistral_variable.hpp>
#include <climits>
#include <iomanip>


using namespace std;
using namespace Mistral;





int main(int argc, char *argv[])
{  
  int i, j, k, N = (argc > 1 ? atoi(argv[1]) : 8);

  Solver s;
  
  Vector< IntVar > X(N);
  for(i=0; i<N; ++i)
    X.add(new VariableWord< unsigned int, 1 >(1,N));

  Vector< IntVar > difference[N-2];
  for(i=0; i<N-2; ++i)
    difference[i].initialise(0,N-i-1);

  Vector< IntVar > scope(3);
  for(i=1; i<N; ++i) 
    for(j=0; j<i; ++j) {
      scope.add(X[i]);
      scope.add(X[j]);
      s.add(new ConstraintNotEqual(scope));
      scope.clear();
    }

  for(i=1; i<N-1; ++i) {
    for(j=i; j<N; ++j) {
      IntVar Y = new VariableWord< unsigned int, 1 >(1-N, N-1);
      difference[i-1].add(Y);
      scope.add(Y);
      scope.add(X[j]);
      scope.add(X[j-i]);
      s.add(new PredicateAdd(scope));
      scope.clear();
    }
    for(j=1; j<N-i; ++j)
      for(k=0; k<j; ++k) {
	scope.add(difference[i-1][j]);
	scope.add(difference[i-1][k]);
	s.add(new ConstraintNotEqual(scope));
	scope.clear();
      }
  }

  s.initialise();
  s.parameters.verbosity = 1;
  s.parameters.time_limit = 10.0;


//   std::cout << "search on all variables but one:" << std::endl;
//   Vector< IntVar > search_vars;
//   for(i=0; i<N-1; ++i) search_vars.add(X[i]);
//   s.depth_first_search(search_vars, NULL, new Geometric(64, 1.333));
//   std::cout << s.statistics << std::endl;

//   for(int i=0; i<N-1; ++i)
//     cout << setw(3) << X[i]->get_value() << " " ;
//   cout << endl;
//   for(i=0; i<N-3; ++i) {
//     for(j=0; j<N-i-2; ++j)
//       cout << setw(3) << difference[i][j]->get_value() << " " ;
//     cout << endl;
//   }
//   //cout << setw(3) << (X[0]->get_value() - X[N-1]->get_value()) << endl << endl;

  std::cout << "\nsearch on everything with min domain" << std::endl;
  s.depth_first_search(X, new GenericDVO<VarSelectorDomain>(&s), new Geometric(256, 1.1));
  std::cout << s.statistics << std::endl;

  for(int i=0; i<N; ++i)
    cout << setw(3) << X[i]->get_value() << " " ;
  cout << endl;
  for(i=0; i<N-2; ++i) {
    for(j=0; j<N-i-1; ++j)
      cout << setw(3) << difference[i][j]->get_value() << " " ;
    cout << endl;
  }
  cout << setw(3) << (X[0]->get_value() - X[N-1]->get_value()) << endl << endl;


}


