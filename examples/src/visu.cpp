

#include "mistral_scheduler.hpp"

#define VERBOSE

using namespace Mistral;



int main( int argc, char** argv )
{

  Solver solver;

  solver.parameters.verbosity = 3;



  std::vector<int> p;
  VarArray start_time;
  VarArray end_time;




  // start_time.add(Variable(0,37));
  // end_time.add(Variable(13,50));
  // solver.add(start_time.back() + 13 <= end_time.back());
  // p.push_back(13);


  // start_time.add(Variable(30,35));
  // end_time.add(Variable(35,40));
  // solver.add(start_time.back() + 5 <= end_time.back());
  // p.push_back(5);

  // start_time.add(Variable(33,40));
  // end_time.add(Variable(38,45));
  // solver.add(start_time.back() + 5 <= end_time.back());
  // p.push_back(5);


  // start_time.add(Variable(35,42));
  // end_time.add(Variable(40,47));
  // solver.add(start_time.back() + 5 <= end_time.back());
  // p.push_back(5);



  start_time.add(Variable(0,37));
  end_time.add(Variable(13,50));
  solver.add(start_time.back() + 13 <= end_time.back());
  p.push_back(13);


  start_time.add(Variable(10,15));
  end_time.add(Variable(15,20));
  solver.add(start_time.back() + 5 <= end_time.back());
  p.push_back(5);

  start_time.add(Variable(13,20));
  end_time.add(Variable(18,25));
  solver.add(start_time.back() + 5 <= end_time.back());
  p.push_back(5);


  start_time.add(Variable(15,22));
  end_time.add(Variable(20,35));
  solver.add(start_time.back() + 5 <= end_time.back());
  p.push_back(5);


  start_time.add(Variable(27,31));
  end_time.add(Variable(31,34));
  solver.add(start_time.back() + 3 <= end_time.back());
  p.push_back(3);

  start_time.add(Variable(28,32));
  end_time.add(Variable(32,35));
  solver.add(start_time.back() + 3 <= end_time.back());
  p.push_back(3);

  start_time.add(Variable(29,33));
  end_time.add(Variable(33,36));
  solver.add(start_time.back() + 3 <= end_time.back());
  p.push_back(3);





  // start_time.add(Variable(0,9));
  // end_time.add(Variable(4,10));
  // solver.add(start_time.back() + 4 <= end_time.back());
  // p.push_back(4);


  // start_time.add(Variable(3,3));
  // end_time.add(Variable(5,5));
  // solver.add(start_time.back() + 2 <= end_time.back());
  // p.push_back(2);




// bool multi{true};
// if(multi) {
//   start_time.add(Variable(0,9));
//   end_time.add(Variable(2,11));
//   solver.add(start_time.back() + 2 <= end_time.back());
//   p.push_back(2);


//   start_time.add(Variable(2,2));
//   end_time.add(Variable(5,5));
//   solver.add(start_time.back() + 3 <= end_time.back());
//   p.push_back(3);


//   start_time.add(Variable(4,9));
//   end_time.add(Variable(6,11));
//   solver.add(start_time.back() + 2 <= end_time.back());
//   p.push_back(2);

//   start_time.add(Variable(6,9));
//   end_time.add(Variable(8,11));
//   solver.add(start_time.back() + 2 <= end_time.back());
//   p.push_back(2);

//   start_time.add(Variable(6,9));
//   end_time.add(Variable(8,11));
//   solver.add(start_time.back() + 2 <= end_time.back());
//   p.push_back(2);

//   start_time.add(Variable(0,16));
//   end_time.add(Variable(4,20));
//   solver.add(start_time.back() + 4 <= end_time.back());
//   p.push_back(4);

// } else {

//   start_time.add(Variable(0,10));
//   end_time.add(Variable(1,11));
//   solver.add(start_time.back() + 1 <= end_time.back());
//   p.push_back(1);

//   start_time.add(Variable(0,10));
//   end_time.add(Variable(1,11));
//   solver.add(start_time.back() + 1 <= end_time.back());
//   p.push_back(1);


//   start_time.add(Variable(2,4));
//   end_time.add(Variable(3,5));
//   solver.add(start_time.back() + 1 <= end_time.back());
//   p.push_back(1);

//     start_time.add(Variable(2,4));
//   end_time.add(Variable(3,5));
//   solver.add(start_time.back() + 1 <= end_time.back());
//   p.push_back(1);

//     start_time.add(Variable(2,4));
//   end_time.add(Variable(3,5));
//   solver.add(start_time.back() + 1 <= end_time.back());
//   p.push_back(1);


//   start_time.add(Variable(4,10));
//   end_time.add(Variable(5,11));
//   solver.add(start_time.back() + 1 <= end_time.back());
//   p.push_back(1);

//   start_time.add(Variable(4,10));
//   end_time.add(Variable(5,11));
//   solver.add(start_time.back() + 1 <= end_time.back());
//   p.push_back(1);

//   start_time.add(Variable(6,10));
//   end_time.add(Variable(7,11));
//   solver.add(start_time.back() + 1 <= end_time.back());
//   p.push_back(1);

//   start_time.add(Variable(6,10));
//   end_time.add(Variable(7,11));
//   solver.add(start_time.back() + 1 <= end_time.back());
//   p.push_back(1);

//   start_time.add(Variable(6,10));
//   end_time.add(Variable(7,11));
//   solver.add(start_time.back() + 1 <= end_time.back());
//   p.push_back(1);

//   start_time.add(Variable(6,10));
//   end_time.add(Variable(7,11));
//   solver.add(start_time.back() + 1 <= end_time.back());
//   p.push_back(1);

//   start_time.add(Variable(0,19));
//   end_time.add(Variable(1,20));
//   solver.add(start_time.back() + 1 <= end_time.back());
//   p.push_back(1);

//   start_time.add(Variable(0,19));
//   end_time.add(Variable(1,20));
//   solver.add(start_time.back() + 1 <= end_time.back());
//   p.push_back(1);

//   start_time.add(Variable(0,19));
//   end_time.add(Variable(1,20));
//   solver.add(start_time.back() + 1 <= end_time.back());
//   p.push_back(1);

//   start_time.add(Variable(0,19));
//   end_time.add(Variable(1,20));
//   solver.add(start_time.back() + 1 <= end_time.back());
//   p.push_back(1);

// }


  solver.add(PreemptiveNoOverlap(start_time, end_time, p));




  solver.consolidate();

  // std::cout << solver << std::endl;

  BranchingHeuristic *heuristic;
  RestartPolicy *restart;

  restart = new Luby();
  restart->base = 128;
  heuristic =
      new LastConflict<GenericDVO<MinDomainOverWeight, 1, ConflictCountManager>,
                       SolutionGuided<MinValue, MinValue>,
                       SolutionGuided<MinValue, MinValue>, 1>(&solver);




  solver.initialise_search(solver.variables, heuristic, restart);


  std::cout << "avant:\n";
  for(auto i{0}; i<start_time.size; ++i) {
    std::cout << start_time[i].get_domain() << " - "  << p[i] << " - " << end_time[i].get_domain() << std::endl;
  }
  std::cout << std::endl;


  solver.propagate();


  std::cout << "apres:\n";
  for(auto i{0}; i<start_time.size; ++i) {
    std::cout << start_time[i].get_domain() << " - "  << p[i] << " - " << end_time[i].get_domain() << std::endl;
  }
  std::cout << std::endl;

 
}
