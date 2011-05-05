
#include <vector>
#include <iomanip>
#include <fstream>

#include <mistral_search.hpp>
#include <mistral_variable.hpp>
#include <mistral_constraint.hpp>


using namespace std;
using namespace Mistral;


// create X random unit tests of each type, 
// with verbosity V and thoroughness T.
// Then run them


#define LOW     1
#define MEDIUM  2
#define HIGH    3
#define EXTREME 4


class UnitTest {

protected:

  void checkDomainIntegrity(Variable X);

public:

  int Verbosity;
  int Quality;
  int Quantity;
  
  UnitTest();
  UnitTest(const int vb, const int ql, const int qt);
  virtual ~UnitTest();

  virtual void run() = 0;

};


class RandomDomainRandomRemove : public UnitTest {

public:
  
  RandomDomainRandomRemove(const int ql=MEDIUM, const int qt=MEDIUM);
  ~RandomDomainRandomRemove();

  virtual void run();
};


class RandomDomainRandomSetDomainAndRestore : public UnitTest {

public:
  
  RandomDomainRandomSetDomainAndRestore(const int ql=MEDIUM, 
					const int qt=MEDIUM);
  ~RandomDomainRandomSetDomainAndRestore();

  virtual void run();

};


class RandomDomainRandomSetMinAndRestore : public UnitTest {

public:
  
  RandomDomainRandomSetMinAndRestore(const int ql=MEDIUM, 
				     const int qt=MEDIUM);
  ~RandomDomainRandomSetMinAndRestore();

  virtual void run();
};


class RandomDomainRandomSetMaxAndRestore : public UnitTest {

public:
  
  RandomDomainRandomSetMaxAndRestore(const int ql=MEDIUM, 
				     const int qt=MEDIUM);
  ~RandomDomainRandomSetMaxAndRestore();

  virtual void run();
};



class RandomDomainRandomSetDomainBitsetAndRestore : public UnitTest {

public:
  
  RandomDomainRandomSetDomainBitsetAndRestore(const int ql=MEDIUM, 
					      const int qt=MEDIUM);
  ~RandomDomainRandomSetDomainBitsetAndRestore();

  virtual void run();
};



class RandomDomainRandomRemoveRangeAndRestore : public UnitTest {

public:
  
  RandomDomainRandomRemoveRangeAndRestore(const int ql=MEDIUM, 
					  const int qt=HIGH);
  ~RandomDomainRandomRemoveRangeAndRestore();

  virtual void run();
};

template< int NUM_HEAD >
class RandomCListRandomEraseAndRestore : public UnitTest {

public:

  int NUM_ELT;
  
  RandomCListRandomEraseAndRestore(const int ql=HIGH, 
				   const int qt=MEDIUM); 
  ~RandomCListRandomEraseAndRestore();

  class NaiveImplementation : public Reversible {

  public:
    Vector<int> elements;
    Vector<int> whichlist;
    Vector<int> whichleveladded;
    Vector<int> whichlevelremoved;

    NaiveImplementation(Solver *s) ;

    void create(int elt, int wl=0);
    void reversible_erase(int elt, int wl);
    void reversible_add(int elt, int wl=0);
    void save();
    void restore();

    void print(int num_lists);
  };


  void checkEquality(ReversibleMultiList<int,NUM_HEAD>& alist, 
		     NaiveImplementation& blist);
  virtual void run();
};



template< class T >
class RandomRevNumAffectations : public UnitTest {

public:
  
  RandomRevNumAffectations(const int ql=LOW, const int qt=LOW); 
  ~RandomRevNumAffectations();

  virtual void run();  
};


class CostasAllDiffAllSolutions : public UnitTest {

public:
  
  int size;
  int consistency;
  int domain;
  
  CostasAllDiffAllSolutions(const int sz, const int ct, const int dt=DYN_VAR);
  ~CostasAllDiffAllSolutions();

  virtual void run();
};



class CostasNotEqualAllSolutions : public UnitTest {

public:
  
  int size;
  
  CostasNotEqualAllSolutions(const int sz);
  ~CostasNotEqualAllSolutions();

  virtual void run();
};


class Reset : public UnitTest {

public:
  
  Reset();
  ~Reset();

  virtual void run();
};

class Pigeons : public UnitTest {

public:
  
  int size;
  
  Pigeons(const int sz);
  ~Pigeons();

  virtual void run();
};

class BoolPigeons : public UnitTest {

public:
  
  int var_type;
  int size;
  
  BoolPigeons(const int sz, const int vt);
  ~BoolPigeons();

  virtual void run();
};


class VarStackDynamicTest : public UnitTest {

public:
  
  VarStackDynamicTest();
  ~VarStackDynamicTest();

  virtual void run();
};

class ModelTest : public UnitTest {

public:
  
  ModelTest();
  ~ModelTest();

  virtual void run();
};

class RewriteTest1 : public UnitTest {

public:
  
  RewriteTest1();
  ~RewriteTest1();

  virtual void run();
};

class OpshopTest : public UnitTest {

public:
  
  OpshopTest();
  ~OpshopTest();

  virtual void run();
};

class RewriteTest : public UnitTest {

public:
  
  RewriteTest();
  ~RewriteTest();

  virtual void run();
};

int main(int argc, char *argv[])
{  

  usrand(12345);

//   showUint((unsigned int)-2, std::cout);
//   std::cout << std::endl;
//   showUint(3, std::cout);
//   std::cout << std::endl;
//   exit(1);

  std::vector<UnitTest*> tests;

  //tests.push_back(new Reset());
  /*

    tests.push_back(new CostasNotEqual2Solutions(atoi(argv[1]))); 
    tests.push_back(new CostasAllDiff2Solutions(atoi(argv[1])));
    tests.push_back(new CostasNotEqualSolutions(atoi(argv[1]))); 

  */
  //tests.push_back(new CostasAllDiffAllSolutions(8));
  //tests.push_back(new CostasAllDiffAllSolutions(atoi(argv[1])));
  //tests.push_back(new CostasNotEqualAllSolutions(atoi(argv[1]))); 

  //tests.push_back(new VarStackDynamicTest());

  int N = 8; //atoi(argv[1]);
  //

  //tests.push_back(new CostasAllDiffAllSolutions(N, BOUND_CONSISTENCY));
  
  //  tests.push_back(new Pigeons(5)); 
  //tests.push_back(new BoolPigeons(5));

  
  //tests.push_back(new ModelTest());

  tests.push_back(new RewriteTest1());
  tests.push_back(new OpshopTest());
  tests.push_back(new BoolPigeons(N+1, EXPRESSION));
  tests.push_back(new BoolPigeons(N+1, BITSET_VAR));
  tests.push_back(new Pigeons(N+2)); 
  tests.push_back(new CostasAllDiffAllSolutions(N, FORWARD_CHECKING));
  tests.push_back(new CostasAllDiffAllSolutions(N, BOUND_CONSISTENCY, RANGE_VAR));
  tests.push_back(new CostasAllDiffAllSolutions(N, BOUND_CONSISTENCY));
  tests.push_back(new CostasNotEqualAllSolutions(N)); 
  tests.push_back(new RandomCListRandomEraseAndRestore<4>());
  tests.push_back(new RandomDomainRandomRemoveRangeAndRestore());
  tests.push_back(new RandomDomainRandomSetDomainBitsetAndRestore());
  tests.push_back(new RandomDomainRandomSetDomainAndRestore());
  tests.push_back(new RandomDomainRandomSetMaxAndRestore());
  tests.push_back(new RandomDomainRandomSetMinAndRestore());
  tests.push_back(new RandomDomainRandomRemove());
  tests.push_back(new RandomRevNumAffectations<int>());


  //tests[0]->Verbosity = HIGH;
  //tests[0]->Quality = HIGH;
  //tests[0]->Quantity = EXTREME;
  //tests[0]->Verbosity = EXTREME;

  while(tests.size() > 0) {
    UnitTest *t = tests.back();

    //t->Verbosity = EXTREME;
    
    double TIME = getRunTime();
    t->run();
    cout << (getRunTime() - TIME) << endl ;
    
    delete t;
    tests.pop_back();
  }

}


UnitTest::UnitTest() {
  Verbosity=LOW;
  Quality=MEDIUM;
  Quantity=MEDIUM;
}
  
UnitTest::UnitTest(const int vb, const int ql, const int qt) 
{
  Verbosity = vb;
  Quality = ql;
  Quantity = qt;
}

UnitTest::~UnitTest() {
}

void UnitTest::checkDomainIntegrity(Variable X) {
  // get all values, min, max and size through iteration
  int nxt = X.get_min();
  int v=nxt;
  
  int xmin = INFTY;
  int xmax = -INFTY;
  Vector<int> values;
  int max_iteration = 10000;

  do {
    if(--max_iteration < 0) {
      cout << "Error while iterating (infinite loop?)!" << endl
	   << X << " in " << X.get_domain() << " / " << values << endl;
      exit(1);
    }
    
    for(int k=v+1; k<nxt; ++k) {
      if(X.contain(k)) {
	cout << "Error while iterating (some values were missed)!" << endl
	     << X << " in " << X.get_domain() << " / " << values << endl;
	exit(1);
      }
    }

    v=nxt;
    if(xmin > v) xmin = v;
    if(xmax < v) xmax = v;
    values.add(v);

    if(!X.contain(v)) {
      cout << "Error while iterating (got non-member values)!" << endl
	   << X << " in " << X.get_domain() << " / " << values << endl;
      exit(1);
    }

    nxt = X.next(v);    
  } while(nxt != v);

  if(xmin != X.get_min()) {
    cout << "lower bound not properly maintained!" << endl
	 << X << " in " << X.get_domain() << " / " << values 
	 << ": " << xmin << " vs " << X.get_min() << endl;
    exit(1);
  }

  if(xmax != X.get_max()) {
    cout << "upper bound not properly maintained!" << endl
	 << X << " in " << X.get_domain() << " / " << values 
	 << ": " << xmax << " vs " << X.get_max() << endl;
    exit(1);
  }

  if(values.size != X.get_size()) {
    cout << "domain size not properly maintained!" << endl
	 << X << " in " << X.get_domain() << " / " << values 
	 << ": " << values.size << " vs " << X.get_size() << endl;
    exit(1);
  }
}


RandomDomainRandomRemove::RandomDomainRandomRemove(const int ql, 
						   const int qt) 
  : UnitTest(LOW, ql, qt) {}
RandomDomainRandomRemove::~RandomDomainRandomRemove() {}

void RandomDomainRandomRemove::run() {

  //Verbosity = EXTREME;

  /// create vars with random domains and remove a random set of values.
  Vector<int> values;
  Vector<int> rand_values;

  int N = (1<<(1<<Quantity));

  int max_size = (1<<(1<<Quality)); // 2, 4, 16, 256, 65536
  int max_min = (1<<(1<<Quality));
  
  int step = (Quantity ? (1<<Quality)/(Quantity) : INFTY);

  if(Verbosity) cout << "Run random removal checks [2.." 
		     << max_size << "] x [2.." << max_min 
		     << "] (" << step << ") ";// << endl;
  if(Verbosity>LOW) cout << endl;

  for(int dom_size = 2; dom_size <= max_size ; dom_size += step) {
    if(Verbosity > LOW) {
      cout << " do " << (N*(max_min-2)/step) 
	   << " checks with domain size = " << dom_size << endl;
    }
      
    for(int dom_min = 2; dom_min <= max_min ; dom_min += step) {
      
      if(Verbosity > LOW) {
	cout << "  do " << N << " checks with domains in [-" << (dom_min/2) 
	     << ".." << (dom_size-(dom_min/2)) << "] to [" << (dom_min/2) 
	     << ".." << (dom_size+(dom_min/2)) << "]" << endl;
      }

      for(int i=0; i<N; ++i) {

	int lb = randint(dom_min) - (dom_min/2);
	int ub = randint(dom_size) + lb;
	
	int nvalues = dom_size;
	for(int j=0; j<nvalues; ++j) 
	  rand_values.add(randint(2*dom_size)-(ub-lb+1));
	
	Solver s;
	Variable X(lb,ub);

	//std::cout << "X.domain_type: " << X.domain_type << " / " << X.variable << std::endl;

	X.add_to(&s);
	//s.add(X);
	//X = X.get_var();
	
	//std::cout << "X.domain_type: " << X.domain_type << " / " << X.variable << std::endl;

	//s.initialise();
	
	if(Verbosity > MEDIUM) {
	  cout << "    " << X << " in " << X.get_domain() << ": " << endl;
	  cout << "    remove " << rand_values << endl; 
	}

	for(int j=0; j<nvalues; ++j) { 
	  
	  if(Verbosity > HIGH) {
	    cout << "      remove " << rand_values[j] << " from " 
		 << X << " in " << X.get_domain() << " => "; 
	  }
	  Event evt = X.remove(rand_values[j]);
	  if( evt != FAIL_EVENT) {
	    if(Verbosity > HIGH) {
	      cout << "    " << X << " in " << X.get_domain() << endl; 
	    }
	  
	    checkDomainIntegrity(X);
	  
	    if(X.contain(rand_values[j])) {
	      cout << "Error on [remove " << rand_values[j]
		   << " from " << X << " in " << X.get_domain() << endl;
	      exit(1);
	    }
	  } else {
	    if(!X.equal(rand_values[j]))  {
	      cout << "Error on [remove " << rand_values[j]
		   << " from " << X << " in " << X.get_domain() << endl;
	      exit(1);
	    }
	      
	    if(Verbosity > HIGH) {
	      cout << "    Wipe out!" << endl;
	    }
	  }
	}
	
	if(Verbosity > MEDIUM) {
	  cout << " OK (" << X << " in " << X.get_domain() << ")" << endl;
	}

	rand_values.clear();
      }      
    }
  }
}

  
RandomDomainRandomSetDomainAndRestore::RandomDomainRandomSetDomainAndRestore(const int ql, 
									     const int qt) 
  : UnitTest(LOW, ql, qt) {}
RandomDomainRandomSetDomainAndRestore::~RandomDomainRandomSetDomainAndRestore() {}

void RandomDomainRandomSetDomainAndRestore::run() {
  int N = (1<<(1<<Quantity));

  if(Verbosity) cout << "Run " << N << " random set_domain checks ";// << endl;
  if(Verbosity>LOW) cout << endl;

  for(int i=0; i<N; ++i) {
    int dom_size = randint(1<<(1<<Quality))+1;

    if(Verbosity > LOW) cout << "  Domain size: " << dom_size << endl;

    for(int j=0; j<dom_size; ++j) {

      Solver s;
      int dom_min = randint(dom_size) - (dom_size/2);
      int lb = dom_min;
      int ub = dom_min+dom_size;
      Variable X(lb,ub);
      s.add(X);
      //s.initialise();

      if(Verbosity > MEDIUM) cout << "    " << X << " in " << X.get_domain() << endl;

      for(int k=0; k<50; ++k) {

	if(Verbosity > HIGH) cout << "      <- " << lb-k-1 ;

	if(X.set_domain(lb-k-1) != FAIL_EVENT) {
	  cout << "Error on [set_domain] " << lb-k-1
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [FAIL] " << X << " in " << X.get_domain() 
	       << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	  exit(1);
	}   

	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;  

	checkDomainIntegrity(X);
      }

      for(int k=0; k<50; ++k) {

	if(Verbosity > HIGH) cout << "      <- " << lb+k+1 ;

	if(X.set_domain(ub+k+1) != FAIL_EVENT) {
	  cout << "Error on [set_domain] " << ub+k+1
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [FAIL] " << X << " in " << X.get_domain() 
	       << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	  exit(1);
	}

	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;       

	checkDomainIntegrity(X);
      }

      for(int k=0; k<dom_size; ++k) {

	if(Verbosity > HIGH) cout << "      <- " << lb+k ;

	X.set_domain(lb+k);
	if(!X.equal(lb+k)) {
	  cout << "Error on [set_domain] " << lb+k
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}

	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;
	  
	checkDomainIntegrity(X);

	X.restore();
	  
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [restore] " << X << " in " << X.get_domain() << " / [" << lb << ".." << ub << "]" << endl;
	  exit(1);
	}
	  
	checkDomainIntegrity(X);
      }
    }
  }
}


  
RandomDomainRandomSetMinAndRestore::RandomDomainRandomSetMinAndRestore(const int ql, 
								       const int qt) 
  : UnitTest(LOW, ql, qt) {}
RandomDomainRandomSetMinAndRestore::~RandomDomainRandomSetMinAndRestore() {}

void RandomDomainRandomSetMinAndRestore::run() {
  int N = (1<<(1<<Quantity));

  //Verbosity = EXTREME;

  if(Verbosity) cout << "Run " << N << " random set_min checks " ;// << endl;
  if(Verbosity>LOW) cout << endl;

  for(int i=0; i<N; ++i) {
    int dom_size = randint(1<<(1<<Quality))+1;

    if(Verbosity > LOW) cout << "  Domain size: " << dom_size << endl;

    for(int j=0; j<dom_size; ++j) {

      Solver s;
      int dom_min = randint(dom_size) - (dom_size/2);
      int lb = dom_min;
      int ub = dom_min+dom_size;
      Variable X(lb,ub);
      s.add(X);
      //s.initialise();

      if(Verbosity > MEDIUM) cout << "    " << X << " in " << X.get_domain() << endl;

      for(int k=0; k<50; ++k) {

	if(Verbosity > HIGH) cout << "      >= " << lb-k ;

	if(X.set_min(lb-k) != NO_EVENT) {
	  cout << "Error on [set_min] " << lb-k
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [FAIL] " << X << " in " << X.get_domain() 
	       << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	  exit(1);
	}
	  
	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;  
	  
	checkDomainIntegrity(X);
      }
	
      for(int k=0; k<50; ++k) {
	  
	if(Verbosity > HIGH) cout << "      >= " << ub+k+1 ;
	  
	if(X.set_min(ub+k+1) != FAIL_EVENT) {
	  cout << "Error on [set_min] " << ub+k+1
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [FAIL] " << X << " in " << X.get_domain() 
	       << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	  exit(1);
	}

	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;       

	checkDomainIntegrity(X);
      }

      for(int k=0; k<dom_size; ++k) {

	if(Verbosity > HIGH) cout << "      >= " << lb+k+1 ;

	X.set_min(lb+k+1);
	if(X.get_min() < lb+k+1) {
	  cout << "Error on [set_min] " << lb+k+1
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}

	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;
	  
	checkDomainIntegrity(X);

	X.restore();
	  
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [restore] " << X << " in " << X.get_domain() << " / [" << lb << ".." << ub << "]" << endl;
	  exit(1);
	}
	  
	checkDomainIntegrity(X);
      }
    }
  }
}

  
RandomDomainRandomSetMaxAndRestore::RandomDomainRandomSetMaxAndRestore(const int ql, 
								       const int qt) 
  : UnitTest(LOW, ql, qt) {}
RandomDomainRandomSetMaxAndRestore::~RandomDomainRandomSetMaxAndRestore() {}

void RandomDomainRandomSetMaxAndRestore::run() {
  int N = (1<<(1<<Quantity));

  //Verbosity = EXTREME;

  if(Verbosity) cout << "Run " << N << " random set_max checks " ; // << endl;
  if(Verbosity>LOW) cout << endl;

  for(int i=0; i<N; ++i) {
    int dom_size = randint(1<<(1<<Quality))+1;

    if(Verbosity > LOW) cout << "  Domain size: " << dom_size << endl;

    for(int j=0; j<dom_size; ++j) {

      Solver s;
      int dom_min = randint(dom_size) - (dom_size/2);
      int lb = dom_min;
      int ub = dom_min+dom_size;
      Variable X(lb,ub);
      s.add(X);
      //s.initialise();

      if(Verbosity > MEDIUM) cout << "    " << X << " in " << X.get_domain() << endl;

      for(int k=0; k<50; ++k) {

	if(Verbosity > HIGH) cout << "      <= " << lb-k-1 ;

	if(X.set_max(lb-k-1) != FAIL_EVENT) {
	  cout << "Error on [set_max] " << lb-k-1
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [FAIL] " << X << " in " << X.get_domain() 
	       << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	  exit(1);
	}   

	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;

	checkDomainIntegrity(X);  
      }

      for(int k=0; k<50; ++k) {

	if(Verbosity > HIGH) cout << "      <= " << ub+k ;

	if(X.set_max(ub+k) != NO_EVENT) {
	  cout << "Error on [set_max] " << ub+k
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [NO EVENT] " << X << " in " << X.get_domain() 
	       << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	  exit(1);
	}

	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;

	checkDomainIntegrity(X);       
      }

      for(int k=0; k<dom_size; ++k) {

	if(Verbosity > HIGH) cout << "      <= " << lb+k ;

	X.set_max(lb+k);
	if(X.get_max() != lb+k) {
	  cout << "Error on [set_max] " << lb+k
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}

	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;
	  
	checkDomainIntegrity(X);

	X.restore();
	  
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [restore] " << X << " in " << X.get_domain() << " / [" << lb << ".." << ub << "]" << endl;
	  exit(1);
	}
	  
	checkDomainIntegrity(X);
      }
    }
  }
}

  
RandomDomainRandomSetDomainBitsetAndRestore::RandomDomainRandomSetDomainBitsetAndRestore(const int ql, 
											 const int qt) 
  : UnitTest(LOW, ql, qt) {}
RandomDomainRandomSetDomainBitsetAndRestore::~RandomDomainRandomSetDomainBitsetAndRestore() {}

void RandomDomainRandomSetDomainBitsetAndRestore::run() {
  int N = (1<<(1<<Quantity));

  if(Verbosity) cout << "Run " << N << " random set_domain (bitset) checks " ; // << endl;
  if(Verbosity>LOW) cout << endl;

  for(int i=0; i<N; ++i) {
    int dom_size = randint(1<<(1<<Quality))+1;

    if(Verbosity > LOW) cout << "  Domain size: " << dom_size << endl;

    for(int j=0; j<dom_size; ++j) {

      int dom_min = randint(dom_size) - (dom_size/2);
      int lb = dom_min;
      int ub = dom_min+dom_size;

      Solver s;
      Variable X(lb,ub);
      s.add(X);
      //s.initialise();

      if(Verbosity > LOW) cout << "    " << X << " in " << X.get_domain() << endl;

      BitSet vals(lb-50, ub+50, BitSet::empt);

	
      for(int k=0; k<dom_size; ++k) {

	int set_size = randint(dom_size+10);
	  
	for(int l=0; l<set_size; ++l) {
	  vals.add(randint(ub-lb+101) + (lb-50));
	}
	  
	if(Verbosity > MEDIUM) cout << "    Set: " << vals << ": ";

	Event evt = X.set_domain(vals);

	if(Verbosity > MEDIUM) cout << X << " in " << X.get_domain() << endl;

	if(evt != FAIL_EVENT) {
	  int z = 0;
	  int nxt = X.get_min();
	  do {
	    z = nxt;
	    if(!vals.contain(z)) {
	      cout << "Error on [set_domain] " 
		   << " on " << X << " in " << X.get_domain() << endl;
	      exit(1);
	    }
	    nxt = X.next(z);
	  } while( nxt != z );
	    
	  
	  vals.flip();
	  z = 0;
	  nxt = vals.min();
	  do {
	    z = nxt;
	    if(X.contain(z)) {
	      cout << "Error on [set_domain] " 
		   << " on " << X << " in " << X.get_domain() << endl;
	      exit(1);
	    }
	    nxt = vals.next(z);
	  } while( nxt != z );
	    
	  checkDomainIntegrity(X);
	  
	  if(evt != NO_EVENT) X.restore();

	  checkDomainIntegrity(X);

	}
	  
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << X.get_min() << " " << lb << endl;
	  cout << X.get_max() << " " << ub << endl;
	  cout << X.get_size() << " " << (ub-lb+1) << endl;
	  cout << "Error on [restore] " << X << " in " << X.get_domain() << " / [" << lb << ".." << ub << "]" << endl;
	  exit(1);
	}
	  
	checkDomainIntegrity(X);
      }
    }
  }
}

  
RandomDomainRandomRemoveRangeAndRestore::
RandomDomainRandomRemoveRangeAndRestore(const int ql, 
					const int qt) 
  : UnitTest(LOW, ql, qt) {}
RandomDomainRandomRemoveRangeAndRestore::
~RandomDomainRandomRemoveRangeAndRestore() {}

void RandomDomainRandomRemoveRangeAndRestore::run() {
  int N = (1<<(1<<Quantity));

  if(Verbosity) cout << "Run " << N << " random set_domain + removeRange checks " ;// << endl;
  if(Verbosity>LOW) cout << endl;

  //Verbosity = EXTREME;


  for(int i=0; i<N; ++i) {
    int dom_size = randint(1<<(1<<Quality))+1;

    if(Verbosity > LOW) cout << "  Domain size: " << dom_size << endl;

    for(int j=0; j<dom_size; ++j) {

      int dom_min = randint(dom_size) - (dom_size/2);
      int lb = dom_min;
      int ub = dom_min+dom_size;

      Solver s;
      Variable X(lb,ub);
      //s.add(X);
      X.add_to(&s);
      //s.initialise();

      if(Verbosity > MEDIUM) cout << "    " << X << " in " << X.get_domain() << endl;

      BitSet vals(lb-50, ub+50, BitSet::empt);
	
      int set_size = randint(dom_size+10);
	  
      for(int l=0; l<set_size; ++l) {
	vals.add(randint(ub-lb+101) + (lb-50));
      }
	  
      if(Verbosity > MEDIUM) cout << "    Set: " << vals << ": " ;
      s.save();
	
      Event evt1 = X.set_domain(vals);

      if(Verbosity > MEDIUM) cout << X << " in " << X.get_domain() << endl;

      if(evt1 != FAIL_EVENT) {	
	int z = 0;
	int nxt = X.get_min();
	do {
	  z = nxt;
	  if(!vals.contain(z)) {
	    cout << "Error on [set_domain] " 
		 << " on " << X << " in " << X.get_domain() << endl;
	    exit(1);
	  }
	  nxt = X.next(z);
	} while( nxt != z );
	  
	  
	vals.flip();
	z = 0;
	nxt = vals.min();
	do {
	  z = nxt;
	  if(X.contain(z)) {
	    cout << "Error on [set_domain] " 
		 << " on " << X << " in " << X.get_domain() << endl;
	    exit(1);
	  }
	  nxt = vals.next(z);
	} while( nxt != z );
	  
	checkDomainIntegrity(X);
	  
	for(int l=lb; l<ub; ++l) {
	  for(int m=l+1; m<=ub; ++m) {
	    s.save();
	    if(Verbosity > HIGH) cout << "        remove [" << l << ".." 
				      << m << "] from " << X << " in " 
				      << X.get_domain() << ": ";
	    Event evt2 = X.removeRange(l,m);
	    if(Verbosity > HIGH) cout << X << " in " << X.get_domain() << endl;
	    if(evt2 != NO_EVENT && evt2 != FAIL_EVENT) {
	      for(int v=l; v<=m; ++v) {
		if(X.contain(v)) {
		  cout << "Error on [removeRange] " 
		       << " on " << X << " in " << X.get_domain() << endl;
		  exit(1);
		}
	      }
	      checkDomainIntegrity(X);
	    }

	    s.restore();
	  }
	}
      }
      
      s.restore();
	
      if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	cout << X.get_min() << " " << lb << endl;
	cout << X.get_max() << " " << ub << endl;
	cout << X.get_size() << " " << (ub-lb+1) << endl;
	cout << "Error on [restore] " << X << " in " << X.get_domain() << " / [" << lb << ".." << ub << "]" << endl;
	exit(1);
      }
	
      checkDomainIntegrity(X);
    }
  }
}



template< int NUM_HEAD >
RandomCListRandomEraseAndRestore< NUM_HEAD >::
RandomCListRandomEraseAndRestore(const int ql, 
				 const int qt) 
  : UnitTest(LOW, ql, qt) { NUM_ELT = 100*(1 << (1 << Quality)); }

template< int NUM_HEAD >
RandomCListRandomEraseAndRestore< NUM_HEAD >::
~RandomCListRandomEraseAndRestore() {}

template< int NUM_HEAD >
void RandomCListRandomEraseAndRestore< NUM_HEAD >::run() {

  Solver s;
  //s.initialise();

  int N = (1<<(1<<Quantity));
  if(Verbosity) cout << "Run " << N << " random ConstraintList checks " ;
  if(Verbosity>LOW) cout << endl;
    
  for(int iteration=0; iteration<N; ++iteration) {
      
    ReversibleMultiList<int,NUM_HEAD> alist(&s);
    //alist.env = &s;
      
    NaiveImplementation blist(&s);
      
    Vector<unsigned int> wl;
      
    Vector<unsigned int> elts;
      
    Vector<unsigned int> idx1;
    Vector<unsigned int> idx2;
    for(int i=0; i<NUM_ELT; ++i)
      elts.add(i+1);
    for(int i=0; i<NUM_ELT; ++i) {
      int j=randint(NUM_ELT-i)+i;
      elts[i] = elts[j];
      elts[j] = i+1;
    }
      
      
    int i_elts = (1 << (1 << Quality));
    int n_elts = i_elts;
    int head;
    for(int i=0; i<n_elts; ++i) {
      head = randint(NUM_HEAD);
      wl.add(head);

      if( i >= (int)(elts.size) ) {
	cout << i << " " << NUM_ELT << endl;
	exit(1);
      }

      unsigned int the_elt = alist.create( elts[i], head );
      idx1.add( the_elt );
      idx2.add( i );
      blist.create( elts[i], head );
    }
      
    checkEquality(alist, blist);
      
      
    int branch_length = randint(2*i_elts);
      
    BitSet isIn(0, NUM_ELT, BitSet::full);
    isIn.set_max(i_elts-1);
      

    if(Verbosity>LOW) cout << "Do a run on a branch of length " << branch_length << endl;

    for(int i=0; i<branch_length; ++i) {
      s.save();
	
      if(!randint(2))
	for(int k=0; k<2; ++k) {
	  if(isIn.empty()) break;
	    
	  int j = randint(n_elts);
	  while(!isIn.contain(j))
	    j = randint(n_elts);
	    
	  alist.reversible_erase(idx1[j], wl[j]);
	  blist.reversible_erase(idx2[j], wl[j]);
	    
	  if(Verbosity>MEDIUM) {
	    for(int l=0; l<s.level; ++l) cout << " ";
	    cout << s.level << " remove " << (elts[idx2[j]]) << endl;
	    for(int l=0; l<s.level; ++l) cout << " ";
	    //alistcout << alist << endl;
	  }

	  isIn.erase(j);
	  
	}

      if(!randint(3) && n_elts<NUM_ELT-1) {
	// add
	head = randint(NUM_HEAD);
	wl.add(head);
	idx1.add( alist.reversible_add( elts[n_elts], head ) );
	idx2.add( n_elts );
	blist.create( elts[n_elts], head );

	if(Verbosity>MEDIUM) {
	  for(int l=0; l<s.level; ++l) cout << " ";
	  cout << s.level << " add " << (elts[idx2[n_elts]]) << " to " 
	       << wl[n_elts] << "th list" << endl;
	  for(int l=0; l<s.level; ++l) cout << " ";
	  //cout << alist << endl;
	}

	isIn.add(n_elts);

	++n_elts;
      }

      checkEquality(alist, blist);
    }
    
    isIn.fill();
    isIn.set_max(i_elts-1);

    for(int i=0; i<branch_length; i++) {
      s.restore();
      //int j1 = (10+(50-i-1))%n_elts;
      //int j2 = (10+(50-i-2))%n_elts;

      //cout << "add " << (elts[idx2[j1]]) << " and " << (elts[idx2[j2]]) << endl;
      if(Verbosity>MEDIUM) {
	for(int l=0; l<s.level; ++l) cout << " ";
	cout << s.level << " restore" << endl;
	for(int l=0; l<s.level; ++l) cout << " ";
	//cout << alist << endl;      
      }
      
      checkEquality(alist, blist);
    }
  }
}

template< int NUM_HEAD >
void RandomCListRandomEraseAndRestore< NUM_HEAD >::
checkEquality(ReversibleMultiList<int,NUM_HEAD>& alist, 
	      RandomCListRandomEraseAndRestore< NUM_HEAD >::
	      NaiveImplementation& blist) {
  BitSet inA(0, NUM_ELT, BitSet::empt);
  BitSet inB(0, NUM_ELT, BitSet::empt);

  for(int i=0; i<NUM_HEAD; ++i) {
    Node<int> nd = alist.first(i);
    while(alist.next(nd)) {
      inA.add((int)nd);
    }
    for(unsigned int j=0; j<blist.elements.size; ++j) {
      if(blist.whichlist[j] >= i && 
	 blist.whichlevelremoved[j] > blist.env->level &&
	 blist.whichleveladded[j] <= blist.env->level) 
	inB.add(blist.elements[j]);
    }

    if(inA != inB) {
      cout << "Discrepancy between the lists!" << endl
	   << inA << endl
	   << inB << endl;
      exit(1);
    }      
    inA.clear();
    inB.clear();
  }
}

template< int NUM_HEAD >
RandomCListRandomEraseAndRestore< NUM_HEAD >::NaiveImplementation::NaiveImplementation(Solver *s) 
  : Reversible(s) {}

template< int NUM_HEAD >
void RandomCListRandomEraseAndRestore< NUM_HEAD >::NaiveImplementation::create(int elt, int wl) {
  elements.add(elt);
  whichlist.add(wl);
  whichleveladded.add(env->level);
  whichlevelremoved.add(INFTY);
}

template< int NUM_HEAD >
void RandomCListRandomEraseAndRestore< NUM_HEAD >::NaiveImplementation::reversible_erase(int elt, int wl) {
  whichlevelremoved[elt] = env->level;
}

template< int NUM_HEAD >
void RandomCListRandomEraseAndRestore< NUM_HEAD >::NaiveImplementation::reversible_add(int elt, int wl) {
  create(elt, wl);
}

template< int NUM_HEAD >
void RandomCListRandomEraseAndRestore< NUM_HEAD >::NaiveImplementation::save() {}
template< int NUM_HEAD >
void RandomCListRandomEraseAndRestore< NUM_HEAD >::NaiveImplementation::restore() {}

template< int NUM_HEAD >
void RandomCListRandomEraseAndRestore< NUM_HEAD >::NaiveImplementation::print(int num_lists) {
  for(int i=0; i<num_lists; ++i) {
    cout << "[";
    for(unsigned int j=0; j<elements.size; ++j)
      if(whichlist[j] >= i)
	cout << " " << elements[j];
    cout << " ] ";
  }
}



template< class T >
RandomRevNumAffectations< T >::
RandomRevNumAffectations(const int ql, const int qt) 
  : UnitTest(LOW, ql, qt) { }
template< class T >
RandomRevNumAffectations< T >::
~RandomRevNumAffectations() {}

template< class T >
void RandomRevNumAffectations< T >::run() {
  if(Verbosity) cout << "Run " << 10 << " random Reversible<int> checks " ;

  Solver s;
  T x = 1000;
  int val = 0;

  ReversibleNum<T> rx(x, &s);

  int x1 = x;
  int x2;

  //s.initialise();

  for(int k=0; k<5; ++k) {
    s.save();
	
    for(int i=0; i<10; ++i) {
      rx += ((i+k)%20)-10;
      for(int j=0; j<5; ++j)
	val += rx;
    }

    x2 = rx;
      
    for(int i=0; i<999990; ++i) {
      if(!(i%10))s.save();
      rx += ((i+k)%20)-10;
      for(int j=0; j<5; ++j)
	val += rx;
    }
      
    for(int i=0; i<999990; i+=10) {
      s.restore();
    }
      
    if(x2 != rx) {
      cout << "Error: discrepancy between the values! " 
	   << x << " " << rx << endl;
      exit(1);
    }

    s.restore();

    if(x1 != rx) {
      cout << "Error: discrepancy between the values! " 
	   << x << " " << rx << endl;
      exit(1);
    }
  }
}
  
  
CostasAllDiffAllSolutions::CostasAllDiffAllSolutions(const int sz, const int ct, const int dt) 
  : UnitTest(LOW, LOW, LOW) { size=sz; consistency=ct; domain=dt; }
CostasAllDiffAllSolutions::~CostasAllDiffAllSolutions() {}

void CostasAllDiffAllSolutions::run() {

  if(Verbosity) cout << "Run costas array of order 8 (alldiff - "
		     << (consistency==FORWARD_CHECKING ? "FC" : "BC" )
		     << "), look for all solutions "; 

  int i, j;
  VarArray X(size, 1, size, domain);
  Solver s;

  s.add( AllDiff(X) );

  Vector< Variable > scope;
  for(i=1; i<size-1; ++i) {
    scope.clear();
    for(j=0; j<size-i; ++j) {
      scope.add(X[j]-X[j+i]);
    }
    s.add( AllDiff(scope, consistency) );
  }


  s.consolidate();
  //std::cout << s << std::endl;
  //cout << s << endl;

  //s.initialise();
  s.initialise_search(X, 
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());
  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
//     for(i=0; i<size; ++i) {
//       cout << " " << setw(2) << X[i].get_solution_value() ;
//     }
//     cout << endl;
  }

  if(Verbosity) cout << "(" << num_solutions << ") " ;
  if(num_solutions != 444) {
    cout << "Error: wrong number of solutions!" << endl;
    //exit(1);
  }
}

  
CostasNotEqualAllSolutions::CostasNotEqualAllSolutions(const int sz) 
  : UnitTest(LOW, LOW, LOW) { size=sz; }
CostasNotEqualAllSolutions::~CostasNotEqualAllSolutions() {}

void CostasNotEqualAllSolutions::run() {
  if(Verbosity) cout << "Run costas array of order 8 (not equals), look for all solutions "; 
		       

  VarArray X(size, 1, size);

  Solver s;
  int i, j, k;
  for(i=0; i<size; ++i)
    for(j=i+1; j<size; ++j)
      s.add( X[i] != X[j] );

  Vector< Variable > distance[size-2];
  for(i=1; i<size-1; ++i) {
    for(j=0; j<size-i; ++j) {
      distance[i-1].add(X[j] - X[j+i]);
    }
    for(j=1; j<size-i; ++j)
      for(k=0; k<j; ++k) {
	s.add( distance[i-1][j] != distance[i-1][k] );
      }
  }

  //cout << s << endl;

  //s.initialise();
  s.consolidate();
  //std::cout << s << std::endl;

  s.initialise_search(X, 
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
//     for(i=0; i<size; ++i) {
//       cout << " " << setw(2) << X[i].get_solution_value() ;
//     }
//     cout << endl;
    ++num_solutions;

    //    if(num_solutions == 2) exit(1);
  }

  if(Verbosity) cout << "(" << num_solutions << ") " ;
  if(num_solutions != 444) {
    cout << "Error: wrong number of solutions!" << endl;
    //exit(1);
  }

}


Reset::Reset() 
  : UnitTest(LOW, LOW, LOW) { }
Reset::~Reset() {}

void Reset::run() {
    
  VarArray X(3,0,4);

  Solver s;
    
  for(int i=0; i<2; ++i) {
    s.add(X[i] < X[i+1]);
  }

  //s.initialise();

  for(int i=0; i<3; ++i) 
    std::cout << X[i] << " in " << (X[i].get_domain()) << std::endl;

  std::cout << s << std::endl;
    
  s.propagate();

  for(int i=0; i<3; ++i) 
    std::cout << X[i] << " in " << (X[i].get_domain()) << std::endl;
    
  std::cout << s << std::endl;

  s.restore();

  for(int i=0; i<2; ++i) {
    s.add(X[i] < X[i+1]);
  }

  for(int i=0; i<3; ++i) 
    std::cout << X[i] << " in " << (X[i].get_domain()) << std::endl;

    
  std::cout << s << std::endl;

  s.initialise_search(X);

  //     std::cout << s.get_next_solution() << std::endl;

  while(s.get_next_solution() == SAT) {

    //s.depth_first_search(X);

    std::cout << X[0].get_solution_value() << "<"
	      << X[1].get_solution_value() << "<"
	      << X[2].get_solution_value() << std::endl;

  }
    
    

}
  
Pigeons::Pigeons(const int sz) 
  : UnitTest(LOW, LOW, LOW) { size=sz; }
Pigeons::~Pigeons() {}

void Pigeons::run() {
  if(Verbosity) cout << "Run pigeon-hole of order 10: "; 

  VarArray X(size, 1, size-1);
  Solver s;
  int i, j;
  for(i=0; i<size; ++i)
    for(j=i+1; j<size; ++j)
      s.add( X[i] != X[j] );

  s.consolidate();
  //std::cout << s << std::endl;
  //s.initialise();

  s.depth_first_search(X, 
		       new GenericHeuristic< NoOrder, MinValue >(&s), 
		       new NoRestart);
  
  if(s.statistics.num_backtracks != 362879) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    exit(1);
  }

}


BoolPigeons::BoolPigeons(const int sz, const int vt) 
  : UnitTest(LOW, LOW, LOW) { size=sz; var_type=vt; }
BoolPigeons::~BoolPigeons() {}

void BoolPigeons::run() {
  if(Verbosity) cout << "Run pigeon-hole (" 
		     << (var_type == BITSET_VAR ? "bitset" : "boolean")
		     << ") of order " << size << ": "; 

  VarArray X(size*(size-1), 0, 1, var_type);
  VarArray domain;
  Solver s;
  int i, j, k;
  for(i=0; i<size; ++i) {
    domain.clear();
    for(j=i*(size-1); j<(i+1)*(size-1); ++j)
      domain.add( X[j] );
    s.add( BoolSum(domain, size-2) );
    //s.add( BoolSum(X[i*(size-1), (i+1)*(size-1)], size-2) );
    for(j=i+1; j<size; ++j)
      for(k=0; k<size-1; ++k)
	s.add( X[i*(size-1)+k] || X[j*(size-1)+k] );
  }
  //s.initialise();

  //cout << s << endl;
  
  //exit(1);

  s.consolidate();
  //  std::cout << s << std::endl;

  s.depth_first_search(X, 
		       new GenericHeuristic< Lexicographic, MinValue >(&s), 
		       new NoRestart);
  
  if(s.statistics.num_backtracks != 40319 /*362879*/) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }

}


// VariableTest::VariableTest() : UnitTest() {}
// VariableTest::~VariableTest() {}

// void VariableTest::run() {

//   Variable X(0,10);
//   Variable Y(-5,5);

//   Solver s;

//   s.add(X);
//   s.add(Y);

//   cout << s << endl;

// }


ModelTest::ModelTest() : UnitTest() {}
ModelTest::~ModelTest() {}

void ModelTest::run() {

  VarArray X(5,1,7);
  Solver s;

  s.add(X);

  cout << s << endl;
  
  s.add( X[0]<=5 );
  s.add( X[1]<=4 );
  s.add( X[2]<=4 );
  s.add( X[3]<=3 );
  s.add( X[4]>=6 );
  s.add( X[0]!=4 );


  cout << s << endl;

  s.add( (X[0]+2)==X[4] || (X[1]+3)==X[4] );

  s.propagate();
  cout << s << endl;

  s.add( (!(X[0]==4) || !(X[1]==4)) == (X[2] >= X[3]) );

  s.propagate();
  cout << s << endl;

  s.add( (X[2] + X[3]) <= X[4] );

//   s.propagate();
//   cout << s << endl;

//   s.add( (X[2] + X[3]) >= X[4] );

  s.propagate();
  cout << s << endl;

  s.save();
  s.rewrite();
  
  cout << s << endl;

  s.restore();

  cout << s << endl;


  s.initialise_search(X, 
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());
  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {    
    ++num_solutions;
    for(unsigned int i=0; i<5; ++i) {
      cout << " " // << setw(2) << X[i] << ":"
	   << X[i].get_solution_value() ;
    }
    cout << endl;
  }

// 

//   cout << s << endl;


}



RewriteTest::RewriteTest() : UnitTest() {}
RewriteTest::~RewriteTest() {}

void RewriteTest::run() {

  VarArray X(4,1,3);
  Solver s;

  s.add(X);

  Variable Y(1,2);

  s.add(Y);


  cout << "Init\n" << s << endl;

  s.add( (Y==1) <= (X[0]==X[2] || X[1]==X[3]) );
  s.add( (Y==1) <= (X[0]!=X[2] || X[1]!=X[3]) );
  s.add( Y!=2 );

  //s.propagate();

  cout << "Add constraints\n" << s << endl;

  s.consolidate();

  cout << "Consolidate\n" << s << endl;  

  s.save();
  s.rewrite();
  
  cout << "Rewrite\n" << s << endl;

  s.save();
  s.add( X[0]+1 < X[1] );

  cout << "Add a new constraint\n" << s << endl;

  s.consolidate();
  s.rewrite();

  cout << "Rewrite\n" << s << endl;

  s.restore();

  cout << "Restore\n" << s << endl;

  s.restore();

  cout << "Restore\n" << s << endl;



//   s.initialise_search(X, 
// 		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
// 		      new NoRestart());
//   int num_solutions = 0;
//   while(s.get_next_solution() == SAT) {    
//     ++num_solutions;
//     for(unsigned int i=0; i<5; ++i) {
//       cout << " " // << setw(2) << X[i] << ":"
// 	   << X[i].get_solution_value() ;
//     }
//     cout << endl;
//   }

// 

//   cout << s << endl;


}


VarStackDynamicTest::VarStackDynamicTest() : UnitTest() {}

VarStackDynamicTest::~VarStackDynamicTest() {}

void VarStackDynamicTest::run() {

  Queue q;
  q.initialise(7,13);

  q.add(8);

  q.add(7);

  q.add(13);

  q.add(11);

  cout << q << endl;

  cout << q.pop() << endl;

  cout << q << endl;

  cout << q.pop() << endl;

  cout << q << endl;


  q.extend(-3);

  cout << q << endl;
  
  q.add(-3);

  cout << q << endl;


  q.declare(20);

  cout << q << endl;

  cout << q.pop() << endl;

  cout << q << endl;

  q.declare(19);

  cout << q << endl;

  cout << q.pop() << endl;

  cout << q << endl;

  cout << q.pop() << endl;

  cout << q << endl;

  cout << q.pop() << endl;

  cout << q << endl;

 cout << q.pop() << endl;

  cout << q << endl;


  q.add(-1);
  cout << q << endl;

  q.add(5);
  cout << q << endl;

  q.add(15);
  cout << q << endl;


//   IntStack witness;
//   witness.initialise(0,10,true);

//   cout << witness << endl;

//   for(int i=0; i<10; ++i) {
//     int elt = 10+(i+1)*5;
//     cout << "add " << elt << ": ";
//     witness.declare(elt);
//     cout << witness << endl;
//   }

//   for(int i=0; i<10; ++i) {
//     int elt = (i+1)*-5;
//     cout << "add " << elt << ": ";
//     witness.declare(elt);
//     cout << witness << endl;
//   }
  
}


OpshopTest::OpshopTest() : UnitTest() {}
OpshopTest::~OpshopTest() {}


void OpshopTest::run() {

  if(Verbosity) cout << "Run open-shop (GP03-01): "; 

  // read file
  int nJobs, nMachines;
  int i=0, j, k, dur=0, lb, opt, bufsize=1000;
  char buf[bufsize];
  std::ifstream infile( "./examples/data/GP03-01.txt", std::ios_base::in );
	
  do {
    infile.getline( buf, bufsize, '\n' );
  } while( buf[0] == '#' );
	
  while( buf[i] != ' ' ) ++i;
  buf[i] = '\0';
  lb = atoi( buf );
	
  while( buf[i] == '\0' || buf[i] == ' ' || buf[i] == '*' ) ++i;
  j = i;
  while( buf[i] != ' ' && buf[i] != '\n' && buf[i] != '\0' ) ++i;
  buf[i] = '\0';
  opt = atoi( &(buf[j]) );
	
  do {
    infile.get( buf[0] );
    if( buf[0] != '#' ) break;
    infile.getline( buf, bufsize, '\n' );
  } while( true );
  infile.unget();
	
  infile >> nJobs;
  infile >> nMachines;

  int *jobs[nJobs];
  for(i=0; i<nMachines; ++i)
    jobs[i] = new int[nMachines];

  infile.getline( buf, bufsize, '\n' );
	
  do {
    infile.get( buf[0] );
    if( buf[0] != '#' ) break;
    infile.getline( buf, bufsize, '\n' );
  } while( true );
  infile.unget();
	
  k = 0;
  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {
      infile >> jobs[i][j];
      dur += jobs[i][j];
    }
  }

  VarArray task;

  int makespan = (int)((double)opt*1.2);

  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {
      Variable t(0, makespan-jobs[i][j]);
      task.add(t);
    }
  }

  Solver s;

  s.add(task);
  
  for(int i=0; i<nJobs; ++i) 
    {
      for(int j=0; j<nMachines-1; ++j)
	for(int k=j+1; k<nMachines; ++k) {
	  s.add(Disjunctive(task[i*nMachines+j], 
			    task[i*nMachines+k], 
			    jobs[i][j],
			    jobs[i][k]
			    ));
	}
    }

  for(int j=0; j<nMachines; ++j)
    for(int i=0; i<nJobs-1; ++i) 
      {
	for(int k=i+1; k<nJobs; ++k) {
	  s.add(Disjunctive(task[i*nMachines+j], 
			    task[k*nMachines+j], 
			    jobs[i][j],
			    jobs[k][j]
			    ));
	}
    }

  s.consolidate();
  makespan = (int)((double)opt*1.0001);
  for(int i=0; i<nJobs; ++i)
    for(int j=0; j<nMachines; ++j) 
      s.add( task[i*nMachines+j] <= makespan-jobs[i][j] );

  //std::cout << s << std::endl;
  s.specialise();
  //std::cout << s << std::endl;

  s.initialise_search(task, 
		      new GenericHeuristic< GenericDVO < MinDomain >, HalfSplit >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
  }


  if(s.statistics.num_backtracks != 774581) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 771316) {
    cout << "Error: wrong number of solutions! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }

}



RewriteTest1::RewriteTest1() : UnitTest() {}
RewriteTest1::~RewriteTest1() {}

void RewriteTest1::run() {

  if(Verbosity) cout << "Run rewritting test ";

  Variable X(0,10);
  Variable Y(-5,5);
  Variable Z(-100,10);

  int N=10;
  VarArray E(N, -50, 50);

  Vector< Variable > all_vars;
  for(int i=0; i<N; ++i)
    all_vars.add(E[i]);
  all_vars.add(X);
  all_vars.add(Y);
  all_vars.add(Z);

  Solver s;
  s.add(X);
  s.add(Y);
  s.add(Z);

  s.add(X != Z);
  s.add(Y <= Z);

  s.add(X == Y);
  for(int i=0; i<N; ++i) s.add(X == E[i]);

  s.rewrite();

  int sum_degree = X.get_degree() + Y.get_degree() + Z.get_degree();
  for(int i=0; i<N; ++i)
    sum_degree += E[i].get_degree();

  if(sum_degree != 26) {
    cout << "Error while rewritting (wrong number of constraints)!" << endl
	 << "degree(X) = " << X.get_degree() << endl
	 << "degree(Y) = " << Y.get_degree() << endl
	 << "degree(Z) = " << Z.get_degree() << endl
	 << "degree(E[5]) = " << E[5].get_degree() << endl;
    
    exit(1);
  }
  
 s.consolidate();
 s.add(E[4] > 1);
 s.propagate();
 s.consolidate();

 s.initialise_search(all_vars, 
		     new GenericHeuristic< Lexicographic, MinValue >(&s), 
		     new NoRestart());
 
 int num_solutions = 0;
 while(s.get_next_solution() == SAT) {
   Solution sol(all_vars);
   ++num_solutions;
   if(num_solutions > 50) break;
 }
 
 cout << "(" << num_solutions << ") " ;

 if(num_solutions != 26) {
   cout << "Error wrong number of solution!" << endl;
   exit(1);
 }

}
