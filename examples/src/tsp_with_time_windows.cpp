 

#include "mistral_scheduler.hpp"


using namespace Mistral;

void get_solution(SAT_Model* model, std::vector<int>& start_time) {
    start_time.clear();
    for(auto& t : model->tasks) {
        start_time.push_back(t.get_solution_min());
    }
}

int main( int argc, char** argv )
{

  ParameterList params; //(argc, argv);
  usrand(12345);

  StatisticList stats;
  stats.start();

    
    SAT_Model *model = new SAT_Model(&params, &stats);
    
    
    model->add_task(10, 0, 130);
    model->add_task(5, 50, 150);
    model->add_task(23, 10, 100);
    model->add_task(12, 30, 150);
    model->add_task(30, 0, 110);
    model->add_task(7, 70, 150);
    
    model->add_transition(0,1,3,8);
    model->add_transition(0,2,10,2);
    model->add_transition(0,3,5,5);
    model->add_transition(0,4,4,10);
    model->add_transition(0,5,13,1);
    model->add_transition(1,2,4,4);
    model->add_transition(1,3,20,3);
    model->add_transition(1,4,1,1);
    model->add_transition(1,5,3,5);
    model->add_transition(2,3,10,5);
    model->add_transition(2,4,15,5);
    model->add_transition(2,5,10,10);
    model->add_transition(3,4,4,8);
    model->add_transition(3,5,7,3);
    model->add_transition(4,5,5,20);
    
    
        
      BranchingHeuristic *heu =
          new SchedulingWeightedDegree<TaskDomOverBoolWeight, Guided<MinValue>, 2>(
              model, model->disjunct_map);
    
    
    Outcome result = model->run(heu);
    
    
    std::vector<int> start_time;
    switch(result) {
        case SAT : std::cout << "OK!\n";
            get_solution(model,start_time);
            for(auto i{0}; i<start_time.size(); ++i) {
                std::cout << "t[" << (i+1) << "] = " << start_time[i] << std::endl;
            }
            break;
        case UNSAT : std::cout << "NO SOLUTION!\n";
            break;
        default: std::cout << "Search limit reached\n";
    }

}
  




