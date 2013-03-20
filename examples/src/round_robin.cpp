//A particular case of round robin tournaments.
//Details came later..

#include <mistral_solver.hpp>
#include <mistral_variable.hpp>
#include <mistral_search.hpp>
#include <tclap/CmdLine.h>

using namespace std;
using namespace Mistral;
using namespace TCLAP;



class RoundRobinModel : public Solver
{


public :


	RoundRobinModel()
{
		nb_players= 6;
		nb_periods= 5;
}
	//the matrix variables with play[i][j]=1 iff. player i plays at period j
	Vector< VarArray > play;

	void setup()
	{
		play.resize(nb_periods);
		for(int i=0; i<nb_periods; ++i)
		{
			play[i].initialise(nb_players, 0, 1);
			add(play[i]);
		}

		//initialisation
		//FREE?????
		//sum = 3
		VarArray subsequence;
		for(int i=0; i<nb_periods; ++i)
		{
			subsequence.clear();
			for(int j=0; j<nb_players; ++j) subsequence.add(play[i][j]);
			add( BoolSum( subsequence, 3, 3) );
		}


		//third constraint
		for(int i=1; i<nb_periods; ++i)
		{
			subsequence.clear();
			for(int j=0; j<nb_players; ++j) subsequence.add(play[i][j]!=play[i-1][j]);
			add( BoolSum( subsequence, 1, INFTY) );
		}

		//forth constraint
		int lower = 3*nb_periods/nb_players;
		int upper = lower+1;
		std::cout << "lower" << lower << " upper" << upper << std::endl;

		for(int j=0; j<nb_players; ++j)
		{
			subsequence.clear();
			for(int i=0; i<nb_periods; ++i) subsequence.add(play[i][j]);
			add( BoolSum( subsequence, lower, upper) );
		}

		//fifth constraint "atleast constraints"
		int n =nb_periods-2;
		for(int j=0; j<nb_players; ++j)
			for(int i=0; i<n; ++i)
			{
				subsequence.clear();
				for(int k=0; k<3; ++k)
					subsequence.add(play[i+k][j]);
				add( BoolSum( subsequence, 1, INFTY) );
			}


		//sixth constraint
		//z[a][b][c] <--> a:period, b and c : players
		int lower_eq = 3*nb_periods/15;
		int upper_eq = lower_eq+1;
		lower_eq--;
		std::cout<< "lower_eq" << lower_eq << " upper_eq" << upper_eq << '\n';

		Vector <Vector <VarArray> > z;

		VarArray __match;
		Vector <VarArray>  __matches;
		for(int period=0; period<nb_periods; ++period)
		{
			__matches.clear();
			//__match.clear();

			for(int b=0; b<nb_players; ++b)
			{
				__match.clear();

				for(int c=0; c<nb_players ; ++c)
				{
					if(c==b)
						__match.add(0);
					else{
						subsequence.clear();
						subsequence.add(play[period][b]);
						subsequence.add(play[period][c]);
						//z.add(BoolSum( subsequence, 2, 2));
						__match.add(BoolSum(subsequence)==2);
						//					z.add(play[period][b]&&play[period][c]);
						//					z.add((play[period][b]+ play[period][c])=2);
					}
				}
				__matches.add(__match);
			}
			z.add(__matches);
		}
		//add(z);
		for(int b=0; b<nb_players; ++b)
			for(int c=0; c<nb_players; ++c)
			{
				subsequence.clear();
				for(int period=0; period<nb_periods; ++period)	subsequence.add(z[period][b][c]);
				add(BoolSum(subsequence, lower_eq, upper_eq));
			}

	}
	void __get_var()
	{

		for(int period=0; period<nb_periods; ++period)
		{
			for(int player=0; player<nb_periods; ++player)
			{
				play[period][player] = play[period][player].get_var();
			}
		}
	}

	int getNbPeriods() const {
		return nb_periods;
	}

	void setNbPeriods(int nbPeriods) {
		nb_periods = nbPeriods;
	}

	int getNbPlayers() const {
		return nb_players;
	}

	void setNbPlayers(int nbPlayers) {
		nb_players = nbPlayers;
	}

private:

	int nb_players, nb_periods;

};


int main(int argc, char **argv)
{
	try
	{

		//	SolverCmdLine cmd("Mistral (round_robin)", ' ', "2.0");
		//	cmd.parse( argc, argv );


		RoundRobinModel * solver = new RoundRobinModel();
		solver->setup();
		solver->__get_var();
		Outcome result = solver->depth_first_search();

		//			Outcome result = solver->depth_first_search(solver->class_at_position, heuristic, restartp);

		if(result)
		{
			//if(cmd.print_solution())
			//{
			for(int i=0; i<solver->getNbPeriods(); ++i)
			{
				Solution sol(solver->play[i]);
				cout << " c  solution: " << sol << endl;
			}
			//}
		}
		//		if(cmd.print_statistics()) {
		cout << solver->statistics << endl;
		//	}

		delete solver;
	}
	catch (ArgException &e)  // catch any exceptions
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
	}

}


