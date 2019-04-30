#include<iostream>
#include<vector>

using namespace std;

#define OBJECTIVE 2
#define POPSIZE 100               /* population size */
#define MAXGENS 500             /* max. number of generations */
#define NVARS 2                 /* no. of problem variables */
#define PXOVER 0.7               /* probability of crossover */
#define PMUTATION 0.07           /* probability of mutation */

struct gene {
	double gene_val;
	double upper;
	double lower;
};

struct genotype /* genotype (GT), a member of the population */
{
	vector<gene> genes;        /* a string of variables */
	vector<double> fitness;            /* GT's fitness */
	vector<double> rfitness;           /* relative fitness */
	vector<double> cfitness;           /* cumulative fitness */
};

struct genotype population[POPSIZE + 1];    /* population */
struct genotype newpopulation[POPSIZE + 1]; /* new population; */

int main() {
	initialize();
	evaluate();
	keep_the_best();
}