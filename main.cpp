#include<iostream>
#include<vector>
#include<ctime>
#include<fstream>

using namespace std;

#define OBJECTIVE 2				/* optimization object number */
#define POPSIZE 100				/* population size */
#define MAXGENS 500				/* max. number of generations */
#define NVARS 30				/* no. of problem variables */
#define PXOVER 0.9				/* probability of crossover */
#define PMUTATION 0.04			/* probability of mutation */

struct gene {
	double gene_val;			/* gene value */
	double upper;				/* gene upper bound */
	double lower;				/* gene lower bound */
};

struct genotype					/* genotype (GT), a member of the population */
{
	vector<gene> genes;			/* a string of variables */
	vector<double> fitness;		/* GT's fitness */
	int rank;					/* pareto rank number */
	double distance;			/* crowding distance */
	vector<genotype *> S;		/* pareto dominance set */
	int n;						/* pareto dominated number */
};

vector<genotype> population;    /* population */
vector<genotype> newpopulation; /* new population */
vector<vector<genotype>> F;		/* non-dominance sorted set */

/* function statement */
void initialize();
void MakeNewPopulation(vector<genotype>& newpop, vector<genotype>& oldpop);
void evaluate(genotype& gt);
bool operator <(const genotype a, const genotype b);
bool operator >(const genotype a, const genotype b);
bool isdominas(genotype p, genotype q);
void FastNonDominatedSort(vector<genotype>& pop);
void SortWithRank(vector<genotype>& I);
void SortWithVal(vector<genotype>& I, int m);
double GetMaxFtns(vector<genotype> I, int m);
double GetMinFtns(vector<genotype> I, int m);
void CrowdingDistanceAssignment(vector<genotype>& I);
inline double random(double a, double b) { return ((double)rand() / RAND_MAX) * (b - a) + a; }
inline int random(int a, int b) { return (rand() % (b - a + 1)) + a; }


/* initialize population */
void initialize() {
	for (int i = 0; i < POPSIZE; ++i) {
		genotype gt_temp;
		for (int j = 0; j < NVARS; ++j) {
			gene g_temp;
			g_temp.lower = 0.0;
			g_temp.upper = 1.0;
			g_temp.gene_val = random(g_temp.lower, g_temp.upper);
			gt_temp.genes.push_back(g_temp);
		}
		evaluate(gt_temp);
		population.push_back(gt_temp);
	}
}

int main() {
	srand((int)time(NULL));
	initialize();

	int gen = 1;	//current generation 

	//fstream fout;
	//fout.open("result.txt", ios::out);

	while (gen < MAXGENS) {

		vector<genotype> Rt;	//put old population and new population into set "Rt"
		Rt.insert(Rt.end(), population.begin(), population.end());
		Rt.insert(Rt.end(), newpopulation.begin(), newpopulation.end());

		FastNonDominatedSort(Rt);	//excute non-dominated sort

		vector<genotype> nextgenpop;
		int _rank = 0;		//non-dominance rank number
		while (_rank < F.size() && nextgenpop.size() + F[_rank].size() <= POPSIZE) { //append F[_rank] into next generation population 
			if (!F[_rank].empty()) 
				CrowdingDistanceAssignment(F[_rank]);	//calculate crowding distance
			for (auto elem : F[_rank]) {
				nextgenpop.push_back(elem);
			}
			++_rank;
		}
		if (_rank == 0) CrowdingDistanceAssignment(F[_rank]);	//if F[0].size > population size
		if (_rank < F.size()) SortWithRank(F[_rank]);
		int rest_num = POPSIZE - nextgenpop.size();
		for (int k = 0; k < rest_num; ++k) {
			nextgenpop.push_back(F[_rank][k]);
		}
		F.clear();

		population = nextgenpop;
		newpopulation.clear();
		MakeNewPopulation(newpopulation, population);	//construct new popuration by conventional GA

		//fout << "-----" << gen << "-----" << endl;
		//for (auto i : population) fout << i.fitness[0] << " " << i.fitness[1] << endl;
		//fout << "--" << endl;
		//for (auto i : newpopulation) fout << i.fitness[0] << " " << i.fitness[1] << endl;

		cout << gen << endl;
		gen++;
	}
	//fout.close();

	/* print result */
	vector<genotype> Rt;
	Rt.insert(Rt.end(), population.begin(), population.end());
	Rt.insert(Rt.end(), newpopulation.begin(), newpopulation.end());
	FastNonDominatedSort(Rt);
	for (auto i : F[0]) cout << i.fitness[0] << " " << i.fitness[1] << endl;
}

/* non-dominated sort, divide population into different rank set:{F} */
void FastNonDominatedSort(vector<genotype>& pop) {
	vector<genotype> Fi;	//first rank set, all the rank 1 individuals are in this set
	/* generate rank 1 and its dominance set:{S} and dominated number "n" */
	for (auto& p : pop) {
		p.S.clear();
		p.n = 0;
		for (auto& q : pop) {
			if (q.fitness != p.fitness) {
				if(p.fitness[0] <= q.fitness[0] && p.fitness[1] <= q.fitness[1]) p.S.push_back(&q);
				else if(q.fitness[0] <= p.fitness[0] && q.fitness[1] <= p.fitness[1]) p.n += 1;
				//if (isdominas(p, q)) p.S.push_back(&q);	
				//else if(isdominas(q, p)) p.n += 1;
			}
		}
		if (p.n == 0) {
			p.rank = 1;
			Fi.push_back(p);
		}
	}
	/* generate others rank throuh dominance set:{S} */
	F.push_back(Fi);
	int i = 0;
	while (!F[i].empty()) {
		vector<genotype> Q;
		for (auto& p : F[i]) {
			for (int q = 0; q < p.S.size(); ++q) {
				p.S[q]->n -= 1;
				if (p.S[q]->n == 0) {
					p.S[q]->rank = (i + 1) + 1;
					Q.push_back(*(p.S[q]));
				}
			}
		}
		i++;
		F.push_back(Q);
	}
}

/* generate new population, include selection, crossover, mutation process */
void MakeNewPopulation(vector<genotype>& newpop, vector<genotype>& oldpop) {
	//selection
	while (newpop.size() < POPSIZE) {
		int ran1 = random(0, POPSIZE - 1);
		int ran2 = random(0, POPSIZE - 1);
		if (oldpop[ran1] < oldpop[ran2]) {
			newpop.push_back(oldpop[ran1]);
		}
		else {
			newpop.push_back(oldpop[ran2]);
		}
	}
	//crossover
	int count = 0;
	genotype* ind[2];
	for (int i = 0; i < POPSIZE; ++i) {
		double randcross = random(0.0, 1.0);
		if (randcross < PXOVER) {
			ind[count] = &newpop[i];
			count++;
			if (count % 2 == 0) {
				int randpoint = random(1, NVARS - 1);
				for (int k = randpoint; k < NVARS; ++k) {
					swap(ind[0]->genes[k], ind[1]->genes[k]);
				}
				count = 0;
			}
		}
	}
	//mutation
	for (int i = 0; i < POPSIZE; ++i) {
		double randmutate = random(0.0, 1.0);
		if(randmutate >= PMUTATION)
			newpop[i].genes[0].gene_val = random(newpop[i].genes[0].lower, newpop[i].genes[0].upper);
		if (randmutate < PMUTATION) {
			int randpos = random(0, NVARS - 1);
			newpop[i].genes[randpos].gene_val = random(newpop[i].genes[randpos].lower, newpop[i].genes[randpos].upper);
		}
		//evaluate
		evaluate(newpop[i]);
	}
}

/* calculate crowding distance in same rank */
void CrowdingDistanceAssignment(vector<genotype>& I) {
	int l = I.size();
	for (auto& i : I) i.distance = 0.0;
	for (int m = 0; m < OBJECTIVE; ++m) {
		SortWithVal(I, m);
		I[0].distance = DBL_MAX, I[l - 1].distance = DBL_MAX;
		double max_ftns = GetMaxFtns(I, m);
		double min_ftns = GetMinFtns(I, m);
		for (int i = 1; i < l - 1; ++i) {
			I[i].distance += (I[i + 1].fitness[m] - I[i - 1].fitness[m]) / (max_ftns - min_ftns);
		}
	}
}

/* evaluate by objective function */
void evaluate(genotype& gt) {
	gt.fitness.clear();

	//ZDT1 func;
	vector<double> f(2);
	f[0] = gt.genes[0].gene_val;
	gt.fitness.push_back(f[0]);

	double gx = 0.0;
	for (int i = 1; i < NVARS; ++i) {
		gx += gt.genes[i].gene_val;
	}
	gx = 1.0 + gx * (9.0 / (NVARS - 1));
	double hx = 1.0 - sqrt(f[0] / gx);
	f[1] = gx * hx;
	gt.fitness.push_back(f[1]);

	//ZDT2 func
	/*vector<double> f(2);
	f[0] = gt.genes[0].gene_val;
	gt.fitness.push_back(f[0]);

	double gx = 0.0;
	for (int i = 1; i < NVARS; ++i) {
		gx += gt.genes[i].gene_val;
	}
	gx = 1.0 + gx * (9.0 / (NVARS - 1));
	double hx = 1.0 - pow(f[0] / gx, 2);
	f[1] = gx * hx;
	gt.fitness.push_back(f[1]);*/
}

/* overload operator */
bool operator <(const genotype a, const genotype b) {	//if a < b, means a.rank < b.rank
	if (a.rank > b.rank) return false;					//or a.rank == b.rank but a.distance > b.distance
	else if (a.rank < b.rank) return true;
	else {
		if (a.distance <= b.distance) return false;
		else return true;
	}
}

bool operator >(const genotype a, const genotype b) {
	if (a.rank < b.rank) return false;
	else if (a.rank > b.rank) return true;
	else {
		if (a.distance > b.distance) return false;
		else return true;
	}
}

/* if p dominate q, return true */
bool isdominas(genotype p, genotype q) {
	for (int i = 0; i < OBJECTIVE; ++i) {
		if (p.fitness[i] > q.fitness[i])
			return false;
	}
	return true;
}

/* quick sort */
int partition(vector<genotype>& vi, int low, int up)
{
	genotype pivot = vi[up];
	int i = low - 1;
	for (int j = low; j < up; j++)
	{
		if (pivot > vi[j])
		{
			i++;
			swap(vi[i], vi[j]);
		}
	}
	swap(vi[i + 1], vi[up]);
	return i + 1;
}
void quickSort(vector<genotype>& vi, int low, int up)
{
	if (low < up)
	{
		int mid = partition(vi, low, up);
		quickSort(vi, low, mid - 1);
		quickSort(vi, mid + 1, up);
	}
}
void SortWithRank(vector<genotype>& I) {
	int size = I.size();
	quickSort(I, 0, size - 1);
}

void SortWithVal(vector<genotype>& I, int m) {
	int size = I.size();
	for (int i = size - 1; i > 0; --i) {
		for (int j = 0; j < i; j++) {
			if (I[j].fitness[m] > I[j + 1].fitness[m]) {
				swap(I[j], I[j + 1]);
			}
		}
	}
}

double GetMaxFtns(vector<genotype> I, int m) {
	double max = I[0].fitness[m];
	for (auto i : I) {
		if (i.fitness[m] > max) max = i.fitness[m];
	}
	return max;
}

double GetMinFtns(vector<genotype> I, int m) {
	double min = I[0].fitness[m];
	for (auto i : I) {
		if (i.fitness[m] < min) min = i.fitness[m];
	}
	return min;
}
