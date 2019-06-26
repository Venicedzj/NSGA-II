#include<iostream>
#include<vector>

using namespace std;

#define OBJECTIVE 2
#define POPSIZE 100               /* population size */
#define MAXGENS 500             /* max. number of generations */
#define NVARS 30                 /* no. of problem variables */
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
	int rank;
	double distance;
	vector<genotype> S;
	int n;
};

vector<genotype> population;    /* population */
vector<genotype> newpopulation; /* new population; */
vector<vector<genotype>> F;

void initialize();
void MakeNewPopulation(vector<genotype>& newpop, vector<genotype>& oldpop);
void evaluate(genotype& gt);
bool operator !=(const genotype a, const genotype b);
bool operator <(const genotype a, const genotype b);
bool operator >(const genotype a, const genotype b);
void FastNonDominatedSort(vector<genotype>& pop);
void SortWithRank(vector<genotype>& I);
void SortWithVal(vector<genotype>& I, int m);
double GetMaxFtns(vector<genotype> I, int m);
double GetMinFtns(vector<genotype> I, int m);
void CrowdingDistanceAssignment(vector<genotype>& I);
inline double random(double a, double b) { return ((double)rand() / RAND_MAX) * (b - a) + a; }
inline int random(int a, int b) { return (rand() % (b - a + 1)) + a; }

int main() {
	initialize();
	//MakeNewPopulation(newpopulation, population);
	int gen = 1;
	while (gen < MAXGENS) {
		vector<genotype> Rt;
		for (auto i : population) Rt.push_back(i);
		for (auto i : newpopulation) Rt.push_back(i);
		FastNonDominatedSort(Rt);
		vector<genotype> nextgenpop;
		int i = 0;
		while (nextgenpop.size() + F[i].size() <= POPSIZE) {
			CrowdingDistanceAssignment(F[i]);
			for (auto elem : F[i]) nextgenpop.push_back(elem);
			++i;
		}
		SortWithRank(F[i]);
		int restnum = POPSIZE - nextgenpop.size();
		for (int k = 0; k < restnum; ++k) nextgenpop.push_back(F[i][k]);
		population = nextgenpop;
		MakeNewPopulation(newpopulation, population);
		gen++;
	}
	vector<genotype> Rt;
	for (auto i : population) Rt.push_back(i);
	for (auto i : newpopulation) Rt.push_back(i);
	FastNonDominatedSort(Rt);
	for (auto i : F[0]) cout << i.fitness[0] << " " << i.fitness[1] << endl;
}

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
void MakeNewPopulation(vector<genotype>& newpop, vector<genotype>& oldpop) {
	//selection
	while (newpop.size() < POPSIZE) {
		int ran1 = random(0, POPSIZE - 1);
		int ran2 = random(0, POPSIZE - 1);
		if (oldpop[ran1] < oldpop[ran2]) newpop.push_back(oldpop[ran1]);
		else newpop.push_back(oldpop[ran2]);
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
					gene temp = ind[0]->genes[k];
					ind[0]->genes[k] = ind[1]->genes[k];
					ind[1]->genes[k] = temp;
				}
			}
		}
	}
	//mutation
	for (int i = 0; i < POPSIZE; ++i) {
		double randmutate = random(0.0, 1.0);
		if (randmutate < PMUTATION) {
			int randpos = random(0, NVARS - 1);
			newpop[i].genes[randpos].gene_val = random(newpop[i].genes[randpos].lower, newpop[i].genes[randpos].upper);
		}
	}
	//evaluate
	for (int i = 0; i < POPSIZE; ++i) {
		evaluate(newpop[i]);
	}
}
void evaluate(genotype& gt) {
	//ZDT func;
	vector<double> f(2);
	f[0] = gt.genes[0].gene_val;
	gt.fitness.push_back(f[0]);

	double gx = 0;
	for (int i = 1; i < NVARS; ++i) {
		gx += gt.genes[i].gene_val;
	}
	gx = 1.0 + gx * (9.0 / (NVARS - 1));
	double hx = 1 - sqrt(f[0] / gx);
	f[1] = gx * hx;
	gt.fitness.push_back(f[1]);
}

bool operator !=(const genotype a, const genotype b) {
	for (int i = 0; i < NVARS; ++i) {
		if(a.genes[i].gene_val != b.genes[i].gene_val) return true;
	}
	return false;
}

bool operator <(const genotype a, const genotype b) {
	if (a.rank > b.rank) return false;
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

bool isdominas(genotype p, genotype q) {
	for (int i = 0; i < OBJECTIVE; ++i) {
		if (p.fitness[i] > q.fitness[i])
			return false;
	}
	return true;
}

void FastNonDominatedSort(vector<genotype>& pop) {
	vector<genotype> Fi;
	int count = 0;
	for (auto &p : pop) {
		count += 1;
		cout << count << endl;
		p.n = 0;
		for (auto &q : pop) {
			if (q != p) {
				if (isdominas(p, q)) 
					p.S.push_back(q); //目前尚未初始化rank，<需重新定义定义
				else if (isdominas(q, p)) 
					p.n += 1;
			}
		}
		if (p.n == 0) {
			p.rank = 1;
			Fi.push_back(p);
		}
	}
	F.push_back(Fi);
	int i = 0;
	while (!F[i].empty()) {
		vector<genotype> Q;
		for (auto &p : F[i]) {
			for (auto &q : p.S) {
				q.n -= 1;
				if (q.n == 0) {
					q.rank = (i + 1) + 1;
					Q.push_back(q);
				}
			}
		}
		i += 1;
		//if(!Q.empty()){
			F.push_back(Q);
		//}
	}
}
void SortWithRank(vector<genotype>& I) {
	int size = I.size();
	for (int i = size - 1; i > 0; --i) {
		for (int j = 0; j < i; j++) {
			if (I[j] > I[j + 1]) {		//>
				genotype temp = I[j];
				I[j] = I[j + 1];
				I[j + 1] = temp;
			}
		}
	}
}

void SortWithVal(vector<genotype>& I, int m) {
	int size = I.size();
	for (int i = size - 1; i > 0; --i) {
		for (int j = 0; j < i; j++) {
			if (I[j].fitness[m] > I[j + 1].fitness[m]) {
				genotype temp = I[j];
				I[j] = I[j + 1];
				I[j + 1] = temp;
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

void CrowdingDistanceAssignment(vector<genotype> &I) {
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