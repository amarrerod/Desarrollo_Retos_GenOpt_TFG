#include <cstdlib>
#include <iostream>
#include <signal.h>
#include <math.h>
#include "OBL_CPSO.h"

//#define DEBUG

// Algorithm initialization
bool OBL_CPSO::init(const vector<string>& params) {
	// Check number of parameters
	bool defaultArgs = false;
	if (params.size() == 1) {
		cout << "seed" << endl;
		defaultArgs = true;
	} else if (params.size() > 1 && params.size() != PARAMS) {
		cout << "Parametros: seed popSize INIT" << endl;
		cout << "INIT" << endl;
		cout << "\t RAND" << endl;
		cout << "\t OBL" << endl;
		cout << "\t QOBL" << endl;
		cout << "\t QROBL" << endl;
		return false;
	} else if (params.size() < 1) {
		cerr << "Parametros: seed" << endl;
		return false;
	}
// Only mono-objective optimization is supported
	if (getSampleInd()->getNumberOfObj() != 1) {
		cout << "Multi-Objective not supported" << endl;
		return false;
	}
	if (defaultArgs) {
		setPopulationSize(DEFAULT_POP);
		setInertia(DEFAULT_W);
		setTypeOfInit(RAND);
	} else {
		setPopulationSize(atoi(params[1].c_str()));
		string init = params[2];
		if (init == "RAND")
			setTypeOfInit(RAND);
		else if (init == "OBL")
			setTypeOfInit(OBL);
		else if (init == "QOBL")
			setTypeOfInit(QOBL);
		else if (init == "QROBL")
			setTypeOfInit(QROBL);
		else {
			cerr << "Not init defined. RAND selected by default" << endl;
			setTypeOfInit(RAND);
		}
	}
	seed = static_cast<uint64_t>(atof(params[0].c_str()));
	srand(seed);
	initRandomSeed(seed);
	mean = getSampleInd()->internalClone();
	best = getSampleInd()->internalClone();
	best->setObj(0, numeric_limits<double>::max()); // iniciamos al maximo
	randomPop.resize(getPopulationSize(), 0);
#ifdef DEBUG
	cout << "Seed: " << seed << endl;
	cout << "Population size: " << getPopulationSize() << endl;
	cout << "INIT: " << getTypeOfInit();
	cout << " (RAND: " << RAND << " OBL: " << OBL
	     << " QOBL: " << QOBL << " QROBL: " << QROBL << ")" << endl;
#endif
	return true;
}



void OBL_CPSO::runGeneration() {
	sortPopulation();
	if (getGeneration() == 0) {
		initVelocity();
		best = (*population)[0]->internalClone();  // Al inicio esta claro
		createRandomPop();
		createMean();
		compete();
	}
	if (getGeneration() >= 1) { // Evitamos re-evaluar en la primera
		evaluateSwarm();
		createRandomPop();
		createMean();
		compete(); // Competir y evaluar
	}
	globalSearch();
}

void OBL_CPSO::initVelocity() {
	const int dimension = (*population)[0]->getNumberOfVar();
	velocity.resize(getPopulationSize(), vector<double>(dimension));
	for (int i = 0; i < getPopulationSize(); i++) {
		for (int j = 0; j < dimension; j++) {
			velocity[i][j] = 0.0;
		}
	}
}

void OBL_CPSO::evaluateSwarm() {
	for (Individual* ind : *population) {
		evaluate(ind);
		if (ind->getInternalOptDirection(0) == MINIMIZE
		    && ind->getObj(0) < best->getObj(0))
			best = ind->internalClone();
		else if (ind->getInternalOptDirection(0) == MAXIMIZE
		         && ind->getObj(0) > best->getObj(0))
			best = ind->internalClone();
	}
}

void OBL_CPSO::createRandomPop() {
	vector<int> temp(getPopulationSize());
	int n = {0};
	generate(temp.begin(), temp.end(), [&n] { return n++; });
	int last = getPopulationSize() - 1;
	for (int i = 0; i < getPopulationSize(); i++) {
		int k = getRandomInteger0_N(last);
		randomPop[i] = temp[i];
		temp[k] = temp[last];
		last--;
	}
}

void OBL_CPSO::createMean() {
	const int dim = (*population)[0]->getNumberOfVar();
	const int pop = getPopulationSize();
	for (int j = 0; j < dim; j++) {
		mean->setVar(j, 0.0);
		for (int i = 0; i < pop; i++) {
			mean->setVar(j, mean->getVar(j) + (*population)[i]->getVar(j));
		}
		mean->setVar(j, mean->getVar(j) / pop);
	}
}

void OBL_CPSO::compete() {
	const int pop = getPopulationSize();
	for (int i = 0; i < pop / 3; i++) {
		int r1 = randomPop[i];
		int r2 = randomPop[i + (int)(pop / 3)];
		int r3 = randomPop[i + 2 * (int)(pop / 3)];
		int winner = compete1(r1, r2, r3);
		int loser = compete2(r1, r2, r3);
		int middle = r1 + r2 + r3 - winner - loser;
		update(winner, middle, loser);
	}
}

int OBL_CPSO::compete1(int r1, int r2, int r3) {
	if ((*population)[r1]->getObj(0) < (*population)[r2]->getObj(0)
	    && (*population)[r1]->getObj(0) < (*population)[r3]->getObj(0))
		return r1;
	if ((*population)[r2]->getObj(0) < (*population)[r1]->getObj(0)
	    && (*population)[r2]->getObj(0) < (*population)[r3]->getObj(0))
		return r2;
	if ((*population)[r3]->getObj(0) < (*population)[r1]->getObj(0)
	    && (*population)[r3]->getObj(0) < (*population)[r2]->getObj(0))
		return r3;
	if ((*population)[r1]->getObj(0) == (*population)[r2]->getObj(0))
		return r1;
	if ((*population)[r1]->getObj(0) == (*population)[r3]->getObj(0))
		return r1;
	if ((*population)[r2]->getObj(0) == (*population)[r3]->getObj(0))
		return r2;
}

/**
 * Determina el peor de todos, el de OBJ mayor
 */
int OBL_CPSO::compete2(int r1, int r2, int r3) {
	if ((*population)[r1]->getObj(0) > (*population)[r2]->getObj(0)
	    && (*population)[r1]->getObj(0) > (*population)[r3]->getObj(0))
		return r1;
	if ((*population)[r2]->getObj(0) > (*population)[r1]->getObj(0)
	    && (*population)[r2]->getObj(0) > (*population)[r3]->getObj(0))
		return r2;
	if ((*population)[r3]->getObj(0) > (*population)[r1]->getObj(0)
	    && (*population)[r3]->getObj(0) > (*population)[r2]->getObj(0))
		return r3;
	if ((*population)[r1]->getObj(0) == (*population)[r2]->getObj(0))
		return r2;
	if ((*population)[r1]->getObj(0) == (*population)[r3]->getObj(0))
		return r3;
	if ((*population)[r2]->getObj(0) == (*population)[r3]->getObj(0))
		return r3;
}

void OBL_CPSO::update(int winner, int middle, int loser) {
	double coeff1, coeff2, coeff3, coeff4;
	inertia = (0.7 - 0.2) * (getEvaluations() - getPerformedEvaluations())
	          / getEvaluations() + 0.2;
	const int dimension = (*population)[0]->getNumberOfVar();
	random_device rd;
	mt19937 mt(rd());
	uniform_real_distribution<double> dist(0.0, std::nextafter(1.0, DBL_MAX));
	for (int j = 0; j <  dimension; j++) {
		coeff1 = dist(mt);
		coeff2 = dist(mt);
		coeff3 = dist(mt);
		coeff4 = dist(mt);
		velocity[loser][j] = coeff1 * velocity[loser][j] + coeff2
		                     * ((*population)[winner]->getVar(j)
		                        - (*population)[loser]->getVar(j))
		                     + coeff3 * getInertia() * (mean->getVar(j)
		                         - (*population)[loser]->getVar(j));
		(*population)[loser]->setVar(j, (*population)[loser]->getVar(j)
		                             + velocity[loser][j]);
		if ((*population)[loser]->getVar(j) > (*population)[loser]->getMaximum(j))
			(*population)[loser]->setVar(j, (*population)[loser]->getMaximum(j));
		if ((*population)[loser]->getVar(j) < (*population)[loser]->getMinimum(j))
			(*population)[loser]->setVar(j, (*population)[loser]->getMinimum(j));
		(*population)[middle]->setVar(j, -(*population)[middle]->getVar(j)
		                              + coeff4 * (*population)[middle]->getVar(j));
		if ((*population)[middle]->getVar(j)
		    > (*population)[middle]->getMaximum(j))
			(*population)[middle]->setVar(j, (*population)[middle]->getMaximum(j));
		if ((*population)[middle]->getVar(j)
		    < (*population)[middle]->getMinimum(j))
			(*population)[middle]->setVar(j, (*population)[middle]->getMinimum(j));
	}
	evaluate((*population)[loser]);
	evaluate((*population)[middle]);
	// Comprobamos si hay mejora con respecto a best
	if ((*population)[loser]->getObj(0) < best->getObj(0))
		best = (*population)[loser]->internalClone();
	if ((*population)[middle]->getObj(0) < best->getObj(0))
		best = (*population)[middle]->internalClone();
}

void OBL_CPSO::globalSearch() {
#ifdef DEBUG
	cout << "Global-Search" << endl;
#endif
	sort(population->begin(), population->end(), sortByObj0);
	Individual* best = (*population)[0];
	vector<bool> explored (getPopulationSize(), false);
	int numberImp;
	int numberExp = 0;
	// Calculates the centroid of the current population
	Individual* centroid = getSampleInd()->internalClone();
	for (int i = 0; i < centroid->getNumberOfVar(); i++) {
		double sum = 0;
		for (int j = 0; j < getPopulationSize(); j++) {
			sum += (*population)[j]->getVar(i);
		}
		centroid->setVar(i, sum / getPopulationSize());
	}
	do {
		bool improvement = false;
		numberImp = 0;
		int k;
		do {
			k = getRandomInteger0_N(getPopulationSize() - 1);
		} while ((explored[k]) && (numberExp < getPopulationSize()));
		explored[k] = true;
		numberExp++;
		do {
			// Three random numbers (a1, a2, a3) belonging to the range [0, 1] are obtained
			// The condition a1 + a2 + a3 == 1 must be satisfied
			double a1, a2, a3;
			do {
				a1 = ((double) rand () / (double) RAND_MAX);
				a2 = ((double) rand () / (double) RAND_MAX);
				a3 = 1.0 - a1 - a2;
			} while  (abs(a1) + abs(a2) + abs(a3) != 1.0);
			int r1;
			// If the population size is 1, this might not work!!!!
			do {
				r1 = getRandomInteger0_N(getPopulationSize() - 1);
			} while (r1 == k);
			// Applies global neighbourhood search strategy considering the centroid of the population
			Individual* v = getSampleInd()->internalClone();
			for (int i = 0; i < v->getNumberOfVar(); i++) {
				v->setVar(i, a1 * (*population)[k]->getVar(i)
				          + a2 * centroid->getVar(i)
				          + a3 * (best->getVar(i) - (*population)[r1]->getVar(i)));
				// Checks lower and upper limits of variables
				if ((v->getVar(i) < v->getMinimum(i))
				    || (v->getVar(i) > v->getMaximum(i))) {
					double r = ((double) rand () / (double) RAND_MAX);
					v->setVar(i, r * (v->getMaximum(i) - v->getMinimum(i)) + v->getMinimum(i));
				}
				v->setVar(i, max(v->getVar(i), v->getMinimum(i)));
				v->setVar(i, min(v->getVar(i), v->getMaximum(i)));
			}
			// Evaluates the new individual
			evaluate(v);
			// Selects the best individual between v and k
			if (((v->getOptDirection(0) == MINIMIZE)
			     && (v->getObj(0) < (*population)[k]->getObj(0)))
			    || ((v->getOptDirection(0) == MAXIMIZE)
			        && (v->getObj(0) > (*population)[k]->getObj(0)))) {
				//delete (*population)[k];
				population->push_back((*population)[k]);
				(*population)[k] = v;
				improvement = true;
				numberImp++;
				// Checks if the new individual v is the best one
				if ((v->getOptDirection(0) == MINIMIZE)
				    && (v->getObj(0) < best->getObj(0)))
					best = v;
				else if ((v->getOptDirection(0) == MAXIMIZE)
				         && (v->getObj(0) > best->getObj(0)))
					best = v;
			} else {
				delete v;
				improvement = false;
			}
		} while (improvement);
	} while ((numberImp > 0) && (numberExp < getPopulationSize()));
	// Selects the best individuals from the current population and the global neighbourhood path
	sort(population->begin(), population->end(), sortByObj0);
	for (int i = 0; i < (population->size() - getPopulationSize()); i++) {
		delete (*population)[population->size() - 1];
		population->pop_back();
	}
}

void OBL_CPSO::fillPopWithNewIndsAndEvaluate() {
	const int type = getTypeOfInit();
	if (type != RAND && type != OBL && type != QOBL && type != QROBL) {
		cerr << "Error al especificar el tipo de inicialización" << endl;
		cerr << "Se empleará una inicialización aleatoria" << endl;
	}
	if (type == RAND)
		fillRandom();
	else if (type == OBL)
		oppositionBasedLearning();
	else if (type == QOBL)
		quasiOppositionBasedLearning();
	else
		quasiReflectedOppositionBasedLearning();
#ifdef DEBUG
	cout << "Poblacion creada" << endl;
	for (int i = 0; i < getPopulationSize(); i++) {
		for (int j = 0; j < (*population)[i]->getNumberOfVar(); j++) {
			cout << "Ind[" << i << ", " << j << "]: "
			     << (*population)[i]->getVar(j) << endl;
		}
		cout << endl;
	}
#endif
#ifdef WAIT
	getchar();
#endif
}

void OBL_CPSO::fillRandom() {
	for (int i = population->size(); i < getPopulationSize(); i++) {
		Individual* ind = getSampleInd()->internalClone();
		ind->restart();
		evaluate(ind);
		population->push_back(ind);
	}
}

void OBL_CPSO::oppositionBasedLearning() {
	//mt19937 generator(seed);
	// Generates a total number of individuals equal to (getPopulationSize() - population->size()) * 2
	vector<Individual*> genPopulation;
	for (int i = 0; i < (getPopulationSize() - population->size()); i++) {
		Individual* ind = getSampleInd()->internalClone();
		Individual* opp_ind = getSampleInd()->internalClone();
		for (int j = 0; j < ind->getNumberOfVar(); j++) {
			// Generates a random population
			//double aux = (double) generator() / (double) generator.max();
			double aux = ((double) rand () / (double) RAND_MAX);
			ind->setVar(j, aux * (ind->getMaximum(j) - ind->getMinimum(j))
			            + ind->getMinimum(j));
			// Generates the oppositional population
			opp_ind->setVar(j, opp_ind->getMinimum(j) + opp_ind->getMaximum(j)
			                - ind->getVar(j));
		}
		evaluate(ind);
		evaluate(opp_ind);
		genPopulation.push_back(ind);
		genPopulation.push_back(opp_ind);
	}
	sort(genPopulation.begin(), genPopulation.end(), sortByObj0);
	for (int i = 0; i < genPopulation.size() / 2; i++) {
		population->push_back(genPopulation[i]);
		delete genPopulation[genPopulation.size() - 1 - i];
	}
	genPopulation.clear();
}

void OBL_CPSO::quasiOppositionBasedLearning() {
	//mt19937 generator(seed);
	// Generates a total number of individuals equal to (getPopulationSize() - population->size()) * 2
	vector<Individual*> genPopulation;
	for (int i = 0; i < (getPopulationSize() - population->size()); i++) {
		Individual* ind = getSampleInd()->internalClone();
		Individual* qopp_ind = getSampleInd()->internalClone();
		for (int j = 0; j < ind->getNumberOfVar(); j++) {
			// Generates a random population
			//double aux = (double) generator() / (double) generator.max();
			double aux = ((double) rand () / (double) RAND_MAX);
			ind->setVar(j, aux * (ind->getMaximum(j) - ind->getMinimum(j))
			            + ind->getMinimum(j));
			// Generates the quasi-oppositional population
			double opp = qopp_ind->getMinimum(j) + qopp_ind->getMaximum(j)
			             - ind->getVar(j);
			double m = (qopp_ind->getMinimum(j) + qopp_ind->getMaximum(j)) / 2.0;
			//double aux2 = (double) generator() / (double) generator.max();
			double aux2 = ((double) rand () / (double) RAND_MAX);
			if (ind->getVar(j) < m) {
				qopp_ind->setVar(j, m + (opp - m) * aux2);
			} else {
				qopp_ind->setVar(j, opp + (m - opp) * aux2);
			}
		}
		evaluate(ind);
		evaluate(qopp_ind);
		genPopulation.push_back(ind);
		genPopulation.push_back(qopp_ind);
	}
	sort(genPopulation.begin(), genPopulation.end(), sortByObj0);
	for (int i = 0; i < genPopulation.size() / 2; i++) {
		population->push_back(genPopulation[i]);
		delete genPopulation[genPopulation.size() - 1 - i];
	}
	genPopulation.clear();
}

void OBL_CPSO::quasiReflectedOppositionBasedLearning() {
	//mt19937 generator(seed);
	// Generates a total number of individuals equal to (getPopulationSize() - population->size()) * 2
	vector<Individual*> genPopulation;
	for (int i = 0; i < (getPopulationSize() - population->size()); i++) {
		Individual* ind = getSampleInd()->internalClone();
		Individual* qropp_ind = getSampleInd()->internalClone();
		for (int j = 0; j < ind->getNumberOfVar(); j++) {
			// Generates a random population
			//double aux = (double) generator() / (double) generator.max();
			double aux = ((double) rand () / (double) RAND_MAX);
			ind->setVar(j, aux * (ind->getMaximum(j) - ind->getMinimum(j))
			            + ind->getMinimum(j));
			// Generates the quasi-reflected-oppositional population
			double m = (qropp_ind->getMinimum(j) + qropp_ind->getMaximum(j)) / 2.0;
			//double aux2 = (double) generator() / (double) generator.max();
			double aux2 = ((double) rand () / (double) RAND_MAX);
			if (ind->getVar(j) < m) {
				qropp_ind->setVar(j, ind->getVar(j) + (m - ind->getVar(j)) * aux2);
			} else {
				qropp_ind->setVar(j, m + (ind->getVar(j) - m) * aux2);
			}
		}
		evaluate(ind);
		evaluate(qropp_ind);
		genPopulation.push_back(ind);
		genPopulation.push_back(qropp_ind);
	}
	sort(genPopulation.begin(), genPopulation.end(), sortByObj0);
	for (int i = 0; i < genPopulation.size() / 2; i++) {
		population->push_back(genPopulation[i]);
		delete genPopulation[genPopulation.size() - 1 - i];
	}
	genPopulation.clear();
}

void OBL_CPSO::sortPopulation() {
	sort(population->begin(), population->end(), sortByObj0);
}

void OBL_CPSO::getSolution(MOFront* p) {
	sortPopulation();
	p->insert((*population)[0]);
}

void OBL_CPSO::printInfo(ostream& os) const {
	os << "Particle Swarm Algorithm"  << endl;
	os << "Number of Evaluations = " << getEvaluations() << endl;
	os << "Population Size = " << getPopulationSize() << endl;
}

void OBL_CPSO::readParam(map<string, string>& args, const string param, int& value) {
	if (args.count(param) == ZERO) {
		cerr << "No se ha especificado correctamente el parametro " << param << endl;
		exit(-1);
	}
	value = atoi(args[param].c_str());
	args.erase(param);
}

void OBL_CPSO::readParam(map<string, string>& args, const string param, double& value) {
	if (args.count(param) == ZERO) {
		cerr << "No se ha especificado correctamente el parametro " << param << endl;
		exit(-1);
	}
	value = atof(args[param].c_str());
	args.erase(param);
}

void OBL_CPSO::readParam(map<string, string>& args, const string param, string& value) {
	if (args.count(param) == ZERO) {
		cerr << "No se ha especificado correctamente el parametro " << param << endl;
		exit(-1);
	}
	value = args[param];
	args.erase(param);
}

bool sortByObj0(const Individual* i1, const Individual* i2) {
	return (i1->getInternalOptDirection(0) == MINIMIZE) ?
	       (i1->getObj(0) < i2->getObj(0)) : (i1->getObj(0) > i2->getObj(0));
}

const int OBL_CPSO::PARAMS = 3;
const int OBL_CPSO::DEFAULT_W = 0.5;
const int OBL_CPSO::DEFAULT_POP = 150;
const int OBL_CPSO::RAND = 0;
const int OBL_CPSO::OBL = 1;
const int OBL_CPSO::QOBL = 2;
const int OBL_CPSO::QROBL = 3;