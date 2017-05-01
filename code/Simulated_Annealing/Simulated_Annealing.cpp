#include <cstdlib>
#include <iostream>
#include <signal.h>
#include <math.h>
#include "Simulated_Annealing.h"

#define DEBUG

// Algorithm initialization
bool Simulated_Annealing::init(const vector<string>& params) {
	// Check number of parameters
	if (params.size() < PARAMS) {
		cout << "Parametros: pop-size seed INIT init_temp [globalSearch]" << endl;
		cout << "INIT" << endl;
		cout << "\t RAND" << endl;
		cout << "\t OBL" << endl;
		cout << "\t QOBL" << endl;
		cout << "\t QROBL" << endl;
		return false;
	}
// Only mono-objective optimization is supported
	if (getSampleInd()->getNumberOfObj() != 1) {
		cout << "Multi-Objective not supported" << endl;
		return false;
	}
	search = false;
	setPopulationSize(atoi(params[0].c_str()));
	setTemp(atof(params[2].c_str()));
	seed = static_cast<uint64_t>(atof(params[1].c_str()));
	srand(seed);
	initRandomSeed(seed);
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
	setTemp(atof(params[3].c_str()));
	setTempVariation(DIFF);
	if (params.size() == PARAMS + 1 && params[PARAMS + 1] == "globalSearch") {
		search = true;
	}
#ifdef DEBUG
	cout << "Seed: " << seed << endl;
	cout << "Population size: " << getPopulationSize() << endl;
	cout << "Init_Temp: " << getTemp() << endl;
	cout << "INIT: " << getTypeOfInit();
	cout << " (RAND: " << RAND << " OBL: " << OBL
	     << " QOBL: " << QOBL << " QROBL: " << QROBL << ")" << endl;
#endif
	return true;
}

void Simulated_Annealing::runGeneration() {
	applyRandomPerturbations();
	evaluateDifference();
	updateTemperature();
	if (isGlobalSearch())
		globalSearch();
}

void Simulated_Annealing::applyRandomPerturbations() {
	offsprings.clear();
	offsprings.resize(getPopulationSize());
	for (int i = 0; i < getPopulationSize(); i++) {
		perturbation = ((double) rand () / (double) RAND_MAX);
		Individual* ind = (*population)[i]->internalClone();
		for (int j = 0; j < (*population)[i]->getNumberOfVar(); j++) {
			ind->setVar(j, getPerturbation() * (ind->getMaximum(j)
			                                    - ind->getMinimum(j))
			            + ind->getMinimum(j));
			if (ind->getVar(j) > ind->getMaximum(j)
			    || ind->getVar(j) < ind->getMinimum(j)) {
				ind->setVar(j, (ind->getMaximum(j) - ind->getMinimum(j))
				            + ind->getMinimum(j));
			}
		}
		offsprings[i] = ind;
	}
}

void Simulated_Annealing::evaluateDifference() {
	for (int i = 0; i < getPopulationSize(); i++) {
		evaluate((*population)[i]);
		evaluate(offsprings[i]);
		double difference = 0;
		difference = (*population)[i]->getObj(0) - offsprings[i]->getObj(0);
		if (isgreater((*population)[i]->getObj(0),
		              offsprings[i]->getObj(0))) {
			(*population)[i] = offsprings[i]->internalClone();
		} else {
			const double randomProb = (double) rand() / (RAND_MAX);
			const double probability = exp(((*population)[i]->getObj(0) -
			                                offsprings[i]->getObj(0)) / getTemp());
			if (isgreater(randomProb, probability)) {
				(*population)[i] = offsprings[i]->internalClone();
			}
		}
	}
}

void Simulated_Annealing::updateTemperature() {
	setTemp(getTemp() * getTempVariation());
}

void Simulated_Annealing::globalSearch() {
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

void Simulated_Annealing::fillPopWithNewIndsAndEvaluate() {
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
}

void Simulated_Annealing::fillRandom() {
	for (int i = population->size(); i < getPopulationSize(); i++) {
		Individual* ind = getSampleInd()->internalClone();
		ind->restart();
		evaluate(ind);
		population->push_back(ind);
	}
}

void Simulated_Annealing::oppositionBasedLearning() {
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

void Simulated_Annealing::quasiOppositionBasedLearning() {
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

void Simulated_Annealing::quasiReflectedOppositionBasedLearning() {
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

void Simulated_Annealing::sortPopulation() {
	sort(population->begin(), population->end(), sortByObj0);
}

void Simulated_Annealing::getSolution(MOFront* p) {
	sortPopulation();
	p->insert((*population)[0]);
}

void Simulated_Annealing::printInfo(ostream& os) const {
	os << "Particle Swarm Algorithm"  << endl;
	os << "Number of Evaluations = " << getEvaluations() << endl;
	os << "Population Size = " << getPopulationSize() << endl;
}

void Simulated_Annealing::readParam(map<string, string>& args, const string param, int& value) {
	if (args.count(param) == ZERO) {
		cerr << "No se ha especificado correctamente el parametro " << param << endl;
		exit(-1);
	}
	value = atoi(args[param].c_str());
	args.erase(param);
}

void Simulated_Annealing::readParam(map<string, string>& args, const string param, double& value) {
	if (args.count(param) == ZERO) {
		cerr << "No se ha especificado correctamente el parametro " << param << endl;
		exit(-1);
	}
	value = atof(args[param].c_str());
	args.erase(param);
}

void Simulated_Annealing::readParam(map<string, string>& args, const string param, string& value) {
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


const int Simulated_Annealing::PARAMS = 4;
const int Simulated_Annealing::RAND = 0;
const int Simulated_Annealing::OBL = 1;
const int Simulated_Annealing::QOBL = 2;
const int Simulated_Annealing::QROBL = 3;
const double Simulated_Annealing::DIFF = 0.9;
const int Simulated_Annealing::DEFAULT_POP = 25;