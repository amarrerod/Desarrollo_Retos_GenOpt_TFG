


#ifndef _Simulated_Annealing_H_
#define _Simulated_Annealing_H_

#include "Individual.h"
#include "EA.h"
#include <vector>
#include <map>
#include <random>
#include <cmath>

#define ZERO 0

using namespace std;

class Simulated_Annealing : public EA {
 public:
	bool init(const vector<string>& params);
	void getSolution(MOFront* p);
	void printInfo(ostream& os) const;
	void setTypeOfInit(const int i) { typeOfInit = i; };
	void setGlobalSearch(bool s) { search = s; };
	void setSeed(const uint64_t s) { seed = s; };
	void setTemp(const double temp) { this->temp = temp; };
	void setTempVariation(const double temp) { tempVariation = temp; };
	void setPerturbation(const double pert) { perturbation = pert; };
	double getPerturbation() { return perturbation; };
	double getTempVariation() { return tempVariation; };
	int getTypeOfInit() { return typeOfInit; };
	uint64_t getSeed() { return seed; };
	bool isGlobalSearch() { return search; };
	double getTemp() { return temp; };
 private:
	// Funciones de inicializacion de parametros
	void readParam(map<string, string>&, const string, string&);
	void readParam(map<string, string>&, const string, double&);
	void readParam(map<string, string>&, const string, int&);
	
	void runGeneration();
	void sortPopulation();
	void applyRandomPerturbations();
	void evaluateDifference();
	void updateTemperature();
	
	// Funciones de inicializacion de la poblacion
	void fillPopWithNewIndsAndEvaluate();
	void oppositionBasedLearning();
	void quasiOppositionBasedLearning();
	void quasiReflectedOppositionBasedLearning();
	void fillRandom();
	
	// Búsqueda global
	void globalSearch();
	
	bool search;
	int typeOfInit;
	double temp;
	double perturbation;
	double tempVariation;
	vector<Individual*> offsprings;
	// Semillar para inicializar la población
	uint64_t seed;
	
	// Constantes
	static const int PARAMS;
	static const int OBL;
	static const int QOBL;
	static const int QROBL;
	static const int RAND;
	static const double TEMP;
	static const double DIFF;
	static const int DEFAULT_POP;
};

bool sortByObj0(const Individual* i1, const Individual* i2);

#endif
