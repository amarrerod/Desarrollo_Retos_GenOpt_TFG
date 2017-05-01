


#ifndef _OBL_PSO_H_
#define _OBL_PSO_H_

#include "Individual.h"
#include "EA.h"
#include <vector>
#include <map>
#include <random>
#include <cmath>

#define ZERO 0

using namespace std;

class OBL_CPSO : public EA {
 public:
	bool init(const vector<string>& params);
	void getSolution(MOFront* p);
	void printInfo(ostream& os) const;
	void setInertia(const double w) { inertia = w; };
	void setTypeOfInit(const int i) { typeOfInit = i; };
	void setGlobalSearch(bool s) { search = s; };
	void setSeed(const uint64_t s) { seed = s; };
	int getTypeOfInit() { return typeOfInit; };
	double getInertia() { return inertia; };
	uint64_t getSeed() { return seed; };
	bool isGlobalSearch() { return search; };
 private:
	// Funciones de inicializacion de parametros
	void readParam(map<string, string>&, const string, string&);
	void readParam(map<string, string>&, const string, double&);
	void readParam(map<string, string>&, const string, int&);
	
	void runGeneration();
	void initVelocity();
	void evaluateSwarm();
	void createRandomPop();
	void createMean();
	void compete();
	int compete1(int, int, int);
	int compete2(int, int, int);
	void update(int, int, int);
	void sortPopulation();
	
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
	// Semillar para inicializar la población
	uint64_t seed;
	// Mejor valor obtenido entre todas las partículas
	Individual* best;
	Individual* mean;
	double inertia;
	vector<int> randomPop;
	vector<vector<double>> velocity;
	vector<vector<double>> position;
	// Constantes
	static const int PARAMS;
	static const int DEFAULT_POP;
	static const int DEFAULT_W;
	static const int OBL;
	static const int QOBL;
	static const int QROBL;
	static const int RAND;
};

bool sortByObj0(const Individual* i1, const Individual* i2);

#endif
