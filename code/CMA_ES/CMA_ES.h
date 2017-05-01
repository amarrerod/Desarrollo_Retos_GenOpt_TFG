
/**
 * Project: TFG
 * File: oplink/algorithms/team/src/plugins/algorithms/CMA_ES/CMA_ES.h
 * Author: Alejandro Marrero Díaz
 * Contact: alu0100825008@ull.edu.es
 * Date: 2017-02-24 13:13:33
 * Description:  Covariance Matrix Adaptation Evolutionary Strategy
  * GENERAL PURPOSE: The CMA-ES (Evolution Strategy with Covariance
  * Matrix Adaptation) is a robust search method which should be
  * applied, if derivative based methods, e.g. quasi-Newton BFGS or
  * conjucate gradient, (supposably) fail due to a rugged search
  * landscape (e.g. noise, local optima, outlier, etc.). On smooth
  * landscapes CMA-ES is roughly ten times slower than BFGS. For up to
  * N=10 variables even the simplex direct search method (Nelder & Mead)
  * is often faster, but far less robust than CMA-ES.  To see the
  * advantage of the CMA, it will usually take at least 30*N and up to
  * 300*N function evaluations, where N is the search problem dimension.
  * On considerably hard problems the complete search (a single run) is
  * expected to take at least 30*N^2 and up to 300*N^2 function
  * evaluations.
  * The strategy parameter lambda (population size, opts.PopSize) is the
  * preferred strategy parameter to play with.  If results with the
  * default strategy are not satisfactory, increase the population
  * size. (Remark that the crucial parameter mu (parentNumber) is
  * increased proportionally to lambda). This will improve the
  * strategies capability of handling noise and local minima. We
  * recomment successively increasing lambda by a factor of about three,
  * starting with initial values between 5 and 20. Casually, population
  * sizes even beyond 1000+100*N can be sensible.
  **/
#ifndef __CMA_ES_H__
#define __CMA_ES_H__

#include "Individual.h"
#include "EA.h"
#include <vector>
#include <memory>
#include <utility>
#include <Eigen/Dense>
#include "eigenmvn.h"
#include <memory>

using namespace std;
using namespace Eigen;

typedef Matrix<double, Dynamic, Dynamic> matrix;
typedef Matrix<double, Dynamic, 1> vectorx;
typedef Matrix<int32_t, Dynamic, 1> ivectorx;

class CMA_ES : public EA {
 public:
	inline void setNumberOfParents(const int p) { numberOfParents = p; };
	inline const int getNumberOfParents() { return numberOfParents; };
	
	inline void setSigma(const double s) { sigma = s; };
	inline const double getSigma() { return sigma; };
	
	inline void setDimension(const int d) { dimension = d; };
	inline int getDimension() { return dimension; };
 public:
	bool init(const vector<string>& params);
	void getSolution(MOFront* p);
	void printInfo(ostream& os) const;
 private:
	void initParams();
	void runGeneration();
	void scale(const vectorx& x, vectorx& scaled);
	void sortPopulation();
	void sampling();
	void rankAndSort();
	void assignNewMean();
	void updateEvolutionPaths();
	void updateWeights();
	void updateStepSize();
	void updateCovarianceMatrix();
	void eigenDecomposition();
	void stoppingCriteria();
	void updateBest();
	void restart();
	bool checkRestart();
	void globalSearch();
 private:
	int dimension;  // Tamaño del problema
	int numberOfParents;    // AKA Mu
	double sigma;
	double sigmaInit;
	uint64_t seed;
	int maxRestarts;
	int restarts;
	bool restarted;
	double fBest;
	vectorx best;
	mt19937 mt;
	normal_distribution<double> dist_normal_real;
	uniform_real_distribution<double> dist_uniform_real;
	EigenSolver<matrix> eigenSolver;
	// Candidatos
	vectorx scaled;
	matrix rankedParents;
	ivectorx keys; // indices de los candidatos
	vectorx objOffSprings;
	matrix xOffSprings;
	matrix zOffSprings;
	matrix yOffSprings;
	matrix yOffSpringsRanked;
	vectorx weights;
	vectorx variableWeights;
	vectorx yMean;
	vectorx xMean;
	vectorx xMeanOld;
	double cc;
	double cs;
	double c1;
	double cmu;
	double chi;
	double ds;
	double muEffective;
	vectorx ps;
	vectorx pc;
	double preCalculatedCFact;
	double preCalculatedSFact;
	vectorx eigenValuesC;
	matrix C;
	matrix CInvSqrt;
	matrix B;
	matrix D;
	bool hSig; // Heaviside
 private:
	static const double DEF_SIGMA;
	static const int DEF_LAMBDA;
	static const int MIN_ARGS;
	static const int MAX_ARGS;
	static const int MAX_IT;  // Nº de iteraciones sin mejora
};

bool sortByObj0(const Individual* i1, const Individual* i2);

#endif
