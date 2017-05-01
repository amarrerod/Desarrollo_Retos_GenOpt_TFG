/*
* @Author: Alejandro Marrero
* @Date:   2017-03-16 16:27:00
* @Last Modified by:   Alejandro Marrero
* @Last Modified time: 2017-03-16 17:48:44
*/

// SOLO DIOS SABE COMO HE CONSEGUIDO QUE ESTO FUNCIONE

#include <cstdlib>
#include <iostream>
#include <signal.h>
#include <math.h>
#include "CMA_ES.h"
#include <algorithm>
#include <chrono>

//#define DEBUG
//#define LOG
//#define WAIT
//#define COUNT

// Algorithm initialization
bool CMA_ES::init(const vector<string>& params) {
	// Check number of parameters
	if (params.size() < MIN_ARGS || params.size() > MAX_ARGS) {
		cout << "Parametros: seed restarts sigma [lambda]" << endl;
		return false;
	}
	// Only mono-objective optimization is supported
	if (getSampleInd()->getNumberOfObj() != 1) {
		cout << "Multi-Objective not supported" << endl;
		return false;
	}
	setDimension(getSampleInd()->getNumberOfVar());
	const int dim = getDimension();
	seed = static_cast<uint64_t>(atof(params[0].c_str()));
	maxRestarts = static_cast<int>(atoi(params[1].c_str()));
	sigmaInit = static_cast<double>(atof(params[2].c_str()));
	restarts = 0;
	// Tamaño de la poblacion
	if (params.size() == MAX_ARGS) {
		int lambda = atoi(params[params.size() - 1].c_str());
		setPopulationSize(lambda);
	} else {
		const int lambda = 4 + floor(3 * log(dim));
		setPopulationSize(lambda);
		cout << "Default popSize based on problem dimension: " << lambda << endl;
		getchar();
	}
	mt.seed(seed);
	sigmaInit = 0.3;
	setNumberOfParents(floor(getPopulationSize() / 2));
	initParams();
	restarted = false;
#ifdef DEBUG
	cout << "init done" << endl;
#endif
	return true;
}

// CMA_ES generation
void CMA_ES::runGeneration() {
#ifdef COUNT
	cout << "Generacion: " << getGeneration() << endl;
	cout << "Evaluaciones: " << getPerformedEvaluations() << endl;
#endif
	sampling();
	sortPopulation();
	rankAndSort();
	updateBest();
	assignNewMean();
	updateEvolutionPaths();
	updateWeights();
	updateCovarianceMatrix();
	eigenDecomposition();
	updateStepSize();
	globalSearch();
	if (checkRestart() && restarts < maxRestarts) {
		restart();
		restarts++;
	}
#ifdef DEBUG
	cout << "Generation done" << endl;
#endif
}

/**
 * Creamos la población y luego la evaluamos
 * Con cada operacion comprobamos que los resultados
 * no se salen de los rangos establecidos para cada
 * una de las variables
 */
void CMA_ES::sampling() {
#ifdef DEBUG
	cout << "Starting sampling!" << endl;
#endif
	const int dim = getDimension();
	for (int j = 0; j < getPopulationSize(); j++) {
		for (int i = 0; i < dim; i++) {
			zOffSprings(i, j) = dist_normal_real(mt);
		}
	}
	matrix bMulD(dim, dim);
	bMulD.block(0, 0, dim, dim) = B * D.block(0, 0, dim, dim);
	yOffSprings.block(0, 0, dim, getPopulationSize()) = bMulD
	    * zOffSprings.block(0, 0, dim, getPopulationSize());
	for (int j = 0; j < getPopulationSize(); j++) {
		xOffSprings.col(j) = xMean + sigma * yOffSprings.col(j);
		scale(xOffSprings.col(j), scaled);
		xOffSprings.col(j) = scaled;
		for (int i = 0; i < dim; i++) {
			(*population)[j]->setVar(i, xOffSprings(i, j));
		}
		evaluate((*population)[j]);
		const double obj = (*population)[j]->getObj(0);
		objOffSprings[j] = isnan(obj) ? numeric_limits<double>::infinity() : obj;
	}
#ifdef DEBUG
	cout << "sampling done" << endl;
#endif
}

void CMA_ES::scale(const vectorx& x, vectorx& scale) {
	for (int i = 0; i < getDimension(); i++) {
		scale[i] = x[i];
		if (x[i] < getSampleInd()->getMinimum(i)
		    || x[i] > getSampleInd()->getMaximum(i)
		    || isnan(x[i])
		    || isinf(x[i])) {
			double randAux = ((double) rand () / (double) RAND_MAX);
			scale[i] = randAux * ((*population)[0]->getMaximum(i)
			                      - (*population)[0]->getMinimum(i))
			           + (*population)[0]->getMinimum(i);
		}
	}
#ifdef DEBUG
	cout << "scaling done" << endl;
#endif
}

void CMA_ES::rankAndSort() {
	iota(keys.data(), keys.data() + getPopulationSize(), 0);
	sort(keys.data(), keys.data() + getPopulationSize(),
	[&](int i, int j) {	return objOffSprings[i] < objOffSprings[j]; });
	for (int i = 0; i < getNumberOfParents(); i++) {
		int key = keys[i];
		rankedParents.col(i) = xOffSprings.col(key);
	}
	for (int i = 0; i < getPopulationSize(); i++) {
		int key = keys[i];
		yOffSpringsRanked.col(i) = yOffSprings.col(key);
	}
#ifdef DEBUG
	cout << "Rank and sort done" << endl;
#endif
}

void CMA_ES::updateBest() {
	double obj = (*population)[0]->getObj(0);
	if (!isnan(obj) && obj < fBest) {
		best = rankedParents.col(0);
		fBest = obj;
	}
}

void CMA_ES::assignNewMean() {
#ifdef DEBUG
	cout << "Starting assignNewMean" << endl;
#endif
	xMeanOld = xMean;
	yMean = yOffSpringsRanked.block(0, 0, getDimension(), getNumberOfParents())
	        * weights.block(0, 0, getNumberOfParents(), 1);
	xMean = xMean + sigma * yMean;
#ifdef DEBUG
	cout << "assignNewMean done" << endl;
#endif
}

void CMA_ES::updateEvolutionPaths() {
#ifdef DEBUG
	cout << "Starting updateEvolutionPaths" << endl;
#endif
	ps = (1.0 - cs) * ps + preCalculatedSFact * CInvSqrt * yMean;
	const double psNorm = ps.norm();
	const double threshold = (1.4 + 2.0 / (getDimension() + 1)) * sqrt(1.0 - pow(1.0 * cs, 2.0 * (getGeneration() + 1))) * chi;
	hSig = psNorm < threshold;
	pc = (1.0 - cc) * pc + hSig * preCalculatedCFact * yMean;
#ifdef DEBUG
	cout << "updateEvolutionPaths done" << endl;
#endif
}

void CMA_ES::updateWeights() {
#ifdef DEBUG
	cout << "Starting updateWeights" << endl;
#endif
	const int dim = getDimension();
	for (int i = 0; i < getPopulationSize(); i++) {
		if (weights[i] < 0) {
			variableWeights[i] = weights[i] * dim / (CInvSqrt
			                     * yOffSpringsRanked.col(i)).squaredNorm();
		}
	}
#ifdef DEBUG
	cout << "updateWeights done" << endl;
#endif
}

void CMA_ES::updateCovarianceMatrix() {
#ifdef DEBUG
	cout << "Starting updateCovarianceMatrix" << endl;
#endif
	const int dim = getDimension();
	double h1 = (1 - hSig) * cc * (2.0 - cc);
	DiagonalMatrix<double, Dynamic, Dynamic> W(variableWeights.block(0, 0, getPopulationSize(), 1));
	double weightsSum = weights.block(0, 0, getPopulationSize(), 1).sum();
	C = (1.0 + c1 * h1 - c1 - cmu * weightsSum) * C + c1 * pc * pc.transpose()
	    + cmu * (yOffSpringsRanked.block(0, 0, dim, getPopulationSize()) * W * yOffSpringsRanked.block(0, 0, dim, getPopulationSize()).transpose());
#ifdef DEBUG
	cout << "updateCovarianceMatrix done" << endl;
#endif
}

void CMA_ES::updateStepSize() {
#ifdef DEBUG
	cout << "Starting updateStepSize" << endl;
#endif
	sigma *= exp(cs / ds * (ps.norm()) / chi - 1.0);
#ifdef DEBUG
	cout << "updateStepSize done" << endl;
#endif
}

void CMA_ES::eigenDecomposition() {
#ifdef DEBUG
	cout << "Starting eigenDecomposition" << endl;
#endif
	eigenSolver.compute(C);
	eigenValuesC = eigenSolver.eigenvalues().real();
	D = eigenValuesC.cwiseSqrt().asDiagonal();
	B = matrix(eigenSolver.eigenvectors().real());
	matrix dInv = eigenValuesC.cwiseSqrt().cwiseInverse().asDiagonal();
	CInvSqrt = B * dInv * B.transpose();
#ifdef DEBUG
	cout << "eigenDecomposition done" << endl;
#endif
}

void CMA_ES::globalSearch() {
#ifdef DEBUG
	cout << "Starting Global-Search" << endl;
#endif
	const int dim = getDimension();
	const int lambda = getPopulationSize();
	sortPopulation();
	vectorx objExtended(lambda);
	objExtended << objOffSprings;
	matrix xOffSpringsExtended(dim, lambda);
	xOffSpringsExtended << xOffSprings;
	int index = lambda;
	Individual* best = (*population)[0];
	vector<bool> explored (lambda, false);
	int numberImp;
	int numberExp = 0;
	// Calculates the centroid of the current population
	Individual* centroid = getSampleInd()->internalClone();
	for (int i = 0; i < dim; i++) {
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
			vectorx ind = vectorx::Zero(dim);
			for (int i = 0; i < v->getNumberOfVar(); i++) {
				v->setVar(i, a1 * (*population)[k]->getVar(i) + a2
				          * centroid->getVar(i)
				          + a3 * (best->getVar(i) - (*population)[r1]->getVar(i)));
				// Checks lower and upper limits of variables
				if ((v->getVar(i) < v->getMinimum(i))
				    || (v->getVar(i) > v->getMaximum(i))) {
					double r = ((double) rand () / (double) RAND_MAX);
					v->setVar(i, r * (v->getMaximum(i) - v->getMinimum(i)) + v->getMinimum(i));
				}
				v->setVar(i, max(v->getVar(i), v->getMinimum(i)));
				v->setVar(i, min(v->getVar(i), v->getMaximum(i)));
				ind[i] = v->getVar(i);
			}
			// Evaluates the new individual
			evaluate(v);
#ifdef DEBUG
			cout << "xOffSpringsExtended " << xOffSpringsExtended.rows()
			     << "x" << xOffSpringsExtended.cols() << endl;
			cout << "objExtended " << objExtended.rows() << "x"
			     << objExtended.cols() << endl;
			cout << "New index: " << index << endl;
#endif
			xOffSpringsExtended.resize(xOffSpringsExtended.rows(), xOffSpringsExtended.cols() + 1);
			xOffSpringsExtended.col(index - 1) = ind;
			objExtended.resize(objExtended.rows() + 1);
			objExtended[index - 1] =  v->getObj(0);
			index ++;
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
#ifdef DEBUG
	cout << "Global-Search almost done. Time to refill" << endl;
#endif
	sortPopulation();
	for (int i = 0; i < (population->size() - getPopulationSize()); i++) {
		delete (*population)[population->size() - 1];
		population->pop_back();
	}
#ifdef DEBUG
	cout << "population reset done" << endl;
#endif
	ivectorx keysExtended(objExtended.rows());
#ifdef DEBUG
	cout << "objExtended rows: " << objExtended.rows() << endl;
#endif
	iota(keysExtended.data(), keysExtended.data()
	     + (objExtended.rows()), 0);
#ifdef DEBUG
	cout << "iota done " << endl << keysExtended << endl;
#endif
	sort(keysExtended.data(), keysExtended.data() + (objExtended.rows()),
	[&](int i, int j) {	return objExtended[i] < objExtended[j]; });
#ifdef DEBUG
	cout << "sort done" << endl;
#endif
	for (int i = 0; i < getNumberOfParents(); i++) {
		int key = keysExtended[i];
		xOffSprings.col(i) = xOffSpringsExtended.col(key);
		rankedParents.col(i) = xOffSpringsExtended.col(key);
	}
#ifdef DEBUG
	cout << "globalSearch done" << endl;
#endif
}

bool CMA_ES::checkRestart() {
	const int dim = getDimension();
	double eigval_min = eigenValuesC[0];
	double eigval_max = eigenValuesC[dim - 1];
	// -> condition of covariance matrix
	if (eigval_max / eigval_min > 1e14) {
#ifdef DEBUG
		cout << "stopping criteria. generation: " << getGeneration()
		     << " reason: bad covariance condition." << std::endl;
#endif
		return true;
	}
	double sigma_fac = sigma / sigmaInit;
	double sigma_up_thresh = 1e12 * std::sqrt(eigval_max);
	if (sigma_fac / sigmaInit > sigma_up_thresh) {
#ifdef DEBUG
		std::cout << "stopping criteria. generation: " << getGeneration()
		          << " reason: sigma up." << std::endl;
#endif
		return true;
	}
	int nea = 0;
	for (int i = 0; i < dim; i++) {
		double ei = 0.1 * sigma * eigenValuesC[i];
		for (int j = 0; j < dim; j++) {
			nea += xMean[i] == xMean[i] + ei * B(j, i);
		}
	}
	if (nea > 0) {
#ifdef DEBUG
		std::cout << "stopping criteria. generation: " << getGeneration()
		          << " reason: no effect axis." << std::endl;
#endif
		return true;
	}
	int nec = 0;
	for (int i = 0; i < dim; i++) {
		nec += xMean[i] == xMean[i] + 0.2 * sigma * std::sqrt(C(i, i));
	}
	if (nec > 0) {
#ifdef DEBUG
		std::cout << "stopping criteria. generation: " << getGeneration()
		          << " reason: no effect coordinate." << std::endl;
#endif
		return true;
	}
	return false;
}


/**
 * Incrementamos el tamaño de la población y definimos un nuevo valor de sigma
 */
void CMA_ES::restart() {
#ifdef DEBUG
	cout << "restart" << endl;
	cout << "Sigma before restart: " << getSigma() << endl;
#endif
	double randExp = - 2 * randomUniform(0, 1);
	sigma = 2 * pow(10, randExp);
	const double fraction = 0.5 * (getPopulationSize() / DEF_LAMBDA);
	const double exp = randomUniform(0, 1);
	int lambda = floor(DEF_LAMBDA * (pow(fraction, exp * exp)));
	if (lambda == 0) {
		setPopulationSize(DEF_LAMBDA);
	} else {
		setPopulationSize(lambda);
	}
	setNumberOfParents(floor(getPopulationSize() / 2));
	population->clear();
	for (int i = 0; i < getPopulationSize(); i++) {
		Individual* ind = getSampleInd()->internalClone();
		population->push_back(ind);
	}
	initParams();
}

/**
 *  Inicialización de todos los parámetros del algoritmo
 *  dado que en init(vector<string>) no tenemos acceso a la
 *  dimensión del problema en cuestión
 */
void CMA_ES::initParams() {
#ifdef DEBUG
	cout << "initParams begins" << endl;
#endif
	const int dim = getDimension();
	const int lambda = getPopulationSize();
	// Definimos el tamaño para todas las matrices y vectores usados
	ps = vectorx::Zero(dim);
	pc = vectorx::Zero(dim);
	scaled = vectorx::Zero(dim);
	xMean = vectorx::Zero(dim);
	xMeanOld.resize(dim);
	yMean.resize(dim);
	eigenValuesC.resize(dim);
	weights.resize(lambda);
	objOffSprings.resize(lambda);
	keys.resize(lambda);
	variableWeights.resize(lambda);
	xOffSprings.resize(dim, lambda);
	zOffSprings.resize(dim, lambda);
	yOffSprings.resize(dim, lambda);
	yOffSpringsRanked.resize(dim, lambda);
	rankedParents.resize(dim, getNumberOfParents());
	C.resize(dim, dim);
	CInvSqrt.resize(dim, dim);
	B.resize(dim, dim);
	D.resize(dim, dim);
	best.resize(dim);
	Individual* ind = getSampleInd();
	for (int i = 0; i < dim; i++) {
		best[i] = ind->getVar(i);
	}
	fBest = ind->getObj(0);
	double weightsNegSum = 0.0, wPosSum = 0.0;
	for (int i = 0; i < lambda; i++) {
		weights[i] = log((lambda + 1.0) / 2.0) - log(i + 1);
		if (weights[i] >= 0)
			wPosSum += weights[i];
		else
			weightsNegSum += weights[i];
	}
	double sumParents = 0.0, weightsSqSumParents = 0.0;
	for (int i = 0; i < getNumberOfParents(); i++) {
		sumParents += weights[i];
		weightsSqSumParents += weights[i] * weights[i];
	}
	muEffective = sumParents * weightsSqSumParents / weightsSqSumParents;
	cs = (muEffective + 2.0) / (dim + muEffective + 5.0);
	cc = (4.0 + muEffective / dim) / (dim + 4.0 + 2.0 * muEffective / dim);
	c1 = 2.0 / (pow((dim + 1.3), 2) + muEffective);
	cmu = 2.0 * (muEffective - 2.0 + 1.0 / muEffective) / (pow(dim + 2.0, 2) + muEffective);
	ds = 1.0 + cs + 2.0 * max(0.0, sqrt((muEffective - 1) / (dim + 1)) - 1);
	chi = sqrt(dim) * (1.0 - 1.0 / (4.0 * dim) + 1.0 / (21.0 * dim * dim));
	preCalculatedCFact = sqrt(cc * (2.0 - cc) + muEffective);
	preCalculatedSFact = sqrt(cs * (2.0 - cs) * muEffective);
	// Active CMA_ES
	double aMu = 1.0 + c1 / cmu;
	double aMueff = 1.0 + 2 * muEffective;
	double aPosDefinited = (1.0 - cs - cmu) / (dim * cmu);
	double aMin = min(min(aMu, aMueff), aPosDefinited);
	for (int i = 0; i < lambda; i++) {
		weights[i] = weights[i] >= 0 ? weights[i] = weights[i] / wPosSum
		             : weights[i] = aMin * weights[i] / abs(weightsNegSum);
	}
	variableWeights = weights;
	sigma = sigmaInit;
	C.setIdentity();
	CInvSqrt.setIdentity();
	B.setIdentity();
	D.setIdentity();
#ifdef DEBUG
	cout << "Lambda: " << getPopulationSize() << endl;
	cout << "Mu: " << getNumberOfParents() << endl;
	cout << "cC: " << cc << endl;
	cout << "c1: " << c1 << endl;
	cout << "cMu: " << cmu << endl;
	cout << "damping: " << ds << endl;
	cout << "sigma: " << getSigma() << endl;
	cout << "initParams done!" << endl;
#endif
}

void CMA_ES::sortPopulation() {
	sort(population->begin(), population->end(), sortByObj0);
}

void CMA_ES::getSolution(MOFront* p) {
	sortPopulation();
	p->insert((*population)[0]);
}

void CMA_ES::printInfo(ostream& os) const {
	os << "Covariance Matrix Adaptation Evolutionary Strategy" << endl;
	os << "PopSize: " << getPopulationSize() << endl;
}

bool sortByObj0(const Individual* i1, const Individual* i2) {
	return (i1->getInternalOptDirection(0) == MINIMIZE) ?
	       (i1->getObj(0) < i2->getObj(0)) : (i1->getObj(0) > i2->getObj(0));
}

const double CMA_ES::DEF_SIGMA = 2;
const int CMA_ES::DEF_LAMBDA = 500;
const int CMA_ES::MIN_ARGS = 3;
const int CMA_ES::MAX_ARGS = 4;
const int CMA_ES::MAX_IT = 10;