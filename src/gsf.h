#include <Rcpp.h>

#include <RcppEigen.h>
#include <iostream>
#include <float.h>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <sstream>

using namespace Eigen;

extern int K, D, N, n, M; 
extern double upsilon;
extern bool arbSigma;

typedef struct PsiStruct {
  VectorXd pii;
  MatrixXd theta;
  MatrixXd sigma;

  PsiStruct() {
    sigma = MatrixXd::Identity(N, N);
  }

  PsiStruct operator - (const PsiStruct psi) const {
    PsiStruct result;
    result.pii = pii - psi.pii;
    result.theta = theta - psi.theta;
    result.sigma = sigma - psi.sigma;

    return result;
  }

  double distance (PsiStruct psi) {
    return sqrt( (pii - psi.pii).squaredNorm() + (theta - psi.theta).squaredNorm() + (sigma - psi.sigma).squaredNorm() );
  }
} Psi;


/*** Functions for computing the permutation. ***/
bool isOrdered (int, int);
MatrixXd reorder (const MatrixXd&);
void alpha (const MatrixXd&, std::vector<int>&);
double thetaDist (const VectorXd&, const VectorXd&);
MatrixXd getDistanceMatrix (const MatrixXd&) ;
bool find (std::vector<int>, int);
MatrixXd reorderTheta (const MatrixXd& theta);
Psi reorderResult (const Psi&);
double fullLogLikFunction(const MatrixXd& y, const MatrixXd& theta, const VectorXd& pii, const MatrixXd& sigma);
double logLikFunction(const MatrixXd& y, const Psi& psi);

int frequency(const MatrixXd& theta);
VectorXd scadUpdate(double u, const Matrix<double, 1, Dynamic>& z, double lambda, double a);
VectorXd mcpUpdate(double u, const Matrix<double, 1, Dynamic>& z, double lambda, double a);
VectorXd adaptiveLassoUpdate(double u, const Matrix<double, 1, Dynamic>& z, double lambda, double a);
VectorXd scadLLAUpdate(double u, const Matrix<double, 1, Dynamic>& z, const Matrix<double, 1, Dynamic>& eta, double lambda, double a);
VectorXd mcpLLAUpdate(double u, const Matrix<double, 1, Dynamic>& z, const Matrix<double, 1, Dynamic>& eta, double lambda, double a);

/*** Auxiliary functions for Normal mixtures in location. ***/

double densityNormalLoc (const Matrix<double, 1, Dynamic>& y,   
                         const Matrix<double, Dynamic, 1>& theta,   
                         const MatrixXd& sigma);
MatrixXd gradBNormalLoc (const MatrixXd& theta, const MatrixXd& sigma);
double bNormalLoc (const VectorXd& theta, const MatrixXd& sigma);
MatrixXd tNormalLoc (const MatrixXd& y); 
MatrixXd transfNormalLoc (const MatrixXd& theta, const MatrixXd& sigma);
MatrixXd invTransfNormalLoc (const MatrixXd& theta, const MatrixXd& sigma);
bool constrCheckNormalLoc (const MatrixXd& theta);

/*** Auxiliary functions for multinomial mixtures. ***/

double densityMultinomial(const Matrix<double, 1, Dynamic>& y,
                          const Matrix<double, Dynamic, 1>& theta,
                          const MatrixXd& sigma);
MatrixXd gradBMultinomial (const MatrixXd& theta, const MatrixXd& sigma);
double bMultinomial (const VectorXd& theta, const MatrixXd& sigma);
MatrixXd tMultinomial (const MatrixXd& y);
MatrixXd transfMultinomial (const MatrixXd&);
MatrixXd invTransfMultinomial (const MatrixXd&, const MatrixXd&);
bool constrCheckMultinomial (const MatrixXd& theta);


/*** Auxiliary functions for Poisson mixtures. ***/

double densityPoisson(const Matrix<double, 1, Dynamic>& y,
                          const Matrix<double, Dynamic, 1>& theta,
                          const MatrixXd& sigma);
MatrixXd gradBPoisson (const MatrixXd& theta, const MatrixXd& sigma);
double bPoisson (const VectorXd& theta, const MatrixXd& sigma);
MatrixXd tPoisson (const MatrixXd& y);
MatrixXd transfPoisson (const MatrixXd&);
MatrixXd invTransfPoisson (const MatrixXd&, const MatrixXd&);
bool constrCheckPoisson (const MatrixXd& theta);


/*** Auxiliary functions for Exponential mixtures. ***/

double densityExponential(const Matrix<double, 1, Dynamic>& y,
                          const Matrix<double, Dynamic, 1>& theta,
                          const MatrixXd& sigma);
MatrixXd gradBExponential (const MatrixXd& theta, const MatrixXd& sigma);
double bExponential (const VectorXd& theta, const MatrixXd& sigma);
MatrixXd tExponential (const MatrixXd& y);
MatrixXd transfExponential (const MatrixXd&);
MatrixXd invTransfExponential (const MatrixXd&, const MatrixXd&);
bool constrCheckExponential (const MatrixXd& theta);





