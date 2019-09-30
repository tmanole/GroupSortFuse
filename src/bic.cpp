#include "gsf.h"

double (*densityFun)(const Matrix<double, 1, Dynamic>&, const Matrix<double, Dynamic, 1>&, const MatrixXd&);

// [[Rcpp::export(.bicLogLik)]]
extern "C" SEXP bicLogLik (SEXP argY, SEXP argTheta, SEXP argPii, SEXP argSigma, SEXP argIndex) {
  MatrixXd y     = Rcpp::as<MatrixXd>(argY);
  MatrixXd theta = Rcpp::as<MatrixXd>(argTheta); 
  MatrixXd sigma = Rcpp::as<MatrixXd>(argSigma);
  VectorXd pii   = Rcpp::as<VectorXd>(argPii);

	double temp, loglikSum = 0.0;
 
  MatrixXd th = theta;

 
  switch(Rcpp::as<int>(argIndex)) {    
                                       
    // User selected a multivariate normal mixture in location. 
    case 1: densityFun = &densityNormalLoc;
            theta      = transfNormalLoc(theta, sigma);
            break;

    // User selected a multinomial mixture.
    case 3: densityFun = &densityMultinomial;
            theta      = transfMultinomial(theta);
            break;
  
    // User selected a Poisson mixture.
    case 5: densityFun = &densityPoisson;
            theta      = transfPoisson(theta);
            break;

    // User selected a mixture of Exponential distributions.
    case 6: densityFun = &densityExponential;
            theta      = transfExponential(theta);
            break;            
  }

  for (int i = 0; i < y.rows(); i++) {
    temp = 0.0;

    for (int k = 0; k < theta.cols(); k++) {
      temp += pii(k) * (*densityFun)(y.row(i), theta.col(k), sigma);
    }

    loglikSum += log(temp);
  }

  return Rcpp::wrap(loglikSum);
}
