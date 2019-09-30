#include "gsf.h"

VectorXd softThresholding(const VectorXd& z, double lambda){
  double c = 1 - (lambda / z.norm());

  if (c > 0) return c * z;
  else return VectorXd::Zero(D);
}

VectorXd scadUpdate(double u, const Matrix<double, 1, Dynamic>& z, double lambda, double a){
  double normZ = z.norm();

  if(normZ <= (u+1) * lambda){
    return softThresholding(z, u * lambda);

  } else if( (u + 1) * lambda <= normZ && normZ < a * lambda) {
    return ((a - 1)/(a - u - 1)) * softThresholding(z, (a * u * lambda)/(a - 1));

  } else {
    return z;
  }
}

VectorXd mcpUpdate(double u, const Matrix<double, 1, Dynamic>& z, double lambda, double a){
  double normZ = z.norm();

  if(normZ <= a * lambda){
    return (a/(a - u)) * softThresholding(z, u * lambda);

  } else {
    return z;
  }
}

// Note: in this case, "a" is the Adaptive Lasso weight. 
VectorXd adaptiveLassoUpdate(double u, const Matrix<double, 1, Dynamic>& z, double lambda, double a){
  return softThresholding(z, u * lambda * a);
}

VectorXd scadLLAUpdate(double u, const Matrix<double, 1, Dynamic>& z, const Matrix<double, 1, Dynamic>& eta, double lambda, double a){
  double scad, normEta = eta.norm();

  if (normEta <= lambda) {
    scad = lambda;

  } else if (lambda < normEta && normEta < a * lambda) {
    scad = (a * lambda - normEta) / (a - 1);

  } else {
    scad = 0;
  }

  return softThresholding(z, u * scad);
}

VectorXd mcpLLAUpdate(double u, const Matrix<double, 1, Dynamic>& z, const Matrix<double, 1, Dynamic>& eta, double lambda, double a){
  double mcp, normEta = eta.norm();

  if (normEta <= lambda * a) {
    mcp = (lambda - (normEta / a));

  } else {
    mcp = 0;
  }

  return softThresholding(z, u * mcp);
}






