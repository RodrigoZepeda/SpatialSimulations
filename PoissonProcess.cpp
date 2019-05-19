#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat rPoissonProcess(double muE, Function lambda) {
  
  //Simulate the number of points on space
  int N = Rcpp::rpois(1, muE)[0];
  arma::mat SimMatrix;
  arma::mat pprocess;
  arma::vec pcriteria;
  arma::vec lambdavec;
  arma::uvec ids;
  
  //Verify that there are points here:
  if (N > 0){
    
    //Distribute the number of points uniformly in space
    SimMatrix = arma::mat(N, 2, arma::fill::none);
    
    //Simulate a homogeneous poisson process
    SimMatrix.col(0) = as<arma::vec>(Rcpp::runif(N, 0, 1));
    SimMatrix.col(1) = as<arma::vec>(Rcpp::runif(N, 0, 1));
    
    //Acceptance and rejection criteria
    pcriteria = as<arma::vec>(Rcpp::runif(N,0,1));
    
    //Create normalized lambda vector
    lambdavec = as<arma::vec>(lambda(SimMatrix.col(0), SimMatrix.col(1)))/muE;
    
    //Check 
    ids = arma::find(lambdavec <= pcriteria); // Find indices
    pprocess = SimMatrix.rows(ids);
    
  }
  return pprocess;
}


/*** R
plot(rPoissonProcess(10, function(x,y){rep(1 , length(x))})) #Aproximadamente 10 pts
plot(pp, xlim = c(0,1), ylim = c(0,1))

plot(rPoissonProcess(10, function(x,y){x^2 + y^2})) #Aproximadamente 10 pts
plot(pp, xlim = c(0,1), ylim = c(0,1))

*/
