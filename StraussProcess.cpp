#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
int sfun (arma::mat v, double r){
  
  //Get number of elements
  int n = v.n_rows;
  
  //NUmber of elements at distance < r
  int count = 0;
  
  //Loop through all elements
  for (int k = 0; k < n; k++){
    for (int j = (k + 1); j < n; j++){
      if ( arma::norm(v.row(k) - v.row(j)) < r ){
        count = count + 1;
      }
    }
  }
  
  return count;
}

// [[Rcpp::export]]
arma::mat rStraussProcess(double gamma, double r, double beta, int n, int niter, arma::vec xlim, arma::vec ylim) {
  
  //Simulate the number of points on space
  int deleterow;
  bool bflag;
  int spairsSimMatrix;
  int spairsSimAux;
  int pairdif;
  double alpha;
  arma::mat SimMatrix;
  arma::mat SimAux;
  arma::vec aux;
  
  //Verify that there are points here:
  if (n > 0){
    
    //Distribute the number of points uniformly in space
    SimMatrix = arma::mat(n, 2, arma::fill::none);
    
    //Simulate a homogeneous poisson process
    SimMatrix.col(0) = as<arma::vec>(Rcpp::runif(n, xlim(0), xlim(1)));
    SimMatrix.col(1) = as<arma::vec>(Rcpp::runif(n, ylim(0), ylim(1)));
    
    //Loop through all
    for (int i = 0; i < niter; i++){
      
      n = SimMatrix.n_rows;
      
      //Birth process
      bflag = (Rcpp::runif(1,0,1)[0] < 0.5);
      if (bflag){
        
        //Create new point
        aux = arma::vec(2);
        aux(0) = Rcpp::runif(1, xlim(0), xlim(1))[0];
        aux(1) = Rcpp::runif(1, ylim(0), ylim(1))[0];
        
        //Add new point
        SimAux = arma::join_cols(SimMatrix, aux.t());
        
      //Death process  
      } else {
        
        //Cannot kill if none alive
        SimAux = SimMatrix;
        
        //Get row to delete
        if (n > 0){
          deleterow = floor(n*Rcpp::runif(1,0,1)[0]);
          SimAux.shed_row(deleterow);
        }
        
      }
      
      //Get number of pairs of points at distance < r
      spairsSimMatrix = sfun(SimMatrix, r);
      spairsSimAux = sfun(SimAux, r);

      //Check number of points at distance 
      pairdif = spairsSimAux - spairsSimMatrix;
      
      //Check pair difference
      if (pairdif < 0){
        alpha = 1;
      //Birth
      } else if (bflag) {
        alpha = beta*pow(gamma, pairdif)/ ( (double) SimAux.n_rows );
      //Death  
      } else {
        alpha = ( (double) SimMatrix.n_rows)*pow(gamma, pairdif)/beta;
      }
      if (Rcpp::runif(1,0,1)[0] < alpha){
        SimMatrix = SimAux;
      }
    }
    
  }
  return SimMatrix;
}